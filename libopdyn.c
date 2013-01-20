// select agents, not opinions, during interaction.
#include <stdio.h>
#include <stdlib.h>
#include <time.h> // to set the random seed
#include <math.h>	// for ceil()
#include "opdyn.h"

opdyn_workspace * opdyn_workspace_alloc(int N, int Q, int tol){
	opdyn_workspace *w = malloc(sizeof(opdyn_workspace));
	w->N = N;
	w->Q = Q;
	w->tol = tol;
	w->x = malloc((N+1)*sizeof(int));
	w->n = malloc((Q+1)*sizeof(int));
	w->a = malloc((Q+1)*sizeof(int));
	w->r = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(w->r, time(0));
	return w;
}

void opdyn_workspace_free( opdyn_workspace *w){
	free(w->n);
	free(w->a);
	free(w->x);
	free(w);
}

void opdyn_workspace_init(opdyn_workspace *w, int *x, double p, void (*biasfun)(int, int, void *, int *, int *))
{
	int c; //counter
	w->x = x;
	
	/* Initialize the histogram of opinions */
	for (c=1; c<= w->Q; c++)
		w->n[c]=0;	
	for (c=1; c <= w->N; c++)
		w->n[ w->x[c] ]++;

	update_active_ops(w);
	w->biasfun = biasfun;
	w->p = p;
}

void update_active_ops( opdyn_workspace *w)
{
	int i,j; // counters
	int activity;
	
	for (i=1; i<=w->Q; i++){
		activity=0;
		if (w->n[i]!=0){
			for (j=MAX(1,i-w->tol); j<=MIN(w->Q,i+w->tol); j++){
				if(j==i) continue;
				if(w->n[j]!=0){
					activity=1;
					break;	// out of j loop
				}
			}
		}
		w->a[i] = activity;
	}
}

int issteady(opdyn_workspace *w)
{
	int i; // counter
	int steady=1;
	for(i=1; i <= w->Q; i++){
		if(w->a[i]==1){
			steady=0;
			break;
		}
	}
	return steady;
}

// BIAS FUNCTIONS
void Unbiased(int a1, int a2, void *vw, int *winner, int *loser)
/* Equally likely. */
{
	/* we have to take the workspace pointer as pointer to void
		and cast to opdyn_workspace because we cannot define the
		fn. prototype in struct with opdyn_workspace *w */
	opdyn_workspace * w = vw;
	
	if(gsl_rng_uniform(w->r) < 0.5){
		*winner = a1;
		*loser = a2;
	}
	else{
		*winner = a2;
		*loser = a1;
	}
	
}

void BiasToMajority_pair(int a1, int a2, void *vw, int *winner, int *loser)
/* Pairwise majority bias. When agents a1, a2 with opinions q1,q2 meet,
 a1 wins with probability n1/(n1+n2) */

{
	opdyn_workspace * w = vw;
	double p;
	int q1 = w->x[a1], q2 = w->x[a2];
	int n1 = w->n[q1], n2 = w->n[q2];
	
	p = 1.0 * n1 / (n1 + n2); // probability that a1 is the winner
	
	if ( gsl_rng_uniform(w->r) < p ){
		*winner = a1; *loser = a2;}
	else {
		*winner = a2; *loser = a1;}
}

void BiasToMajority_global(int a1, int a2, void *vw, int *winner, int *loser)
/* Global majority bias. When agents a1, a2 with opinions q1 and q2 meet,
  the one with the higher popularity n wins with probability n/N */
{
	opdyn_workspace * w = vw;
	double p;
	int q1=w->x[a1], q2 = w->x[a2];
	int n1=w->n[q1], n2 = w->n[q2];
	
	p = 0.5 + fabs(0.5*(n1 - n2)/ w->N);	
	
	if ( n1 == n2 )
		Unbiased(a1, a2, vw, winner, loser);
		
	if ( n1 > n2 ){
		if ( gsl_rng_uniform(w->r) < p ){
			*winner = a1; *loser = a2;}
		else {
			*winner = a2; *loser = a1;}
	}
	
	if ( n1 < n2 ){
		if ( gsl_rng_uniform(w->r) < p ){
			*winner = a2; *loser = a1;}
		else {
			*winner = a1; *loser = a2;}
	}
}

/* END OF BIAS FUNCTIONS */

int interact( opdyn_workspace *w)
{
	int a1, a2;	// Interacting agents
	int q1, q2; // Opinions of a1, a2, resp.
	int winner, loser; // agents
	int retval = 0; // return value
	
	// select one agent randomly. Make sure its opinion is active.
	do {
		a1 = gsl_rng_uniform_int(w->r, w->N) + 1;
		q1 = w->x[a1];
	}
	while(w->a[q1]==0);
			
	// Select a compatible agent with opinion between max(1, q1-d) and min(Q, q1+d).
	do{
		a2 = gsl_rng_uniform_int( w->r, w->N) +1;
		q2 = w->x[a2];
	} while (a2==a1 || q2==q1 || q2 > MIN(w->Q, q1+w->tol) || q2 < MAX(1,q1-w->tol));
	
	(w->biasfun)(a1, a2, w, &winner, &loser);

	w->n[w->x[winner]]++;
	w->n[w->x[loser]]--;
	
	// If the popularity of the losing opinion has decreased to zero,
	// update the activity vector and check if steady state is achieved.
	if (w->n[w->x[loser]]==0){
		update_active_ops(w);
		if(issteady(w)){
			retval = 1;
		}
	}
	
	// Update agent opinions last.
	w->x[loser] = w->x[winner];
	
	return retval; /* 0 for normal, 1 for steady state */
}

void disp_nops( FILE *stream, opdyn_workspace *w){
	for(int c=1; c <= (w->Q - 1) ; c++)
		fprintf(stream, "%d ", w->n[c]);
	fprintf(stream,"%d\n", w->n[w->Q]);
}

int run( opdyn_workspace *w, int maxsteps, int display)
{
	int s;	// counters
	int steadyflag=0;	// becomes 1 when steady state is achieved.	
	// Time iteration loop:
	for(s=1; s<=maxsteps || maxsteps==0; s++){
		if(display){
			disp_nops(stdout, w);
		}				
		if(steadyflag) break;		
		steadyflag=interact(w);
	}
	return s;
}

double get_mean_op( opdyn_workspace *w)
{
	double sum=0.0;
	for(int q=1; q <= w->Q; q++)
		sum += q * w->n[q];
	return sum/w->N;
}
