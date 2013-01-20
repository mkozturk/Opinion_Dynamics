#include <stdio.h>
#include <gsl/gsl_rng.h>

#define MAX(x,y) ((x)>(y) ? (x) : (y) )
#define MIN(x,y) ((x)<(y) ? (x) : (y) )

typedef struct {
	int N;	// Number of agents
	int Q;	// Number of opinions
	int tol;	// Tolerance (confidence interval)
	int *x;		// Array of agent opinions (size N+1, offset 1)
	int *n;		// Array of popularity of opinions (size Q+1, offset 1)
	int *a;		// Array of activity of opinions (size Q+1, offset 1)
	gsl_rng *r;	// GSL random number generator
	void (*biasfun)(int, int, void *, int *, int *); //Pointer to bias func.
	double p;
} opdyn_workspace;

opdyn_workspace * opdyn_workspace_alloc(int N, int Q, int tol);

void opdyn_workspace_free( opdyn_workspace *w);

void update_active_ops( opdyn_workspace *w);

void opdyn_workspace_init(
	opdyn_workspace *w, int *x, double p,
	void (*biasfun)(int, int, void *, int *, int *));
	
int issteady(opdyn_workspace *w);

void Unbiased(int a1, int a2, void *vw, int *strong, int *weak);

void BiasToMajority_pair(int a1, int a2, void *vw, int *strong, int *weak);

void BiasToMajority_global(int a1, int a2, void *vw, int *strong, int *weak);

int interact( opdyn_workspace *w);

void disp_nops( FILE *stream, opdyn_workspace *w);

int run( opdyn_workspace *w, int maxsteps, int display);

double get_mean_op( opdyn_workspace *w);
