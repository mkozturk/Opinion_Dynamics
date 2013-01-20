/* Agent-based opinion dynamics program.
Exactly one interaction per time step.
Fully connected homogeneous agents.

Prints the state (opinions) at each step.

Compile with:
gcc singlerun.c libopdyn.c -std=gnu99 -Wall -lm -lgsl -lgslcblas -o singlerun

Usage:
singlerun <number of agents> <number of opinions> <confidence bound> [biastype]

Output: Popularity vector at each step.
*/

#include <stdio.h>
#include <stdlib.h>
#include "opdyn.h"

int main( int argc, char *argv[])
{
	double bias = 0.5;
	int biastype = 0;

	int c; // counter
	
	if(argc<4){
		fprintf(stderr,
		"Usage: singlerun N Q tol [biastype]\n"
		"biastype 1: p = 0.5 + 0.5*|n-n'|/(n+n')\n"
		"biastype 2: p = 0.5 + 0.5*|n-n'|/N\n");
		exit(1);
	}
	
	int N = atoi(argv[1]);
	int Q = atoi(argv[2]);
	int tol = atoi(argv[3]);
	if(argc>=5)
		biastype = atoi(argv[4]);
	
	opdyn_workspace *w = opdyn_workspace_alloc(N, Q, tol);

	int x[N+1];
	int endtime;

	/* distribute opinions uniformly on agents */
	for(c=1; c<=N; c++)
		x[c] = (c % Q) + 1;
	/* Initialize workspace struct */
	if(biastype==0)
		opdyn_workspace_init(w, x, bias, Unbiased);
	if(biastype==1)
		opdyn_workspace_init(w, x, bias, BiasToMajority_pair);
	if(biastype==2)
		opdyn_workspace_init(w, x, bias, BiasToMajority_global);
						
	endtime = run(w, 0, 1);
		
	return 0;
}
