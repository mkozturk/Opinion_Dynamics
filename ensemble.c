/* Agent-based opinion dynamics program.
Exactly one interaction per time step.
Fully connected homogeneous agents.

Makes an ensemble of runs. For each run, prints the absorption time and absorbing state.

Compile with:
gcc ensemble.c libopdyn.c -std=gnu99 -Wall -lm -lgsl -lgslcblas -o ensemble

Usage:
ensemble <number of agents> <number of opinions> <confidence bound> <number of runs> [biastype]

Output: Absorption time and final popularity vector of opinions (absorbing state) for each run.
*/

#include <stdio.h>
#include <stdlib.h>
#include "opdyn.h"

int main( int argc, char *argv[])
{
	double bias = 0.5;
	int biastype = 0;

	int c; // counter
	
	if(argc<5){
		fprintf(stderr,
		"Usage: ensemble N Q tol ensemblesize [biastype]\n"
		"biastype 1: p = 0.5 + 0.5*|n-n'|/(n+n')\n"
		"biastype 2: p = 0.5 + 0.5*|n-n'|/N\n");
		exit(1);
	}
	
	int N = atoi(argv[1]);
	int Q = atoi(argv[2]);
	int tol = atoi(argv[3]);
	int enssize = atoi(argv[4]);
	if(argc>=6)
		biastype = atoi(argv[5]);
	
	opdyn_workspace *w = opdyn_workspace_alloc(N, Q, tol);

	int x[N+1];
	int endtime;
	for (int k=1; k<=enssize; k++){
		/* distribute opinions uniformly on agents */
		for(c=1; c<=N; c++)
			x[c] = (c % Q) + 1;
		// Initialize
		if(biastype==0)
			opdyn_workspace_init(w, x, bias, Unbiased);
		if(biastype==1)
			opdyn_workspace_init(w, x, bias, BiasToMajority_pair);
		if(biastype==2)
			opdyn_workspace_init(w, x, bias, BiasToMajority_global);
						
		endtime = run(w, 0, 0);
		
		printf("%d ", endtime);
		disp_nops(stdout, w);
	}
	return 0;
}
