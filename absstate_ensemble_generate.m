% Runs the C program "ensemble" for various N,Q,d and bias.
% For each run:
%	- Generates data file using the program.
%	- Prints the probability vector of having 1,2,3,... different opinions
%		at absorbing state.
%	- Plots the cumulative distribution function for absorption time.

ensemblesize = 1e5;
N0 = 100; % the integer divisible by Q and <=N0 will be used.
for bias=0:1
	for Q=3:7
		for d=1:Q-2
			N = N0 - mod(N0,Q); % make N divisible by Q
			outfile = sprintf('N%d_Q%d_d%d_ens%.0e_bias%d.dat',...
				N, Q, d, ensemblesize, bias); % file to store output
			execfile = './ensemble';
			command = sprintf('%s %d %d %d %d %d > %s',...
				execfile, N, Q, d, ensemblesize, bias, outfile);
			tic
			system(command);
			toc
			disp(sprintf('N=%d, Q=%d, d=%d ensemble size= %.0e, bias=%d',...
				N, Q, d, ensemblesize, bias));
			M=dlmread(outfile);
			
			% Display the distribution of the number of states
			nstates=zeros(Q,1);
			for i=1:ensemblesize
				nstates(i)=sum(M(i, 2:end)>0);
			end
			ns = hist(nstates, 1:5);
			disp('Frequency of the number of opinions in absorbing state:');
			disp(ns/ensemblesize);
			

		end
	end
end
