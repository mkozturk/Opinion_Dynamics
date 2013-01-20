function steady = issteady(nx, bci)
% Checks whether the distribution of opinions nx is a steady state.
% bci is the bounded confidence interval.
% nx is an array of length M (number of opinions).
% nx(k) is the number of followers of the opinion k.

% Requires that all agents have the same bci.

	k=1; steady=1; M=size(nx,2);
	% Move to the first opinion with nonzero population.
	while nx(k)==0, k=k+1; end
	i1=k; k=k+1;

	while 1
		% Move to the next opinion with nonzero population.
		while k<=M && nx(k)==0, k=k+1;end
		if k>M, break, end
		i2=k; k=k+1;
		% If the opinion difference is less than bci
		% steady state is not achieved.
		if i2-i1<=bci
			steady=0;
			break
		end
		i1=i2;
	end
end