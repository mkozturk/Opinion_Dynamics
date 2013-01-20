function [x nx] = Unbiased(i, j, x, nx, p, params)
% Pairwise interaction without bias.
% The opinion of either agent can be accepted with equal probability.
% i,j : Indices of interacting agents.
% x : array of opinions (size = number of agents)
% p and params required for compatibility, it is ignored.

% Output: Updated x.

if rand()<0.5
	% Update the histogram
	nx(x(i)) = nx(x(i))-1;
	nx(x(j)) = nx(x(j))+1;
	% Convert i to j's opinion.
	x(i)=x(j);
else
	% Update the histogram
	nx(x(i)) = nx(x(i))+1;
	nx(x(j)) = nx(x(j))-1;
	% Convert j to i's opinion.
	x(j)=x(i);
end
