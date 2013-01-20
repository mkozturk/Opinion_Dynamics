function [x nx] = BiasToMajority(i,j,x,nx,p,param)
% Pairwise interaction with bias.
% The opinion with more followers is accepted with probability p.
% i,j : Indices of interacting agents.
% x : array of opinions (size = number of agents)
% nx : array of the population of opinions (size = number of opinions)
% param required for compatibility; ignored.

% Output: Updated x and nx.

% stronger: The agent that is favored with probability p > 0.5.
% weaker: The agent that is not favored, probability 1-p
if nx(x(i)) > nx(x(j))
	stronger = i; weaker = j;
else
	stronger = j; weaker = i;
end

p = 0.5 + 0.5*abs(nx(x(i))-nx(x(j)))/(nx(x(i))+nx(x(j)));
if rand()<p
	% update the histogram
	nx(x(weaker)) = nx(x(weaker)) - 1;
	nx(x(stronger)) = nx(x(stronger)) + 1;
	x(weaker)=x(stronger);
else
	nx(x(weaker)) = nx(x(weaker)) + 1;
	nx(x(stronger)) = nx(x(stronger)) - 1;
	x(stronger)=x(weaker);
end
end