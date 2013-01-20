function [x nx] = BiasToMean(i,j,x,nx,p, param)
% Pairwise interaction with bias.
% The opinion which is closer to the overall average opinion
% is accepted with probability p.
% i,j : Indices of interacting agents.
% x : array of opinions (size = number of agents)
% nx : array of the population of opinions (size = number of opinions)
% param required for compatibility, ignored.

% Output: Updated x and nx.

avg = sum(x)/length(x);

% stronger: The agent that is favored with probability p > 0.5.
% weaker: The agent that is not favored, probability 1-p
if abs(x(i)-avg) < abs(x(j)-avg)
	% i is closer to the mean than j
	stronger = i; weaker = j;
else
	% j is closer to mean than i
	stronger = j; weaker = i;
end

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