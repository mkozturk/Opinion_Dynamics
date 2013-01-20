% Plots the probability of absorbing state [i, 0, N-i] vs. i
% Runs the C program "ensemble"

N = 30;	% number of agents
Q = 3;	% range of opinions, from 1 to Q
d = 1; % bounded confidence interval
ensemblesize = 1e4;
bias = 0; % 0 for unbiased, 1 for bias to pairwise majority.

% If "ensemble" is not in current directory, write the correct path below.
path = './';
command = sprintf('%sensemble %d %d %d %d %d > temp',...
	path,N, Q, d, ensemblesize, bias);
tic
system(command);
toc
M=dlmread('temp');
system('rm temp');
c = zeros(N+1,1); % c(i+1) is the frequency of state [i, 0, N-i]
for e=1:ensemblesize
	s = M(e,2:end);
	if s(2)==0
		c(s(1)+1) = c(s(1)+1) + 1;
	end
end
c = c/ensemblesize;
% Plot c vs i
figure(); plot(linspace(0,1,N+1), c, '.-');
xlabel('i/N'); ylabel('Pr{ X_s = (i, 0, N-i) }');
