% Visualize dynamics of opinions with given parameters

% Initialize
N = 20;	% number of agents
Q = 5;	% range of opinions, from 1 to Q
bci = 1; % bounded confidence interval
p = 0.5; % Probability of opinion change for biased interactions.
interaction = @Unbiased;
% Can be @BiasToMean, @BiasToMajority, @Unbiased
displaystep = 1; % time interval between displays
maxstep = 1e5;
tracer = N/2; % choose this agent as tracer on display

%x = ceil(Q.*rand(N,1)); % assign opinions randomly from 1 to 100.
x = mod(1:N, Q)+1; % assign opinions uniformly.
nx = hist(x, 1:Q); % histogram of opinion values.
s = 0;	% step counter
steady = 0;
%rand('twister',231); % use same random sequence, for testing.

% Run simulation until steady state is achieved or maxsteps is reached.
while s <= maxstep && ~steady

	% Display the data
	if s==0 || mod(s, displaystep)==0		
		% show opinions
		subplot(2,1,1)
		plot(1:N, x,'.', tracer,x(tracer),'ro','LineWidth',2,'MarkerSize',10);
		xlabel('agent'); ylabel('opinion');
		ylim([0,Q+1]);
		%line([0 N],mean(x)*[1 1])
		title(sprintf('confidence interval = %d',bci));
		grid on;

		% show the distribution of opinion values
		subplot(2,1,2);
		bar(1:Q, nx);
		xlabel('opinion'); ylabel('# of followers');
		xlim([0,Q+1]); ylim([0,N+1]);
		text(0.8*Q,0.8*N,sprintf('step %d',s));
		grid on;

		drawnow;
	end
	 
	% Pick two distinct agents
	i = ceil(N*rand());
	j = ceil(N*rand());
	c=1;
	while ~steady && (j==i || x(i)==x(j) || abs(x(i)-x(j)) > bci)
		i = ceil(N*rand());
		j = ceil(N*rand());
		c = c+1;
		% If too many tries do not give a compatible pair, check for steady
		% state.
		if c>100
			steady = issteady(nx, bci);
			c=0;
		end
	end
	
	if ~steady
		% Interact:
		[x nx] = interaction(i,j,x,nx,p);
	end
	s=s+1;
end
if steady
	disp(sprintf('Absorbing state %s achieved in %d steps.',...
		mat2str(nx), s));
else
	sprintf('Steady state not achieved after %e steps.',maxstep)
end