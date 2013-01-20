% Markov matrix for opinion dynamics with Q=3, d=1
% Enforcing exactly one conversion per step.
tic
Q=3;
N=30; % Takes 700 seconds for N=99
sz = (N+1)*(N+2)/2; %  markov matrix size
A=zeros(sz);
% M(n1,n2,n3) gives the total number of all compatible pairs.
M = @(n1,n2,n3) n1*n2 + n2*n3;
% r(n1,n2,n3) gives the 1-based index corresponding to state (n1,n2,n3).
r = @(n1,n2, n3) (N+2)*n1 - n1*(n1+1)/2 + n2 + 1;
% Initialize lookup table for states
states = zeros(sz, Q);

for n1=0:N
	for n2=0:(N-n1)
		n3 = N-n1-n2;
		i = r(n1,n2,n3);
		states( r(n1, n2, n3), : ) = [n1,n2,n3];
		if n1>1
			A(i, r(n1-1, n2+1, n3)) = 0.5*(n1-1)*(n2+1) / M(n1-1,n2+1,n3);
		end
		if n2>1
			A(i, r(n1+1, n2-1, n3)) = 0.5*(n1+1)*(n2-1) / M(n1+1, n2-1, n3);
			A(i, r(n1, n2-1, n3+1)) = 0.5*(n2-1)*(n3+1) / M(n1, n2-1, n3+1);
		end
		if n3>1
			A(i, r(n1, n2+1, n3-1)) = 0.5*(n2+1)*(n3-1) / M(n1,n2+1,n3-1);
		end
		if M(n1,n2,n3)==0
			A(i, i) = 1;
		end
	end
end

[V, D] = eig(A);

p0 = zeros(sz, 1);
p0( r(N/3, N/3, N/3) ) = 1; % Initial state: n1=n2=n3=N/3

ind = find(diag(D)==1); % indices of 1 eigenvalues (1-based)

c = V\p0; % initial state wrt eigenvector basis (c=inv(V)*p0)
% print the end states (eigenvectors with eigenvalue 1)
% and the probability of that end state.
disp('Absorbing states and their probabilities');
for i=ind'
	i1 = find(V(:,i)==1);
	disp(sprintf('Pr{%s} = %f', mat2str(states(i1,:)), c(i)));
end

% Plot the probability of state (n,0,N-n) vs. n
ps = zeros(N+1,1);
for i=ind'
	i1 = find(V(:,i)==1);
	if(states(i1,2)==0)
		ps(states(i1,1)+1) = c(i);
	end
end
figure();
plot(linspace(0,1,N+1), ps, 'r.-')
title('Probability of state (n,0,N-n)');
xlabel('i/N'); ylabel('Pr[ X_s = (i, 0, N-i) ]');

% Plot the probability T(t) that system stops at time t.
a = zeros(sz, 1)';
for i = ind'
	a( V(:,i)==1 ) = 1;
end
tmax = 1500;
p=zeros(sz, tmax+1);
p(:,1) = inv(V)*p0; % wrt eigenvector basis
for t=2:tmax+1
	p(:,t) = diag(D).*p(:,t-1);
end
p = V*p; % convert back to original basis
figure();
plot(0:tmax, a*p);
title('Cumulative Probability Distribution of absorption time');
xlabel('time t'); ylabel('Pr[ t_{abs} < t ]');
toc