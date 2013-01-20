% Markov matrix for opinion dynamics with Q=3, d=1
% Using only possible states in the probability vector.
Q=3;
N=30;
sz = (N+1)*(N+2)/2; %  markov matrix size
A=zeros(sz);
% r(n1,n2,n3) gives the 0-based index corresponding to state (n1,n2,n3).
r = @(n1,n2, n3) (N+2)*n1 - n1*(n1+1)/2 + n2;
% Create lookup table for states
states = zeros(sz, 3);

for n1=0:N
	for n2=0:(N-n1)
		n3 = N-n1-n2;
		states( r(n1, n2, n3) + 1, : ) = [n1,n2,n3]; 
		for j=0:(sz-1)
			a = 0;
				if n2>1 && r(n1+1, n2-1, n3)==j
					a = a + (n1+1)*(n2-1);
				end
				if n3>1 && r(n1, n2+1, n3-1)==j
					a = a + (n2+1)*(n3-1);
				end
				if n1>1 && r(n1-1, n2+1, n3)==j
					a = a + (n1-1)*(n2+1);
				end
				if n2>1 && r(n1, n2-1, n3+1)==j
					a = a + (n2-1)*(n3+1);
				end
				if r(n1,n2,n3)==j
					a = a + 2*(N*(N-1)-n1*n2-n2*n3);
				end
				A( r(n1, n2, n3)+1, j+1 ) = a / (2*N*(N-1));
		end
	end
end

[V, D] = eig(A);

p0 = zeros(sz, 1);
p0( r(N/3, N/3, N/3) + 1) = 1; % Initial state: n1=n2=n3=N/3

ind = find(diag(D)==1); % indices of 1 eigenvalues (1-based)

c = inv(V)*p0; % initial state wrt eigenvector basis
% print the end states (eigenvectors with eigenvalue 1)
% and the probability of that end state.
disp('End states:');
for i=ind'
	i1 = find(V(:,i)==1);
	disp([states(i1,:), c(i)]);
end

% Plot the probability T(t) that system stops at time t.
a = zeros((N+1)*(N+2)/2, 1)';
for i = ind'
	a( V(:,i)==1 ) = 1;
end
tmax = 1000;
p=zeros((N+1)*(N+2)/2, tmax+1);
p(:,1) = inv(V)*p0; % wrt eigenvector basis
for t=2:tmax+1
	p(:,t) = diag(D).*p(:,t-1);
end
p = V*p; % convert back to original basis
plot(0:tmax, a*p);