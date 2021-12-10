function W = initialnet_gamma(N, p, gamma)

%parameters from Guillaumes paper

R =10;
M=N/2;

%Create random matrix with zeros on the diagonal and density p elsewhere
NN = round(p*N*(N-1));
fill = [ones(1,NN), zeros(1,N*(N-1) - NN)];
fill = reshape(fill(randperm(N*(N-1))),N,N-1);

W1 = zeros(N);
W1(1:end-1,2:end) = fill(1:end-1,:);

W2 = zeros(N);
W2(2:end,1:end-1) = fill(2:end,:);

W = triu(W1,1) + tril(W2,-1);

%Give some structure to the E-I connections
%MM = round(p/2*M*M);
%fill = [ones(1,MM), zeros(1,M^2 - MM)];
%fill = reshape(fill(randperm(M*(M))),M,M);
%W(M+1:N, 1:M)=fill;


%Create synaptic strengths
w0 = sqrt(2)*R/(sqrt(p*(1-p)*(1 + gamma^2)));

W = W * (w0/sqrt(N)); %E synapses
W(:,(N/2 + 1):end) = -gamma*W(:,(N/2 +1):end); %I synapses

