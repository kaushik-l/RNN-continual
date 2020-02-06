N = 100;
gamma = 0.005:0.005:1; Ng = numel(gamma);
beta = 0.005:0.005:1; Nb = numel(beta);
alpha = 1; Na = numel(alpha);
for i=1:Ng
    for j=1:Nb
        W = DesignNonNormal(N,N,gamma(i),beta(j),alpha);
        [~,D] = eig(W);
        lambdas(i,j,:) = sort(diag(D));
    end
end

%% plot
figure; imagesc(beta,gamma,squeeze(max(lambdas,[],3))<1);

%%
function W = DesignNonNormal(N,chainlen,gamma,beta,alpha)
    W = gamma*eye(N) + (beta*alpha)*circshift(eye(N),-1) + (beta/alpha)*circshift(eye(N),1);
    W(1,end) = 0; W(end,1) = 0; % no end-to-end connections
    for i=chainlen:chainlen:N-1, W(i,i+1) = 0; W(i+1,i) = 0; end
end