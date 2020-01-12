networkparams.N = 2^3;
networkparams.chainlen = 2.^(0:log2(networkparams.N));
networkparams.gamma = 0.9;
networkparams.beta = [0 1];
networkparams.alpha = 1:10;

kmax = 50;
signalparams.density = (1:2:networkparams.N)/networkparams.N; % signal density
signalparams.rho = 0:0.25:1; % input signal correlation
signalparams.loc = 1;
signalparams.var = 1; % input signal variance
noiseparams.var = 1; % noise variance

%% construct combinations of parameters
[gamma,beta,alpha,N,chainlen,density,rho,loc] = ...
    ndgrid(networkparams.gamma,networkparams.beta,networkparams.alpha,...
    networkparams.N,networkparams.chainlen,...
    signalparams.density,signalparams.rho,signalparams.loc);
gamma = gamma(:);
beta = beta(:);
alpha = alpha(:);
N = N(:);
chainlen = chainlen(:);
density = density(:);
rho = rho(:);
loc = loc(:);

%%
duplicates = (beta == 0 & alpha > 1);
gamma(duplicates) = []; beta(duplicates) = []; alpha(duplicates) = [];
N(duplicates) = []; chainlen(duplicates) = [];
density(duplicates) = []; rho(duplicates) = []; loc(duplicates) = [];
numcases = numel(gamma);

%% simulate each combination
for k=1:numcases
    W = DesignNonNormal(N(k),chainlen(k),gamma(k),beta(k),alpha(k));
    [inputstats,wts_decode,var_decode,mse_decode,corr_decode] = ...
        DecodeMultivariateSignal(W,density(k),rho(k),loc(k),signalparams.var,noiseparams.var,kmax);
end

%% function to design connectivity
function W = DesignNonNormal(N,chainlen,gamma,beta,alpha)
    W = gamma*eye(N) + beta*circshift(eye(N),1) + alpha*circshift(eye(N),-1);
    W(1,end) = 0; W(end,1) = 0; % no end-to-end connections
    for i=chainlen:chainlen:N-1, W(i,i+1) = 0; W(i+1,i) = 0; end
end

%% function to decode input pattern
function [inputstats,wts_decode,var_decode,mse_decode,corr_decode] = ...
        DecodeMultivariateSignal(W,density,rho,loc,var_signal,var_noise,kmax)
    
    N = size(W,1);
    %% generate input signal
    tmax = 3e5;
    Sig = zeros(N);
    Sig(1:(density*N),1:(density*N)) = rho;
    Sig = Sig - diag(diag(Sig)) + double(diag(1:N <= (density*N)));
    Sig = Sig*var_signal;
    s = mvnrnd(zeros(1,N),Sig,tmax);
    
    %% generate noise
    z = normrnd(0,sqrt(var_noise),tmax,N);  
    
    %% simulate the dynamics 
    X = zeros(tmax,N);
    for t = 2:tmax, X(t,:) = X(t-1,:)*W + s(t,:) + z(t,:); end
end