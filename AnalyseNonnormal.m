networkparams.N = 2^3;
networkparams.chainlen = 1:networkparams.N; %2.^(0:log2(networkparams.N));
networkparams.gamma = 0.9;
networkparams.beta = [0 0.05];
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

%% flatten
gamma = gamma(:);
beta = beta(:);
alpha = alpha(:);
N = N(:);
chainlen = chainlen(:);
density = density(:);
rho = rho(:);
loc = loc(:);

%% keep only alpha=1 for beta=0
duplicates = (beta == 0 & alpha > 1);
gamma(duplicates) = []; beta(duplicates) = []; alpha(duplicates) = [];
N(duplicates) = []; chainlen(duplicates) = [];
density(duplicates) = []; rho(duplicates) = []; loc(duplicates) = [];

%% keep only beta=0 for chainlen=1
duplicates = (chainlen == 1 & beta > 0);
gamma(duplicates) = []; beta(duplicates) = []; alpha(duplicates) = [];
N(duplicates) = []; chainlen(duplicates) = [];
density(duplicates) = []; rho(duplicates) = []; loc(duplicates) = [];

%% count cases
numcases = numel(gamma);

%% simulate each combination
wts_decoder = cell(numcases,1); 
var_decoder = cell(numcases,1);
r2_decoder = cell(numcases,1); 
corr_decoder = cell(numcases,1);
wts_pcadecoder = cell(numcases,1); 
var_pcadecoder = cell(numcases,1);
r2_pcadecoder = cell(numcases,1); 
corr_pcadecoder = cell(numcases,1);

parfor k=1:numcases
    fprintf(['simulating case #' num2str(k) ' of ' num2str(numcases) '\n'])
    W = DesignNonNormal(N(k),chainlen(k),gamma(k),beta(k),alpha(k));
    [wts_decoder{k},var_decoder{k},r2_decoder{k},corr_decoder{k},...
        wts_pcadecoder{k},var_pcadecoder{k},r2_pcadecoder{k},corr_pcadecoder{k}] = ...
        DecodeMultivariateSignal(W,density(k),rho(k),loc(k),signalparams.var,noiseparams.var,kmax);
end

%% save
params.gamma = gamma;
params.beta = beta;
params.alpha = alpha;
params.chainlen = chainlen;
params.N = N;
params.density = density;
params.rho = rho;
params.loc = loc;
decoder.wts = wts_decoder;
decoder.var = var_decoder;
decoder.r2 = r2_decoder;
decoder.corr = corr_decoder;
pcdecoder.wts = wts_pcadecoder;
pcdecoder.var = var_pcadecoder;
pcdecoder.r2 = r2_pcadecoder;
pcdecoder.corr = corr_pcadecoder;
signalparams.maxlag = kmax;
save('AnalyseNonnormal.mat','params','signalparams','noiseparams','decoder','pcdecoder');

%% plot results
plotresults = false;
if plotresults
    decoder.corr{i}
end

%% function to design connectivity
function W = DesignNonNormal(N,chainlen,gamma,beta,alpha)
    W = gamma*eye(N) + (beta*alpha)*circshift(eye(N),-1) + (beta/alpha)*circshift(eye(N),1);
    W(1,end) = 0; W(end,1) = 0; % no end-to-end connections
    for i=chainlen:chainlen:N-1, W(i,i+1) = 0; W(i+1,i) = 0; end
end

%% function to decode input pattern
function [wts_decoder,var_decoder,r2_decoder,corr_decoder,...
    wts_pcadecoder,var_pcadecoder,r2_pcadecoder,corr_pcadecoder] = ...
        DecodeMultivariateSignal(W,density,rho,loc,var_signal,var_noise,kmax)
    
    N = size(W,1);
    numfeatures = N*density;
    
    %% generate input signal
    tmax = 3e5;
    Sig = zeros(N);
    Sig(1:numfeatures,1:numfeatures) = rho;
    Sig = Sig - diag(diag(Sig)) + double(diag(1:N <= numfeatures));
    Sig = Sig*var_signal;
    s = mvnrnd(zeros(1,N),Sig,tmax);
    
    %% specify source node
    s = circshift(s,loc-1,2);
    
    %% generate noise
    z = normrnd(0,sqrt(var_noise),tmax,N);  
    
    %% simulate the dynamics 
    X = zeros(tmax,N);
    for t = 2:tmax, X(t,:) = X(t-1,:)*W + s(t,:) + z(t,:); end
    
    %% decode features
    wts_decoder = nan(kmax+1,N,numfeatures);
    var_decoder = nan(kmax+1,numfeatures);
    r2_decoder = nan(kmax+1,numfeatures);
    corr_decoder = nan(kmax+1,numfeatures);
    
    minsampleindx = round(0.1*tmax); % let the response statistics stabilize
    samplesize = 1e5; % number of observations used for regression
    for m=1:numfeatures
        inputfeature = loc + m - 1;
        for k=0:kmax
            trainsamples = randsample(minsampleindx:tmax,samplesize);
            s_train = s(trainsamples-k,inputfeature);
            X_train = X(trainsamples,:);
            wts = regress(s_train,X_train); % obtain readout weights by regressing
            
            wts_decoder(k+1,:,m) = wts;
            var_decoder(k+1,m) = var(X_train*wts);
            r2_decoder(k+1,m) = 1 - (mean((s_train - X_train*wts).^2)/var(s_train));
            corr_decoder(k+1,m) = corr(s_train , X_train*wts);
        end
    end
    
    %% decode PCs
    [~,s_pca,lambda] = pca(s);
    numpcs = sum(lambda ~= 0);
    wts_pcadecoder = nan(kmax+1,N,numpcs);
    var_pcadecoder = nan(kmax+1,numpcs,numpcs);
    r2_pcadecoder = nan(kmax+1,numpcs);
    corr_pcadecoder = nan(kmax+1,numpcs);
    
    minsampleindx = round(0.1*tmax); % let the response statistics stabilize
    samplesize = 1e5; % number of observations used for regression
    for m=1:numpcs
        pcafeature = m;
        for k=0:kmax
            trainsamples = randsample(minsampleindx:tmax,samplesize);
            s_train = s_pca(trainsamples-k,pcafeature);
            X_train = X(trainsamples,:);
            wts = regress(s_train,X_train); % obtain readout weights by regressing
            
            wts_pcadecoder(k+1,:,m) = wts;
            var_pcadecoder(k+1,m) = var(X_train*wts);
            r2_pcadecoder(k+1,m) = 1 - (mean((s_train - X_train*wts).^2)/var(s_train));
            corr_pcadecoder(k+1,m) = corr(s_train , X_train*wts);
        end
    end            
end