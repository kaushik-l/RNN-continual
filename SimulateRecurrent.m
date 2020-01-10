function [inputstats,wts_decode,var_decode,mse_decode,corr_decode,...
    var_optimal,mse_optimal,corr_optimal] = ...
    SimulateRecurrent(v,W,tau_signal,var_signal,tau_noise,var_noise,kmax,wts_optimal)

N = numel(v); % number of units

%% build filter
sz = 2*kmax; %filter size
t2 = linspace(-sz/2, sz/2, sz+1);
h = @(tau) exp(-t2.^2/(2*tau^2))/sum(exp(-t2.^2/(2*tau^2)));

%% simulate the dynamics
tmax = 3e5;
s = normrnd(0,sqrt(var_signal),tmax,1);
if tau_signal > 0, s = conv(s,h(tau_signal),'same'); s = s*sqrt(var_signal)/sqrt(var(s)); end
z = normrnd(0,sqrt(var_noise),N,tmax);
if tau_noise > 0, for i=1:N, z(i,:) = conv(z(i,:),h(tau_noise),'same'); end; end
[inputstats.temporalcorr_signal, inputstats.lags] = xcorr(s(:),kmax,'coeff');
inputstats.temporalcorr_noise = xcorr(z(:),kmax,'coeff');
X = zeros(N,tmax);
for t = 2:tmax, X(:,t) = W*X(:,t-1) + v*s(t) + z(:,t); end

%% decode the signal by linear regression
minsampleindx = round(0.1*tmax); % let the response statistics stabilize
samplesize = 1e5; % number of observations used for regression
for k=0:kmax
    trainsamples = randsample(minsampleindx:tmax,samplesize);
    s_train = s(trainsamples-k);
    X_train = X(:,trainsamples);
    wts = regress(s_train,X_train'); % obtain readout weights by regressing
    
    wts_decode(k+1,:) = wts;
    var_decode(k+1) = var(wts'*X_train);
    mse_decode(k+1) = mean((s_train' - wts'*X_train).^2);
    corr_decode(k+1) = corr(s_train , (wts'*X_train)');
    
    var_optimal(k+1) = var(wts_optimal(k+1,:)*X_train);
    mse_optimal(k+1) = mean((s_train' - wts_optimal(k+1,:)*X_train).^2);
    corr_optimal(k+1) = corr(s_train , (wts_optimal(k+1,:)*X_train)');
end

%% decode the signal by linear regression (cross-validated)
% minsampleindx = round(0.1*tmax); % let the response statistics stabilize
% samplesize = 1e5; % number of observations used for regression
% for k=0:kmax
%     trainsamples = randsample(minsampleindx:tmax,samplesize);
%     s_train = s(trainsamples-k);
%     X_train = X(:,trainsamples);
%     wts = regress(s_train,X_train'); % obtain readout weights by regressing
%     testsamples = randsample(setdiff(minsampleindx:tmax,trainsamples),samplesize);
%     s_test = s(testsamples-k);
%     X_test = X(:,testsamples);
%     
%     wts_decode(k+1,:) = wts;
%     var_decode(k+1) = var(wts'*X_test);
%     mse_decode(k+1) = mean((s_test' - wts'*X_test).^2);
% end