function [inputstats,J_decode] = SimulateFMC(v,W,tau_signal,var_signal,tau_noise,var_noise,kmax)

N = numel(v); % number of units

%% build filter
sz = 4*max(tau_signal,tau_noise); %filter size
t2 = linspace(-sz/2, sz/2, sz);
h = @(tau) exp(-t2.^2/(2*tau^2))/sum(exp(-t2.^2/(2*tau^2)));

%% simulate the dynamics
tmax = 3e5;
s = normrnd(0,sqrt(var_signal),tmax,1);
if tau_signal > 0, s = conv(s,h(tau_signal),'same'); end
z = normrnd(0,sqrt(var_noise),N,tmax);
if tau_noise > 0, for i=1:N, z(i,:) = conv(z(i,:),h(tau_noise),'same'); end; end
[inputstats.temporalcorr_signal, inputstats.lags] = xcorr(s(:),kmax,'coeff');
inputstats.temporalcorr_noise = xcorr(z(:),kmax,'coeff');
X = zeros(N,tmax); 
for t = 2:tmax, X(:,t) = W*X(:,t-1) + v*s(t) + z(:,t); end

%% decode the signal by linear regression
minsampleindx = round(0.1*tmax); % let the response statistics stabilize
samplesize = 1e5; % number of observations used for regression
J_decode = nan(1,1+kmax);
for k=0:kmax
    trainsamples = randsample(minsampleindx:tmax,samplesize);
    s_train = s(trainsamples-k);
    X_train = X(:,trainsamples);
    w_opt = regress(s_train,X_train'); % obtain readout weights by regressing
    w_opt = w_opt/(w_opt'*((W^k)*v)); % normalise for unbiasedness
    testsamples = randsample(setdiff(minsampleindx:tmax,trainsamples),samplesize);
    s_test = s(testsamples-k);
    X_test = X(:,testsamples);
    
    J_decode(k+1) = 1/var(s_test' - w_opt'*X_test);
end