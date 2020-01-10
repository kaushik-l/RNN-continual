function [inputs,outputs,errflag] = ComputeFisherMemory(globalparams,networkparams)

%% network parameters
W = networkparams.W; kmax = networkparams.kmax; N = networkparams.N;
infinity = 1e3; % define large time for estimating response covariance

%% input properties
inputs.v = globalparams.v;
inputs.tau_signal = globalparams.tau_signal;
inputs.var_signal = globalparams.var_signal;
inputs.tau_noise = globalparams.tau_noise;
inputs.var_noise = globalparams.var_noise;

v = inputs.v; % direction of input
tau_signal = inputs.tau_signal; % timescale of temporal correlation in input
var_signal = inputs.var_signal; % input variance
tau_noise = inputs.tau_noise; % timescale of temporal correlation in noise
var_noise = inputs.var_noise; % noise variance

%% estimate noise covariance
Cnoise = arrayfun(@(k) (W^k)*(W^k)', 0:infinity, 'UniformOutput', false);
Cnoise = sum(cat(3,Cnoise{:}),3);
Cnoiseinv = inv(Cnoise);
outputs.Cnoise = Cnoise;
outputs.Cnoiseinv = Cnoiseinv;

%% estimate signal covariance
% build filter
sz = infinity + 1; %filter size
trange = linspace(0, infinity, sz);
h = @(tau) exp(-trange.^2/(2*tau^2));
% estimate signal covariance
for i=1:numel(tau_signal)
    if tau_signal(i)>0, g = h(tau_signal(i)); else, g = zeros(1,sz); g(1) = 1; end
    S = toeplitz(g);
    P = cell2mat(arrayfun(@(k) (W^k)*v,0:infinity,'UniformOutput',false));
    Csig = P*S*P';
    Csiginv = inv(Csig);
    Cinv = inv(Cnoise + Csig);
    
    outputs.S{i} = S;
    outputs.Csig{i} = Csig;
    outputs.Csiginv{i} = Csiginv;
    outputs.C{i} = Cnoise + Csig;
    outputs.Cinv{i} = Cinv;
end

%% calculate theoretical fisher memory matrix FMM (this does not depend on Csig)
for k=0:kmax, for l=0:kmax, J(k+1,l+1) = v'*(W^k)'*Cnoiseinv*(W^l)*v; end; end
outputs.J = J;

%% calculate Cramer-Rao Bound for correlated signal 
for i=1:numel(tau_signal)
    S = outputs.S{i}; S = S(1:kmax+1, 1:kmax+1);
    J_tilde = real(S^.5)*J*real(S^.5);
    B = real(S^.5)*(J_tilde*((eye(kmax+1) + J_tilde))^-1)*real(S^-.5) - eye(kmax+1);
    M = real(S^.5)*(J_tilde*((eye(kmax+1) + J_tilde))^-1)*real(S^.5);
    outputs.CRB{i} = real(S^.5)*(J_tilde*((eye(kmax+1) + J_tilde))^-2)*real(S^.5);
    outputs.MSE{i} = outputs.CRB{i} + B'*B;
    outputs.J_tilde{i} = J_tilde;
    outputs.B{i} = B;
    outputs.M{i} = M;
end

%% calculate optimal readout weights
for i=1:numel(tau_signal)
    P = cell2mat(arrayfun(@(k) (W^k)*v,0:infinity,'UniformOutput',false));
    S = outputs.S{i};
    U = S*P'*outputs.Cinv{i};
    outputs.Uopt{i} = U(1:kmax+1,:);
end

%% simulate fisher info (by decoding)
FMC_decode = nan(numel(tau_signal),kmax+1);
for i=1:numel(tau_signal)
    fprintf(['Simulating tau_signal = ' num2str(tau_signal(i)) '\n']);
    [inputs.stats(i),wts_decode{i},var_decode{i},mse_decode{i},corr_decode{i},...
        var_optimal{i},mse_optimal{i},corr_optimal{i}] = ...
        SimulateRecurrent(v,W,tau_signal(i),var_signal,tau_noise,var_noise,kmax,outputs.Uopt{i});
end
outputs.wts_decode = wts_decode;
outputs.var_decode = var_decode;
outputs.mse_decode = mse_decode;
outputs.corr_decode = corr_decode;
outputs.var_optimal = var_optimal;
outputs.mse_optimal = mse_optimal;
outputs.corr_optimal = corr_optimal;

errflag = any(isnan(outputs.var_decode{1})); % happens if random matrix is not well designed

%% function
function sigcov = SignalCov(timelag,g,t,kmax)
if abs(timelag) > kmax, sigcov = 0;
else, sigcov = g(abs(t - timelag) < 1e-10); end