function [inputs,outputs] = ComputeFisherMemory(globalparams,networkparams)

%% network parameters
W = networkparams.W; kmax = networkparams.kmax; N = networkparams.N;
infinity = 1e2; % define large time for estimating response covariance

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
% % build filter
% sz = 2*kmax; %filter size
% trange = linspace(-sz/2, sz/2, sz+1);
% h = @(tau) exp(-trange.^2/(2*tau^2));
% % estimate signal covariance
% for i=1:numel(tau_signal)
%     Csig = zeros(N);
%     if tau_signal(i)>0, g = h(tau_signal(i)); else, g = zeros(1,sz+1); g(sz/2 + 1) = 1; end
%     for k=0:infinity
%         for l=0:infinity
%             Csig = Csig + (W^k)*v*SignalCov(k-l,g,trange,kmax)*v'*(W^l)';
%         end
%     end
%     Csiginv = inv(Csig);
%     Ctotinv = inv(Cnoise + Csig);
%     
%     outputs.Csig{i} = Csig;
%     outputs.Csiginv{i} = Csiginv;
%     outputs.Ctot{i} = Cnoise + Csig;
%     outputs.Ctotinv{i} = Ctotinv;
% end

%% calculate theoretical fisher info using noise cov only (Ganguli 2008)
FMC_theory = arrayfun(@(k) v'*(W^k)'*Cnoiseinv*(W^k)*v, 0:kmax);
outputs.FMC_theory = FMC_theory;

for k=0:kmax, for l=0:kmax, FMM_theory(k+1,l+1) = v'*(W^k)'*Cnoiseinv*(W^l)*v; end; end
outputs.FMM_theory = FMM_theory;

%% calculate theoretical fisher info using total cov
% for i=1:numel(tau_signal)
%     J_theory2 = arrayfun(@(k) v'*(W^k)'*outputs.Ctotinv{i}*(W^k)*v, 0:kmax);
%     outputs.J_theory2(i,:) = J_theory2;
% end

%% simulate fisher info (by decoding)
FMC_decode = nan(numel(tau_signal),kmax+1);
for i=1:numel(tau_signal)
    fprintf(['Simulating tau_signal = ' num2str(tau_signal(i)) '\n']);
    [inputs.stats(i),FMC_decode(i,:)] = SimulateFMC(v,W,tau_signal(i),var_signal,tau_noise,var_noise,kmax);
end
outputs.FMC_decode = FMC_decode;

%% function
function sigcov = SignalCov(timelag,g,t,kmax)
if abs(timelag) > kmax, sigcov = 0;
else, sigcov = g(abs(t - timelag) < 1e-10); end