%% define network parameters
globalparams.N = 81;
globalparams.alpha = 0.9;
globalparams.kmax = 50;

%% define input parameters
v = zeros(globalparams.N,1); v(1) = 1; v = v/(v'*v); 
globalparams.v = v; % input direction
globalparams.tau_signal = 0:8:32; % timescale of temporal correlation in input
globalparams.var_signal = 1; % input variance
globalparams.tau_noise = 0; % timescale of temporal correlation in noise
globalparams.var_noise = 1; % noise variance

%% analyse various architectures
architectures = {'Delay Ring','2D Lattice','Random Symmetric','Delay Line','Random'};
N_archs = numel(architectures);
for k=1:N_archs
    fprintf(['**********Analysing ' architectures{k} '**********\n']);
    repeat = 1;
    while repeat
        networkparams(k) = BuildNetwork(globalparams,architectures{k});
        [inputs(k),outputs(k),repeat] = ComputeFisherMemory(globalparams,networkparams(k));
    end
    PlotRecurrentMemory(networkparams(k),inputs(k),outputs(k));
end

%% compare memory capacity of various architectures
figure; set(gcf,'Position',[100 100 1400 250]);
hold on;
for k=1:N_archs
    subplot(1,N_archs,k); hold on;
    plot(0:globalparams.kmax,1-cell2mat(outputs(k).mse_decode'),'linewidth',2);
    xlabel('Time lag'); ylabel('Mean R^2');
    title(architectures{k});
    axis([0 globalparams.kmax 0 1]);
end
leg = legend(arrayfun(@(k) num2str(k), globalparams.tau_signal, 'UniformOutput', false));
title(leg,'corr timescale');