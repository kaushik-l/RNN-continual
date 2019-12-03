%% define network parameters
globalparams.N = 36;
globalparams.alpha = 0.9;
globalparams.kmax = 25;

%% define input parameters
v = zeros(globalparams.N,1); v(1) = 1; v = v/(v'*v); 
globalparams.v = v; % input direction
globalparams.tau_signal = 0:4; % timescale of temporal correlation in input
globalparams.var_signal = 1; % input variance
globalparams.tau_noise = 0; % timescale of temporal correlation in noise
globalparams.var_noise = 1; % noise variance

%% analyse various architectures
architectures = {'Delay Ring','2D Lattice','Random Symmetric','Delay Line','Random'};
N_archs = numel(architectures);
for k=1:N_archs
    fprintf(['**********Analysing ' architectures{k} '**********\n']);
    networkparams(k) = BuildNetwork(globalparams,architectures{k});
    [inputs(k),outputs(k)] = ComputeFisherMemory(globalparams,networkparams(k));
    PlotRecurrentMemory(networkparams(k),inputs(k),outputs(k));
end

%% compare memory capacity of various architectures
figure; set(gcf,'Position',[100 100 1400 700]);
hold on;
for k=1:N_archs
    subplot(1,2,1); hold on;
    plot(globalparams.tau_signal,sum(outputs(k).FMC_decode,2),'-s','MarkerFaceColor','k');
    legend(architectures,'FontSize',14);
    xlabel('Timescale of input corr, \tau'); ylabel('Total Fisher Information');
    
    subplot(1,2,2); hold on;
    plot(globalparams.tau_signal,sum(outputs(k).FMC_decode,2)/max(sum(outputs(k).FMC_decode,2)),...
        '-s','MarkerFaceColor','k');
    legend(architectures,'FontSize',14);
    xlabel('Timescale of input corr, \tau'); ylabel('Normalized total Fisher Information');
end