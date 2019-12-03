function PlotRecurrentMemory(networkparams,inputs,outputs)

%% load parameters
networkname = networkparams.name;
N = networkparams.N;
W = networkparams.W;
kmax = networkparams.kmax;
eigen_vals = networkparams.eigen_vals;

%% load inputs
v = inputs.v;
tau_signal = inputs.tau_signal;
tau_noise = inputs.tau_noise;

%% load ouputs
Cnoise = outputs.Cnoise;
FMM_theory = outputs.FMM_theory;
FMC_theory = outputs.FMC_theory;
FMC_decode = outputs.FMC_decode;
% J_theory2 = outputs.J_theory2;

%% plot
figure; set(gcf,'Position',[100 100 1400 700]);
hold on; colormap hot;
sgtitle(networkname,'FontWeight','Bold');

subplot(2,4,1); hold on;
imagesc(W); axis tight;
xlabel('Unit'); ylabel('Unit'); title('Weight matrix');

subplot(2,4,2); hold on;
if ~isreal(eigen_vals)
    plot(real(eigen_vals),imag(eigen_vals),'ok');
    axis([-2 2 -2 2]);
    vline(0,'--r'); hline(0,'--r');
    xlabel('Re (\lambda)'); ylabel('Im (\lambda)'); title('Eigenspectrum');
else
    plot(eigen_vals,'ok');
    ylim([-2 2]);
    xlabel('Eigenmode'); ylabel('\lambda'); title('Eigenspectrum');
end

subplot(2,4,3); hold on;
imagesc(FMM_theory); axis tight;
xlabel('Unit'); ylabel('Unit'); title('Fisher Memory Matrix (FMM)');

subplot(2,4,4); hold on;
plot(0:kmax,FMC_decode(tau_signal==0,:)); plot(0:kmax,FMC_theory,'--k');
ylim([0 max(FMC_theory)]);
legend('Simulated','Theoretical','Fontsize',14);
xlabel('Time lag'); ylabel('Fisher Information'); title('Fisher Memory Curve (FMC)');

subplot(2,4,5); hold on;
plot(inputs.stats(1).lags,cell2mat({inputs.stats.temporalcorr_signal}));
xlabel('Time lag, k - l'); ylabel('Correlation, s_{kl}'); title('Autocorrelogram of input signal');

subplot(2,4,6); hold on;
plot(0:kmax,FMC_decode); set(gca,'YScale','log');
xlabel('Time lag'); ylabel('Fisher Information'); title('Simulated FMC');

subplot(2,4,7); hold on;
plot(tau_signal, sum(FMC_decode,2),'-sk','MarkerFaceColor','r');
xlabel('Timescale of input corr, \tau'); ylabel('Total Fisher Information');

% subplot(2,4,8); hold on;
% plot(sum(J_theory2,2), sum(J_decode,2),'-sk','MarkerFaceColor','r');
% xlabel('Timescale of input corr, \tau'); ylabel('Total Fisher Information');
