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

J = outputs.J;
Jdiag = diag(J);
CRB = outputs.CRB;
MSE = outputs.MSE;

var_decode = outputs.var_decode;
mse_decode = outputs.mse_decode;

var_optimal = outputs.var_optimal;
mse_optimal = outputs.mse_optimal;

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
imagesc(J); axis tight;
xlabel('Unit'); ylabel('Unit'); title('Fisher Memory Matrix (FMM)');

subplot(2,4,4); hold on;
plot(0:kmax,var_decode{tau_signal==0}); plot(0:kmax,diag(CRB{tau_signal==0}),'--k');
ylim([0 max(var_decode{tau_signal==0})]);
legend('Decoded','Cramer-Rao Lower Bound (CRB)','Fontsize',14);
xlabel('Time lag'); ylabel('Variance'); title('Estimator variance');

% subplot(2,4,5); hold on;
% plot(0:kmax,cell2mat(arrayfun(@(k) diag(CRB{k})',1:5,'UniformOutput',false)'),'--');
% xlabel('Time lag'); ylabel('Variance'); t = title({'Cramer-Rao'; 'Lower Bound (CRB)'});
% set(t, 'horizontalAlignment', 'right');

% subplot(2,4,5); hold on;
% plot(0:kmax,cell2mat(var_decode'));
% xlabel('Time lag'); ylabel('Variance'); title('Decoded');

subplot(2,4,5); hold on;
plot(0:kmax,1-cell2mat(mse_decode'));
xlabel('Time lag'); ylabel('R^2'); t = title({'Decoded by'; 'regression'});
set(t, 'horizontalAlignment', 'right');

subplot(2,4,6); hold on;
plot(tau_signal, mean(1-cell2mat(mse_decode'),2),'-sk','MarkerFaceColor','b');
ylabel('Mean R^2');
yyaxis right;
[~,bestlag] = max(1-cell2mat(mse_decode'),[],2); bestlag = bestlag-1;
plot(tau_signal, bestlag,'-sk','MarkerFaceColor','r');
legend({'Mean R^2','Best lag'},'Fontsize',14)
xlabel('Timescale of input corr, \tau'); ylabel('Best lag');
title({'Decoded by'; 'regression'});

subplot(2,4,7); hold on;
for k=1:numel(tau_signal), scatter(1-mse_optimal{k},1-mse_decode{k},'s','Filled'); end
plot(0:.1:.5,0:.1:.5,'--k'); axis([0 .5 0 .5]);
title('R^2'); ylabel({'Decoded by'; 'regression'}); xlabel({'Optimal decoder'});

subplot(2,4,8); hold on;
[~,bestlag_opt] = max(1-cell2mat(mse_optimal'),[],2); bestlag_opt = bestlag_opt-1;
for k=1:numel(tau_signal), scatter(bestlag_opt(k),bestlag(k),100,'s','Filled'); end
plot(0:max(tau_signal),0:max(tau_signal),'--k'); axis([0 max(tau_signal) 0 max(tau_signal)]);
title('Best lag'); ylabel({'Decoded by'; 'regression'}); xlabel({'Optimal decoder'});

axes('Position',[.26 .38 .05 .1]); box on; hold on;
plot(inputs.stats(1).lags,cell2mat({inputs.stats.temporalcorr_signal}));
xlabel('Time lag, k - l'); ylabel('Correlation, s_{kl}');
