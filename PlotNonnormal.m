function PlotNonnormal(params,signalparams,noiseparams,decoder,pcdecoder)


%% load network params
chainlen = params.chainlen; chainlens = unique(chainlen); numchainlens = numel(chainlens); 
minchainlen = min(chainlen); maxchainlen = max(chainlen);
beta = params.beta; numbetas = numel(unique(beta));  minbeta = min(beta); maxbeta = max(beta);
alpha = params.alpha; alphas = unique(alpha); numalphas = numel(alphas); 
minalpha = min(alpha); maxalpha = max(alpha);
gamma = unique(params.gamma);
N = unique(params.N);

%% load signal parameters
density = params.density; densities = unique(density); numdensities = numel(densities);
mindensity = min(density); maxdensity = max(density);
rho = params.rho; rhos = unique(rho); numrhos = numel(rhos);
minrho = min(rho); maxrho = max(rho);
loc = params.loc; minloc = min(loc);

%% identify cases correspond to minimum values of different parameters
minchainnets = (chainlen==minchainlen); maxchainnets = (chainlen==maxchainlen);
minbetanets = (beta==minbeta); maxbetanets = (beta==maxbeta);
minalphanets = (alpha == minalpha); maxalphanets = (alpha == maxalpha);
mindensitynets = (density == mindensity); maxdensitynets = (density == maxdensity);
minrhonets = (rho == minrho);
minlocnets = (loc == minloc);

%% baseline case
baselinenet = find(minchainnets & minbetanets & minalphanets & mindensitynets & minrhonets & minlocnets);

%% effect of network parameters
figure; hold on;
cmap = [[1 0 0]; [1 0.5 0.5]; winter(numalphas-1)];
for i=1:numchainlens
    % select networks of different chain lengths
    nets = find(chainlen==chainlens(i) & mindensitynets & minrhonets & minlocnets);  
    % plot each chain length
    subplot(2,ceil(numchainlens/2),i); hold on;
    for j=1:numel(nets), plot(0:signalparams.maxlag(:),decoder.corr{nets(j)},'linewidth',2); 
    lgnd{j} = ['\it \gamma=' num2str(gamma)...
        '\it \beta=' num2str(beta(nets(j))) ', \it \alpha=' num2str(alpha(nets(j)))]; 
    end
    set(gca,'Colororder',cmap); axis([0 signalparams.maxlag 0 0.5]);
    % add title
    if i==1, title(['chainlength, \it l = ' num2str(chainlens(i))],'Interpreter','Latex','Fontsize',14);
    else, title(['\it l = ' num2str(chainlens(i))],'Interpreter','Latex','Fontsize',14); end
    % add legend
    if i==numchainlens, legend(lgnd,'Fontsize',10); end
    % add axis labels
    xlabel('Time lag, \it $$\Delta$$', 'Interpreter', 'Latex','Fontsize',14); 
    ylabel('Corr(\it s,$$\hat{s}$$)', 'Interpreter', 'Latex','Fontsize',14);
end
clear lgnd;

%% effect of network parameters
figure; hold on;
cmap = summer(numchainlens);
for i=1:numalphas
    % select networks of different alphas
    nets = find(beta>0 & alpha==alphas(i) & mindensitynets & minrhonets & minlocnets);  
    % plot each alpha
    subplot(2,ceil(numalphas/2),i); hold on;
    for j=1:numel(nets), plot(0:signalparams.maxlag(:),decoder.corr{nets(j)},'linewidth',2); 
    lgnd{j} = ['\it l=' num2str(chainlen(nets(j)))]; 
    end
    set(gca,'Colororder',cmap); axis([0 signalparams.maxlag 0 0.5]);
    % add title
    title(['\gamma = ' num2str(gamma) ', \beta = ' num2str(beta(nets(j))) ...
        ', \alpha = ' num2str(alphas(i))],'Fontsize',14);
    % add legend
    if i==numalphas, legend(lgnd,'Fontsize',10,'Interpreter','Latex'); end
    % add axis labels
    xlabel('Time lag, \it $$\Delta$$', 'Interpreter', 'Latex','Fontsize',14); 
    ylabel('Corr(\it s,$$\hat{s}$$)', 'Interpreter', 'Latex','Fontsize',14);
end
clear lgnd;

%% effect of network parameters
figure; hold on;
% select networks of different params
nets = find(mindensitynets & minrhonets & minlocnets);
r2 = mean(log10([decoder.r2{nets}])); r2 = r2(:);
r2 = [r2(1)*ones(numalphas,1) ; r2];
r2Mat = reshape(r2,[numalphas+numbetas-1 numchainlens]);
imagesc(chainlens,[0 ; alphas],r2Mat); colormap(parula); axis tight;
% add title
title('Decoding precision, $$log(R^2) = g(k,\alpha)$$','Interpreter','Latex','Fontsize',18);
xlabel('chainlength, \it k','Interpreter','Latex','Fontsize',18); 
ylabel('degree of non-normality, $$\alpha$$','Interpreter','Latex','Fontsize',18);
h = colorbar; h.Label.String = 'log\it(R^2)'; h.FontSize = 14; hlim = h.Limits;
clear lgnd;

%% effect of signal properties
figure; hold on;
cmap = [brewermap(numrhos,'Greys') ; brewermap(numrhos,'Reds') ; ...
    brewermap(numrhos,'Greens') ; brewermap(numrhos,'Blues')];
% select signals of different params
for i=1:numalphas
    for j=1:numchainlens
        if chainlens(j)==1, nets = find(beta==0 & alpha==alphas(i) & chainlen==chainlens(j) & minlocnets);
        else, nets = find(beta>0 & alpha==alphas(i) & chainlen==chainlens(j) & minlocnets); end
        if ~isempty(nets)
            [~,indx] = sort(density(nets)); nets = nets(indx); % sort by density
            % plot each d,rho
            subplot(numalphas,numchainlens,numchainlens*(i-1)+j); hold on;
%             figure; hold on;
            for k=1:numel(nets), plot(0:signalparams.maxlag(:),mean(decoder.corr{nets(k)},2),'linewidth',0.25);
                lgnd{k} = ['\it d=' num2str(density(nets(k))) '$$\rho$$=' num2str(rho(nets(k)))];
            end
            set(gca,'Colororder',cmap); axis([0 signalparams.maxlag 0 0.5]);
            set(gca,'XScale','log'); set(gca,'YScale','log');
        end
    end
    % add legend
    if i==numalphas && j==numchainlens, legend(lgnd,'Fontsize',10,'Interpreter','Latex'); end
end
clear lgnd;

%% effect of signal properties
% cmap = [[1 0 0]; [1 0.5 0.5]; winter(numalphas-1)];
% % select signals of different params
% for i=1:numdensities
%     for j=1:numrhos
%         figure; hold on;
%         for k=1:numchainlens
%             % select networks of different chain lengths
%             nets = find(chainlen==chainlens(k) & density==densities(i) & rho==rhos(j) & minlocnets);
%             % plot each chain length
%             subplot(2,ceil(numchainlens/2),k); hold on;
%             for m=1:numel(nets), plot(0:signalparams.maxlag(:),mean(decoder.corr{nets(m)},2),'linewidth',2);
%                 lgnd{m} = ['\it \gamma=' num2str(gamma)...
%                     '\it \beta=' num2str(beta(nets(m))) ', \it \alpha=' num2str(alpha(nets(m)))];
%             end
%             set(gca,'Colororder',cmap); axis([0 signalparams.maxlag 0 0.5]);
%             % add title
%             if k==1, title(['chainlength, \it l = ' num2str(chainlens(i))],'Interpreter','Latex','Fontsize',14);
%             else, title(['\it l = ' num2str(chainlens(k))],'Interpreter','Latex','Fontsize',14); end
%             % add legend
%             if k==numchainlens, legend(lgnd,'Fontsize',10); end
%             % add axis labels
%             xlabel('Time lag, \it $$\Delta$$', 'Interpreter', 'Latex','Fontsize',14);
%             ylabel('Corr(\it s,$$\hat{s}$$)', 'Interpreter', 'Latex','Fontsize',14);
%         end
%         sgtitle(['\it d = ' num2str(densities(i)) ', $$\rho$$ = ' num2str(rhos(j))]);
%     end
% end
% clear lgnd;

%% effect of signal properties
figure; hold on;
% select signals of different params
for i=1:numdensities
    numinputs = densities(i)*N;
    for j=1:numrhos
        nets = find(density==densities(i) & rho==rhos(j) & minlocnets);
        r2 = mean(real(log10([decoder.r2{nets}]))); % average across timelags
        r2 = reshape(r2,[numinputs numel(r2)/numinputs]); r2 = mean(r2,1); % average across signal components
        r2 = [r2(1)*ones(numalphas,1) ; r2(:)];
        r2Mat = reshape(r2,[numalphas+numbetas-1 numchainlens]);
        subplot(numdensities,numrhos,numrhos*(i-1)+j);
        imagesc(chainlens,[0 ; alphas],r2Mat,hlim); colormap(parula); axis tight;
        set(gca,'YDir','normal')
        % add title
        title(['\it d = ' num2str(densities(i)) ', $$\rho$$ = ' num2str(rhos(j))],...
            'Interpreter','Latex','Fontsize',12);
        xlabel('\it k','Interpreter','Latex','Fontsize',14);
        ylabel('$$\alpha$$','Interpreter','Latex','Fontsize',14);
    end
end
% add global title
sgtitle('Decoding precision, $$log(R^2) = g(k,\alpha,\it d, \rho)$$',...
    'Interpreter','Latex','Fontsize',18);

%% effect of signal properties
figure; hold on;
% select signals of different params
for i=1:numdensities
    numinputs = densities(i)*N;
    for j=1:numrhos
        nets = find(density==densities(i) & rho==rhos(j) & minlocnets);
        corr_decode = mean([decoder.corr{nets}]); % average across timelags
        corr_decode = reshape(corr_decode,[numinputs numel(corr_decode)/numinputs]); 
        corr_decode = mean(corr_decode,1); % average across signal components
        corr_decode = [corr_decode(1)*ones(numalphas,1) ; corr_decode(:)];
        corr_decodeMat = reshape(corr_decode,[numalphas+numbetas-1 numchainlens]);
        subplot(numdensities,numrhos,numrhos*(i-1)+j);
        imagesc(chainlens,[0 ; alphas],corr_decodeMat,[0 0.2]); colormap(parula); axis tight;
        set(gca,'YDir','normal')
        % add title
        title(['\it d = ' num2str(densities(i)) ', $$\rho$$ = ' num2str(rhos(j))],...
            'Interpreter','Latex','Fontsize',12);
        xlabel('\it k','Interpreter','Latex','Fontsize',14);
        ylabel('$$\alpha$$','Interpreter','Latex','Fontsize',14);
    end
end
% add global title
sgtitle('Decoding precision, $$Corr(\it s, \hat{s}) = g(k,\alpha,\it d, \rho)$$',...
    'Interpreter','Latex','Fontsize',18);

%% effect of signal properties
figure; hold on;
% select signals of different params
for i=1:numdensities
    for j=1:numrhos
        nets = find(density==densities(i) & rho==rhos(j) & minlocnets);
        [~,bestnet(i,j,:)] = ...
            max(cell2mat(cellfun(@(x) mean(x,2), decoder.corr(nets), 'UniformOutput', false)'),[],2);
        subplot(numdensities,numrhos,numrhos*(i-1)+j); hold on;
        plot(chainlen(squeeze(bestnet(i,j,:))),'-s'); ylim([0 8]);
        yyaxis right; plot(alpha(squeeze(bestnet(i,j,:))),'-s'); ylim([0 10]);     
    end
end