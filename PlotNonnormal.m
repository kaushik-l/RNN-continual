function PlotNonnormal(params,signalparams,noiseparams,decoder,pcdecoder)


%% load network params
chainlen = params.chainlen; numchainlens = numel(unique(chainlen)); minchainlen = min(chainlen); maxchainlen = max(chainlen);
beta = params.beta; numbetas = numel(unique(beta));  minbeta = min(beta); maxbeta = max(beta);
alpha = params.alpha; numalphas = numel(unique(alpha)); minalpha = min(alpha); maxalpha = max(alpha);

%% load signal parameters
density = params.density; mindensity = min(density); maxdensity = max(density);
rho = params.rho; minrho = min(rho);
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

%% plot effect of chain length
nets = find(maxbetanets & maxalphanets & mindensitynets & minrhonets & minlocnets);
figure; hold on;
plot(0:signalparams.maxlag(:),decoder.r2{baselinenet});
for i=1:numel(nets), plot(0:signalparams.maxlag(:),decoder.corr{nets(i)}); end

%% plot effect of nonnormality (beta/alpha)
nets = find(maxchainnets & mindensitynets & minrhonets & minlocnets);
figure; hold on;
for i=1:numel(nets), plot(0:signalparams.maxlag(:),decoder.corr{nets(i)}); end

%% plot effect of network parameters
nets = find(mindensitynets & minrhonets & minlocnets);
nrows = numchainlens; ncols = numalphas + numbetas - 1;

figure(1); hold on;
for i=1:numel(nets)
    subplot(nrows,ncols,(chainlen(nets(i))>1)*(ncols-1) + i);
    plot(0:signalparams.maxlag(:),decoder.corr{nets(i)},'linewidth',2); 
    axis([0 signalparams.maxlag 0 0.5]);        
end

figure(2); hold on;
for i=1:numel(nets)
    subplot(nrows,ncols,(chainlen(nets(i))>1)*(ncols-1) + i);
    plot(0:signalparams.maxlag(:),log10(decoder.r2{nets(i)}),'linewidth',2);
%     axis([0 signalparams.maxlag 1 1.5]);
end