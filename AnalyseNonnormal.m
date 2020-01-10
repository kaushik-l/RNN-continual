kmax = 50;
signalparams.density = (1:N)/N; % signal density
signalparams.corr = 0:0.2:1;
signalparams.location = 1:N;
noiseparams.var = 1; % noise variance

networkparams.N = 2^3;
networkparams.chainlength = 2.^(0:log2(N));
networkparams.gamma = 0.9;
networkparams.beta = [0 1];
networkparams.alpha = 1:10;

for beta = networkparams.beta
   if beta == 0, alpha = 1;
       W = DesignNonNormal(gamma,beta,alpha);
       [inputstats,wts_decode,var_decode,mse_decode,corr_decode] = ...
           DecodeSignal(v,W,var_signal,var_noise,kmax,[]);
   else
       for alpha = networkparams.alpha
           W = DesignNonNormal(gamma,beta,alpha);
           [inputstats,wts_decode,var_decode,mse_decode,corr_decode] = ...
               DecodeSignal(v,W,var_signal,var_noise,kmax,[]);
       end
   end
end