function networkparams = BuildNetwork(globalparams,architecture)

%% network parameters
networkparams.name = architecture;
networkparams.N = globalparams.N; 
networkparams.alpha = globalparams.alpha; 
networkparams.kmax = globalparams.kmax;

N = networkparams.N; alpha = networkparams.alpha;

%% construct weight matrix
switch lower(networkparams.name)
    case 'delay ring'
        W = sqrt(alpha)*circshift(eye(N),1);
    case '2d lattice'
        m = round(sqrt(N));
        if abs(sqrt(N)-m)<1e-10
            [XX,YY] = meshgrid(1:m,1:m);
            x_neighbour = ismember(abs(XX(:) - XX(:)') , [1 m-1]);
            x_identical = XX(:) - XX(:)' == 0;
            y_neighbour = ismember(abs(YY(:) - YY(:)') , [1 m-1]);
            y_identical = YY(:) - YY(:)' == 0;
            W = sqrt(alpha)*((x_neighbour & y_identical) | (y_neighbour & x_identical));
            W = W/mean(sum(W>0)); % normalise by mean #connections for stability
        else, error(['N must be a perfect square for building ' networkparams.name]);
        end
    case 'random symmetric'
        W = normrnd(0,sqrt(alpha/(4*N)),[N,N]);
        W = triu(W,1) + triu(W,1)' + diag(diag(W)); % symmetrise
    case 'delay line'
        W = sqrt(alpha)*circshift(eye(N),1);
        W(end,:) = 0;
    case 'random'
        W = normrnd(0,sqrt(alpha/(1.1*N)),[N,N]);
end
networkparams.W = W;

%% eigen decomposition
[V,D] = eig(W);
networkparams.eigen_vecs = V;
networkparams.eigen_vals = diag(D);