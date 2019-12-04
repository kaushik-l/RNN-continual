%% define network parameters
globalparams.N = 36; % network graph visualization is readily interpretable for N=4
globalparams.alpha = 0.9; % connection strength

architecture = '2D Lattice'; % 'Delay Ring','2D Lattice','Random Symmetric','Delay Line','Random'

%% build network
networkparams = BuildNetwork(globalparams,architecture);

%% plot network graph, connectivity matrix and eigenspectrum
PlotNetworkParams(networkparams);

%% functions
function networkparams = BuildNetwork(globalparams,architecture)

%% network parameters
networkparams.name = architecture;
networkparams.N = globalparams.N;
networkparams.alpha = globalparams.alpha;

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

end

function PlotNetworkParams(networkparams)

%% load parameters
networkname = networkparams.name;
N = networkparams.N;
W = networkparams.W;
eigen_vals = networkparams.eigen_vals;

%% plot
figure; set(gcf,'Position',[50 50 1400 400]);
hold on; colormap hot;
sgtitle(networkname,'FontWeight','Bold');

% plot graph
subplot(1,3,1); hold on;
plot(digraph(W),'Layout','subspace3');
axis off; title('Network Graph');

% plot connectivity matrix
subplot(1,3,2); hold on;
imagesc(W); axis tight;
xlabel('Unit'); ylabel('Unit'); title('Weight matrix');

% plot eigenspectrum
subplot(1,3,3); hold on;
if ~isreal(eigen_vals)      % for networks with complex eigenvalues
    plot(real(eigen_vals),imag(eigen_vals),'ok');
    axis([-2 2 -2 2]);
    line([0 0],[-2 2],'Color','r'); line([-2 2],[0 0],'Color','r');
    xlabel('Re (\lambda)'); ylabel('Im (\lambda)'); title('Eigenspectrum');
else                        % for networks with real eigenvalues
    plot(eigen_vals,'ok');
    ylim([-2 2]);
    xlabel('Eigenmode'); ylabel('Eigenvalue, \lambda'); title('Eigenspectrum');
end

end