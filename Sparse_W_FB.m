% Set random generators for reproducibility
seed = 192017 ;
rand('state',seed) ;
randn('state',seed) ;

% Enable toolbox
getd = @(p)path(p,path); 
getd('toolbox_signal/');
getd('toolbox_general/');

%% Hyperparameter
N = 128; % size of images: N*N
m = N*N;

% Options
opts.gamma = .01;
opts.lambda = 10;
opts.nb_it_max = 10;
opts.verbose = true;

%% Read an image
name = 'images/Valladolid.jpg';
X_init = imread(name);
N_init = size(X_init,1);

%% Operators
SoftThresh = @(x,gamma)x.*max( 0, 1-gamma./max(abs(x),1e-10) );
HardThresh = @(x,gamma)x.*(abs(x)>gamma);

% Wavelet dictionaries and operators
% 1. Orthogonal wavelet
Jmax = log2(N)-1; Jmin = Jmax-3; options.ti = 0;
PsiS_ortho = @(f)perform_wavelet_transf(f, Jmin, +1,options);
Psi_ortho = @(a)perform_wavelet_transf(a, Jmin, -1,options);
%   - Thresholding in this dictonary
SoftThreshPsi_ortho = @(f,T)Psi_ortho(SoftThresh(PsiS_ortho(f),T));
HardThreshPsi_ortho = @(f,T)Psi_ortho(HardThresh(PsiS_ortho(f),T));

% 2. Translation invariant wavelet dictionary
option.ti = 1;
Xi = @(a)perform_wavelet_transf(a, Jmin, -1,options);
PsiS = @(f)perform_wavelet_transf(f, Jmin, +1,options);
J = Jmax-Jmin+1;  u = [4^(-J) 4.^(-floor(J+2/3:-1/3:1)) ];
U = repmat( reshape(u,[1 1 length(u)]), [N N 1] );
Psi = @(a)Xi(a./U);

%% Random filter
try 
    load('L.mat');
    load('L_inv.mat');
catch
    L = randn(m,N*N) / sqrt(m);
    fprintf('Inverting acquisition matrix: '); save('L.mat', 'L');
    L_inv = pinv(L); save('L_inv.mat', 'L_inv'); 
    fprintf('done.\n');
end
opts.A = L_inv;

%% Acquisition
% Choose a patch on the image as signal
ind1 = min(ceil(N_init*rand), N_init-N);
ind2 = min(ceil(N_init*rand), N_init-N);
X = X_init(ind1+1:ind1+N,ind2+1:ind2+N);

% Vectorialization
opts.size_init = size(X);
X = double(X(:));

% Acquisition
tmp = L*X;
b = abs(tmp);
phase_init = tmp./b;

%% Phase retrieval
% First on orthogonal wavelet
op.proxG = HardThreshPsi_ortho;
op.gradientF = @(f,y) conj(L') * ( L * f - y );
op.PsiS = PsiS_ortho;
opts.nb_it_max = floor(opts.nb_it_max/10)+1;
[ f, phi, ~, ~ ] = FW_descent_ortho( L, b, op, opts );
opts.phase_init = phi;
imwrite(uint8(vec2im(X, opts.size_init)), 'FW_init.jpg')
imwrite(uint8(vec2im(f, opts.size_init)), 'FW_retrieve_ortho.jpg')

% Refined with translation invariant wavelet
op.PsiS = PsiS; op.U = U; op.Psi = Psi; 
op.Thresh = SoftThresh;
opts.nb_it_max = opts.nb_it_max*10;
[ f, phi, energy, norm_a ] = FW_descent( L, b, op, opts );
imwrite(uint8(vec2im(f, opts.size_init)), 'FW_retrieve.jpg')

