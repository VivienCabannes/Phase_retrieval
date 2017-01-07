% Set random generators for reproducibility
seed = 192017 ;
rand('state',seed) ;
randn('state',seed) ;

% Enable toolbox
getd = @(p)path(p,path); 
getd('toolbox_signal/');
getd('toolbox_general/');

%% Hyperparameter
% Options
opts.gamma = 1;
opts.lambda = 1;
opts.nb_it_max = 50;
opts.verbose = true;

%% Read an image
name = 'images/Valladolid.jpg';
X = imread(name);
N = size(X,1);

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

%% Acquisition
tmp = fft2(X);
b = abs(tmp);
phase_init = tmp./b;

%% Phase retrieval
% First on orthogonal wavelet
op.proxG = SoftThreshPsi_ortho;
op.gradientF = @(f,y) ifft2( fft2(f) - y );
op.PsiS = PsiS_ortho;
%opts.nb_it_max = floor(opts.nb_it_max/1000)+1;
[ f, phi, E, ~ ] = FW_descent_fourier_ortho( b, op, opts );

imwrite(uint8(X), 'FW_fourier_init.jpg')
imwrite(uint8(real(f)), 'FW_fourier_retrieve.jpg')

% Refined with translation invariant wavelet
% To be done