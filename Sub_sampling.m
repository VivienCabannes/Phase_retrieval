%% Preprocessing
% Options
name = 'images/Ek_Balam.jpg';
% load('opts.mat')
opts.nb_it_max = 10000;
opts.eps = 10^(-7);
opts.calcul_obj = true;
% opts.init_guess = phi;

% Set random generators for reproducibility
seed = 192017 ;
rand('state',seed) ;
randn('state',seed) ;

% Read the two images
X_init = imread(name);
N_init = size(X_init,1);

%% Fourier transform
% % Simulate measurements
% tmp = fft2(X_init);
% b = abs(tmp);
% phase_init = tmp./b;
% 
% % Recover a relative phase up
% [ f, phi, obj ] = fourier_AM( b, opts );
% 
% figure(1);clf;
% imagesc(reshape(real(f), size_init));
% 
% figure(2);clf;
% imagesc(reshape(X, size_init));
% 
%% Random filter
N = 16; % Size of testing sample
nb_ex = 100;
success = zeros(nb_ex+1,11);

for k = 1:11
m = 100*(k-1); display(m);
success(end,k) = m/(N*N); % Compression factor

for j = 1:nb_ex
% Take a sample of size N*N
ind1 = max(min(ceil(N_init*rand), N_init-N),1);
ind2 = max(min(ceil(N_init*rand), N_init-N),1);
X = X_init(ind1:ind1+N,ind2:ind2+N);

% Vectorialization
size_init = size(X);
X = double(X(:));

% Acquisition
L = randn(m,size(X,1)) + 1i * randn(m,size(X,1)); opts.L = L;
opts.A = pinv(L);
tmp = L*X;
b = abs(tmp);
phase_init = tmp./b;

% Retrieval
[ ~, phi, ~ ] = general_AM( L, b, opts );
success(j,k) = mean(abs(phi-phase_init) < .05);
end
end

%% Display results
suc=mean(success(1:end-1,:)>.9,1);
figure(3); clf;
plot(success(end,:),suc);
xlabel('Sub-sampling factor')
ylabel('Average recovery')
title('Probability of recovery regarding compression factor')
saveas(3,'proba_recovery','epsc')

% save('opts.mat', 'opts') % useful if A was slow to compute