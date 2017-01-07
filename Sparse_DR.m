%% Hyperparameter
N = 400; m = 100;

opts.mu = 1;
opts.gamma = 1;
opts.nb_it_max = 10000;
opts.verbose = true;

% Set random generators for reproducibility
seed = 192017 ;
rand('state',seed) ;
randn('state',seed) ;

%% Douglas-Rachford sparse recovery
% Operators
L = randn(m,N) / sqrt(m); 
op.pL = L'*pinv( L * L');
op.proxG = @(x,gamma)max(0,1-gamma./max(1e-15,abs(x))).*x; 
op.rproxG = @(x,tau)2*op.proxG(x,tau)-x;
op.proxF = @(x,y)x + op.pL*(y-L*x);
op.rproxF = @(x,y)2*op.proxF(x,y)-x;

% Sparse signals
q = 1000; slist = 1:1:10; % number of signals, sparsity list
Slist = slist(mod(0:q-1,length(slist))+1);
U = rand(N,q); v = sort(U); v = v( (0:q-1)*N + Slist );
a_init = U <= repmat( v, [N 1] ); % |a_init(:,j)| has sparsity |Slist(j)|.

% Acquisition
tmp = L*a_init;
b = abs(L*a_init);
phase_init = tmp./b;

% Douglas-Rachford algorithm
[ a, phi, err, norm_a ] = DR_descent( L, b, op, opts );

% Recovering probabilty
proba = zeros(1,length(slist));
E = mean(abs(a-a_init))<.05;
for j=1:length(slist)
    s = slist(j);
    proba(j) = mean(E(Slist==s));
end

% Display
figure(1); clf;
h = plot(slist, proba, 'k');
axis([min(slist) max(slist) -.05 1.05]);
title('Probability to recover a signal given it sparsity')
xlabel('Sparsity')
ylabel('Proba of recovery')

%% Cross-validation 
% norms_a = zeros(opts.nb_it,9); errors = zeros(opts.nb_it,9);
% legends = {}; k=0;
% for gamma = [.01 1 10]
% for mu = [.001 .01 1]
%     k=k+1; opts.gamma = gamma; opts.mu = mu;
%     [ a, phi, err, norm_a ] = DR_descent( L, b, op, opts );
%     norms_a(:,k) = norm_a';
%     legends{end+1} = ['\gamma=' num2str(gamma) ', \mu=' num2str(mu)];
% end
% end
% figure(2); clf;
% h = plot(log10(norms_a(1:end/2,:)-min(norms_a(:))));
% set(h, 'LineWidth', 2); legend(h, legends); 