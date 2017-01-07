function [ a, phi, energy, norm_a ] = FW_descent( L, b, op, opts )
% This is a first version, 
% It was abandon since sparsity doesn't really improved the practical
% problems we were looking at
    
    % Collecting arguments
    PsiS = op.PsiS; U = op.U; Psi = op.Psi; Thresh = op.Thresh;
    niter = opts.nb_it_max; gamma = opts.gamma; lambda = opts.lambda;
    im_size = opts.size_init; A = opts.A;
    
    % Initialize the phase and a first guess
    try
        phi = opts.init_guess;
    catch
        if isreal(L)
            phi = sign(randn(size(b)));
        else
            phi = exp(2*pi*(1i)*rand(size(b)));
        end
    end
    y = b.*phi; f = A * y; a = U.*PsiS(vec2im(f, opts.size_init)); tmp = L*f - y;

    energy = zeros(size(b,2), niter); norm_a = zeros(size(b,2), niter); 
    if opts.verbose
        fprintf('Iteration: 0, ')
    end
    lambdas = linspace(lambda, 0, niter);
    for i=1:niter
        if opts.verbose && ~mod(i,500)
            fprintf('%d, ', i)
        end
        
        a = Thresh( a - gamma*PsiS(vec2im((L')*tmp,opts.size_init)), lambdas(i)*gamma );
        f = Psi(a); 
   
        % Cast the function toward real phase
        tmp = phase(f);
        phase_mean = phase(mean(tmp(:)));
        if ~isnan(phase_mean)
            f = f./phase_mean;
        end
          
        % Cast the function toward real intensity
%         ind_neg = f<0;
%         ind_pos = f>255;
%         f(ind_neg) = .5*f(ind_neg);
%         f(ind_pos) = 255 + .5*(f(ind_pos)-255);
        
        % Update
        tmp = L*(f(:)); phi = phase(tmp);  y = b.*phi;
        
        % Register energy
        norm_a(i) = sum( abs(a(:)) );
        energy(i) = 1/2*sqrt(sum((tmp-y).^2)) + lambdas(i) * norm_a(i);   
    end
    if opts.verbose
        fprintf('done.\n')
    end
end

