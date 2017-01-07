function [ f, phi, energy, norm_f ] = FW_descent_ortho( L, b, op, opts )
    
    % Collecting arguments
    proxG = op.proxG; gradientF = op.gradientF; PsiS = op.PsiS;
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
    y = b.*phi; f = A * y; 
    energy = zeros(size(b,2), niter); norm_f = zeros(size(b,2), niter); 
    if opts.verbose
        fprintf('Iteration: 0, ')
    end
    lambdas = linspace(lambda, 0, niter);
    for i=1:niter
        if opts.verbose && ~mod(i,500)
            fprintf('%d, ', i)
        end
        f = proxG( vec2im(f - gamma*gradientF(f,y), im_size), gamma*lambdas(i) ); f = f(:);
        %f = proxG(vec2im(A*y,im_size), gamma*lambda); f = f(:);
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
        tmp = L*f; phi = phase(tmp); y = b.*phi;
        
        % Register energy
        norm_f(i) = norm(PsiS(vec2im(f, im_size)), 1);
        energy(i) = 1/2*sum(abs(tmp - y).^2,1) + lambda*norm_f(i);     
    end
    if opts.verbose
        fprintf('done.\n')
    end
end

