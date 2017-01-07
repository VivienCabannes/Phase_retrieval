function [ f, phi, energy, norm_f ] = FW_descent_fourier_ortho( b, op, opts )
    
    % Collecting arguments
    proxG = op.proxG; gradientF = op.gradientF; PsiS = op.PsiS;
    niter = opts.nb_it_max; gamma = opts.gamma; lambda = opts.lambda;
    
    % Initialize the phase and a first guess
    try
        phi = opts.init_guess;
    catch
        phi = sign(randn(size(b)));
    end
    y = b.*phi; f = fft2(y); 
    energy = zeros(1, niter); norm_f = zeros(1, niter); 
    if opts.verbose
        fprintf('Iteration: 0, ')
    end
    lambdas = linspace(lambda, 0, niter);
    for i=1:niter
        if opts.verbose && ~mod(i,500)
            fprintf('%d, ', i)
        end
        %f = proxG( f - gamma*gradientF(f,y), gamma*lambdas(i) ); 
        f = proxG(fft2(y), gamma*lambdas(i));
        % Cast the function toward real phase
        tmp = phase(f);
        phase_mean = phase(mean(tmp(:)));
        if ~isnan(phase_mean)
            f = f./phase_mean;
        end
          
        % Cast the function toward real intensity
        f = ceil(255*(f-min(abs(f(:))))/max(abs(f(:))));
%        ind_neg = f<0;
%        ind_pos = f>255;
%        f(ind_neg) = .5*f(ind_neg);
%        f(ind_pos) = 255 + .5*(f(ind_pos)-255);
        
        % Update
        tmp = fft2(f); phi = phase(tmp); y = b.*phi;
        
        % Register energy
        norm_f(i) = norm(PsiS(f), 1);
        energy(i) = 1/2*sum(abs(tmp(:) - y(:)).^2) + lambdas(i)*norm_f(i);     
    end
    if opts.verbose
        fprintf('done.\n')
    end
end

