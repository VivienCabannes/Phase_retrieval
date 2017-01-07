function [ f, phi, obj ] = fourier_AM( b, opts )
% Reconstruct a signal from measurement b = |Lf|
    try
        nb_it = opts.nb_it_max;
        eps = opts.eps;
    catch
        nb_it = 5000;
        eps = 10^(-7);
    end
    
    obj = zeros(1,nb_it);
    
    % Initialization
    try
        phi = opts.init_guess;
    catch
        phi = exp(2*pi*(1i)*rand(size(b)));
    end
    
    % While a convergence criterion is not reached
    k = 0; current_obj = Inf; fprintf('Iterations: ')
    while current_obj > eps && k < nb_it 
        k=k+1;
        if mod(k,50)==0
            fprintf('%d, ', k);
        end
        % Update signal
        f = ifft2(b.*phi);
        
        % Cast the function toward real
        tmp = phase(f);
        phase_mean = phase(mean(tmp(:)));
        if ~isnan(phase_mean)
            f = f./phase_mean;
        end
            
        % Cast the function toward real intensity
        ind_neg = f<0;
        ind_pos = f>255;
        f(ind_neg) = .5*f(ind_neg);
        f(ind_pos) = 255 + .5*(f(ind_pos)-255);

        % Update phase
        tmp = fft2(f);
        phi = phase(tmp);
        
        % Calcul objectif if ask
        if opts.calcul_obj
            tmp = abs(tmp - b.*phi);
            current_obj = sum(tmp(:).^2);
            obj(k) = current_obj;
        end
    end
    fprintf(' done.\n')
end