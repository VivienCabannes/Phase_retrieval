function [ a, phi, err, norm_a ] = DR_descent( L, b, op, opts )

    % Recover usual object
    niter = opts.nb_it_max; mu = opts.mu; gamma = opts.gamma ;
    rproxG = op.rproxG ; rproxF = op.rproxF ; proxF = op.proxF;
    N = size(L, 2);
    
    % Initialize the phase
    try
        phi = opts.init_guess;
    catch
        if isreal(L)
            phi = sign(randn(size(b)));
        else
            phi = exp(2*pi*(1i)*rand(size(b)));
        end
    end
    % Main loop
    err = zeros(size(b,2), niter); norm_a = zeros(size(b,2), niter); 
    tmp_a = zeros(N,1);
    if opts.verbose
        fprintf('Iteration: 0, ')
    end
    for i=1:niter
        if opts.verbose && ~mod(i,50)
            fprintf('%d, ', i)
        end
        y = b.*phi;
        tmp_a = (1-mu/2)*tmp_a + mu/2*rproxG(rproxF(tmp_a, y), gamma);
        a = proxF(tmp_a, y);
        tmp = L*a;
        err(:,i) = sum(abs(b-abs(tmp)).^2,1); norm_a(:,i) = sum(abs(a),1);
        phi = phase(tmp);
    end
    if opts.verbose
        fprintf('done.\n')
    end
end

