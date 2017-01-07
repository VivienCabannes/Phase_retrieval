function [ phi ] = phase( X )
    phi = X./abs(X);
    phi(isnan(phi)) = exp(2*pi*(1i)*rand);
end

