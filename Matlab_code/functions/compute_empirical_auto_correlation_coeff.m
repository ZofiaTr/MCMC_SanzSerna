function rho_nu = compute_empirical_auto_correlation_coeff(X, lag)


    rho_nu = zeros(1, length(lag));
    
    count=0;
    
    for nu = lag
        
        count=count+1;
        
        rho_nu(count) = autocovariance(X, nu);
    end
    
    rho_nu = rho_nu ./ rho_nu(1);


end