function param_val_prop = make_mh_proposal(param_val_curr,param_bounds)    
    prop_sigma = 0.25*param_val_curr;
%     prop_sigma = 0.1*diff(param_bounds);
    % 'Z' from the non-standard Gaussian N(m,s^2)
    % X=trandn((l-m)/s,(u-m)/s) and set Z=m+s*X;
    prop_temp = trandn((param_bounds(1)-param_val_curr)/prop_sigma,...
                (param_bounds(2)-param_val_curr)/prop_sigma);
    param_val_prop = param_val_curr + prop_temp*prop_sigma;