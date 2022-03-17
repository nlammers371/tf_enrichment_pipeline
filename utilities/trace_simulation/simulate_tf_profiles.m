function tf_profile_array = simulate_tf_profiles(frac_init,frac_final,simInfo)

    KD = simInfo.KD;
    HC = 5;%simInfo.HC;
    deltaT = simInfo.deltaT;
    seq_length = simInfo.seq_length;
    n_traces = simInfo.n_traces;
    
    % generate hypothetical temporal TF input profile (assume linear)        
    r_vec = [KD*(frac_init/(1-frac_init))^(1/HC) KD*(frac_final/(1-frac_final))^(1/HC)];
    if frac_final > frac_init     
        r0 = max(r_vec);
        r1 = min(r_vec);        
    elseif frac_final < frac_init
        r0 = min(r_vec);
        r1 = max(r_vec);
    end    
    
    t_shift = 10*60/deltaT;
    center_time = seq_length/2;
    tf_profile_array = NaN(seq_length+1,n_traces);
    tf_profile_array(1:center_time-t_shift/2,:) = r0;
    tf_profile_array(center_time+2+t_shift/2:end,:) = r1;
    tf_profile_array(center_time+1-t_shift/2:center_time+1+t_shift/2,:) = repmat(linspace(r0,r1,t_shift+1)',1,n_traces);
    tf_profile_array(tf_profile_array<0) = 1e-6;