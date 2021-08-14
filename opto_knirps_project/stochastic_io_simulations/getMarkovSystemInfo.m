function simInfo = getMarkovSystemInfo(simInfo)

    simType = simInfo.simType;
    
%     simInfo = struct;
    
    simInfo.memory = 7;
    simInfo.deltaT = 20;
    simInfo.t_MS2 = 1.4;
    
    % specify 2 state network architecture (eventually this will be drawn from
    % actual fits)
    simInfo.R2 = [-.92  1/1.07; 
                        .92 -1/1.07]/60;    

%     % estimate r for now
%     pon = systemInfo.R2(2,1) / (systemInfo.R2(2,1) + systemInfo.R2(1,2));
%     f_mean = nanmean(nanmean(io_ref_struct.fluo_raw_array(1:30,:),2),1);
%     fluo_kernel = ms2_loading_coeff (systemInfo.t_MS2, systemInfo.memory);
%     f_mean / sum(fluo_kernel) / pon;
    r2 = 4.6784e4;%
    % systemParams.r2 = [0 4]*1e4; % loading rate for each state (in units of au per time step)
    simInfo.r2 = [0 r2];
    simInfo.pi0 = [0.5 0.5];
    if contains(simType,'_on') && ~contains(simType,'2')
        simInfo.pi0 = [0 0];
    end
    simInfo.noise = 5e3; % NL: this isn't really doing anything atm