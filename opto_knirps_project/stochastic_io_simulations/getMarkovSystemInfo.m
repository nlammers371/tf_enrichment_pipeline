function simInfo = getMarkovSystemInfo(simInfo)

    simType = simInfo.simType;

    simInfo.memory = 7;
    simInfo.deltaT = 20;
    simInfo.t_MS2 = 1.4;
    
    % specify 2 state network architecture (eventually this will be drawn from
    % actual fits)
    simInfo.R2 = [-.92  1/1.07; 
                   .92 -1/1.07]/60;    

    r2 = 4.6784e4;    
    simInfo.r2 = [0 r2];
    simInfo.pi0 = [0.5 0.5];
    simInfo.noise = 5e3; % NL: this isn't really doing anything atm