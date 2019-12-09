function [fluo_mean_vec, fluo_se_vec, offset_mean_vec, offset_se_vec, nc14_flag, quantiles, tlims] = ...
                mcp_saturation_analysis(Prefix,varargin)

% specify defaults
quantiles = [1 25 50 75 98]/100; % quantiles to calculate
nBoots = 100;
nc14_flag = false;
tlims = [0 Inf];
% check for optional arguments
for i = 1:numel(varargin)
    if ischar(varargin{i}) && i < numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

% get path info, spot channel, and movie database
[rawDataPath,ProcPath,DropboxFolder,MS2CodePath, PreProcPath,...
    rawDataFolder, Prefix, ExperimentType,Channel1,Channel2,OutputFolder,...
    Channel3, spotChannels, MovieDataBaseFolder, movieDatabase]= readMovieDatabase(Prefix);

% load Particle info
cp_path = [DropboxFolder, filesep, Prefix, filesep, 'CompiledParticles.mat'];
if exist(cp_path)
    load(cp_path);
    if iscell(CompiledParticles)
        CompiledParticles = CompiledParticles{1};
    end

    % extract fluorescence and offset values
    fluo_vec = [CompiledParticles.Fluo];
    offset_vec = [CompiledParticles.Off];
    time_vec = ElapsedTime([CompiledParticles.Frame]);
    time_ft = true(size(time_vec));
    if nc14_flag
        if isempty(nc14)
            warning('no nuclear cycle info found')
            nc14 = 1;
        end
        nc14_time = ElapsedTime(nc14);
        time_vec_adjusted = time_vec - nc14_time;
        time_ft = time_vec_adjusted >= tlims(1) & time_vec_adjusted < tlims(2);    
    end

    if sum(time_ft) > 100
        fluo_vec = fluo_vec(time_ft);
        offset_vec = offset_vec(time_ft);

        % initialize arrays
        fluo_boot_array = NaN(nBoots,numel(quantiles));
        offset_boot_array = NaN(nBoots,numel(quantiles));
        index_vec = 1:numel(fluo_vec);
        N = numel(index_vec);
        % perform bootstrapping
        for n = 1:nBoots
            % take bootstrap sample
            boot_ids = randsample(index_vec,N,true);
            fluo_vec_boot = fluo_vec(boot_ids);
            offset_vec_boot = offset_vec(boot_ids);
            % calculate quantiles
            for q = 1:numel(quantiles)
                fq1 = quantile(fluo_vec_boot,quantiles(q));
                fq2 = quantile(fluo_vec_boot,quantiles(q)+.01);
                % find nearest match
                q_ind = fluo_vec_boot>=fq1&fluo_vec_boot<fq2;
                fluo_boot_array(n,q) = nanmean(fluo_vec_boot(q_ind));
                offset_boot_array(n,q) = nanmean(offset_vec_boot(q_ind));
            end
        end

        % calculate bootstrap means and standard errors
        fluo_mean_vec = nanmean(fluo_boot_array);
        fluo_se_vec = nanstd(fluo_boot_array);

        offset_mean_vec = nanmean(offset_boot_array);
        offset_se_vec = nanstd(offset_boot_array);
    else
        warning('Too few data points met criteria. Skipping...')
        fluo_mean_vec = NaN(1,numel(quantiles));
        fluo_se_vec = NaN(1,numel(quantiles));

        offset_mean_vec = NaN(1,numel(quantiles));
        offset_se_vec = NaN(1,numel(quantiles));
    end
else
    warning('No CompiledParticles structure found...')
    fluo_mean_vec = NaN(1,numel(quantiles));
    fluo_se_vec = NaN(1,numel(quantiles));

    offset_mean_vec = NaN(1,numel(quantiles));
    offset_se_vec = NaN(1,numel(quantiles));
end