% script to extract stripe positions over time and build lookup tables
clear
close all

% list of projects to process
projectList = {'eveWT','EveGtSL','EveS1Null','EveGtSL-S1Null'};

% pointers to subfolders with desired stripe info
customFolder = 'DataFeb2021_AP';
customSubFolderList = {'eveS1wt_eveS2wt_AP','eveS1wt_eveS2Gt_AP','eveS1Null_eveS2wt_AP','eveS1Null_eveS2Gt_AP'};

% cell arrays indicating reference pairings for each stripe
ref_struct = struct;
% WT
ref_struct(1).from_stripes = 0;
ref_struct(1).to_stripe = 1;
% Delta Gt
ref_struct(2).from_stripes = [0 1.5];
ref_struct(2).to_stripe = 1;
% S1 Null
ref_struct(3).from_stripes = [0 1 -1];
ref_struct(3).to_stripe = 2;
% both
ref_struct(4).from_stripes = [0 -1 1 1.5];
ref_struct(4).to_stripe = 2;

for p = 1:length(projectList)
    
    % load the compiled spots structure
    liveProject = LiveEnrichmentProject(projectList{p});
    load([liveProject.dataPath filesep 'spot_struct.mat']);
    
    % extract time grid for interpolation
    interp_time_vec = [spot_struct.timeInterp];
    Tres = nanmedian(diff(interp_time_vec));
    interp_grid = 0:Tres:60*60;%unique(interp_time_vec(~isnan(interp_time_vec)));
    
    % intialize lookup table
    set_index = unique([spot_struct.setID]);
    stripe_array_master = NaN(length(set_index)*length(interp_grid),9); % one column for each stripe. Col 1 := time, col 2:= setID    
    stripe_array_ref = NaN(length(set_index)*length(interp_grid),9); % this tracks position of reference stripes over time
            
    xy_to_ap = [];
    
    if liveProject.hasTracesCompiled
      % pull normalized AP positions from custom-segmented data sets
      numExperiments = length(liveProject.includedExperiments);
      resultsDir = liveProject.includedExperiments{1}.userResultsFolder;
      
      customDataPath = [resultsDir filesep customFolder filesep customSubFolderList{p} filesep];

      % iterate through each individual experiment
      for e = 1:numExperiments
        currExperiment = liveProject.includedExperiments{e};

        % search for corresponding custom data file
        customFile = dir([customDataPath filesep '*' currExperiment.Prefix '*']);
        if length(customFile) == 1
          
          % load custom file
          load([customFile.folder filesep customFile.name]);
          if iscell(CompiledParticles)
            CompiledParticles = CompiledParticles{1};
          end          
                    
          % extract stripe positions
          temp_set_array = [];
          xy_ap_array = [];
          for c = 1:length(CompiledParticles)
              ap_vec = CompiledParticles(c).centerPseudoAP;
              ap_vec_particle = CompiledParticles(c).APPos;
              ap_vec_ref = CompiledParticles(c).centerPseudoAP_ALL;
              if ~isempty(ap_vec)
                frame_vec = CompiledParticles(c).Frame;
                time_vec = ElapsedTime(frame_vec) - ElapsedTime(nc14);
                stripe_vec = CompiledParticles(c).Stripe;
                
                temp_set_array = [temp_set_array ; [60*time_vec' stripe_vec' ap_vec' ap_vec_ref']];  
                xy_ap_array = [xy_ap_array ; [CompiledParticles(c).xPos' CompiledParticles(c).yPos' ap_vec_particle']];  
              end
          end
          
          % build linear model for AP inference
          mdl = fitlm(xy_ap_array(:,1:2),xy_ap_array(:,3));
          if mdl.Rsquared.Ordinary < 0.99
              error('problem with pixel-to-ap conversion')
          end
          xy_to_ap(e,:) = mdl.Coefficients.Estimate;
%           if e == 5
%             error('asfa')
%           end
          
          % take unique set
          temp_array_u = unique(temp_set_array,'rows');
          stripe_index = 1:7;%unique(temp_set_array(:,2));
          
          % interpolate
          set_pos_array = NaN(length(interp_grid),9);
          set_ref_array = NaN(length(interp_grid),9);
          for s = 1:length(stripe_index)
              s_ind = stripe_index(s) + 2;
              stripe_filter = temp_array_u(:,2)==stripe_index(s);
              if any(stripe_filter)
                  ap_vec_stripe = temp_array_u(stripe_filter,3);

                  set_pos_array(:,s_ind) = interp1(temp_array_u(stripe_filter,1),temp_array_u(stripe_filter,3),interp_grid);              

                  set_ref_array(:,s_ind) = interp1(temp_array_u(stripe_filter,1),temp_array_u(stripe_filter,4),interp_grid); 
              end
          end
          
          % stripe-specific
          set_pos_array(:,1) = interp_grid;
          set_pos_array(:,2) = e;
          
          % "average" stripe
          set_ref_array(:,1) = interp_grid;
          set_ref_array(:,2) = 0;
          
          % add to master set
          start_i = 1 + (e-1)*length(interp_grid);
          stop_i = e*length(interp_grid);
          stripe_array_master(start_i:stop_i,:) = set_pos_array;
          
          % add to ref structure
          stripe_array_ref(start_i:stop_i,:) = set_ref_array;
          
        else
          error('Duplicate CompiledParticles files found')
        end
      end
      
      % convert set-specific array to table
      stripe_position_table = array2table(stripe_array_master,'VariableNames',...
                        {'time','setID','stripe_1_ap','stripe_2_ap','stripe_3_ap',...
                         'stripe_4_ap','stripe_5_ap','stripe_6_ap','stripe_7_ap'});
      % save
      writetable(stripe_position_table,[liveProject.dataPath filesep 'stripe_position_key.csv']);
      
      % generate condensed output reference table
      stripe_ref_array_cd = NaN(length(interp_grid),8);
      for  t = 1:length(interp_grid)
%           t_table = stripe_array_ref(stripe_array_ref(:,1)==interp_grid(t),3:end);
          t_table = stripe_array_master(stripe_array_master(:,1)==interp_grid(t),3:end);
          stripe_ref_array_cd(t,2:end) = nanmean(t_table);
      end
      stripe_ref_array_cd(:,1) = interp_grid;
      stripe_ref_table = array2table(stripe_ref_array_cd,'VariableNames',...
                        {'time','stripe_1_ap','stripe_2_ap','stripe_3_ap',...
                         'stripe_4_ap','stripe_5_ap','stripe_6_ap','stripe_7_ap'});
      % save
      writetable(stripe_ref_table,[liveProject.dataPath filesep 'stripe_ref_key.csv']);
      
      % update normalized AP positions
      for s = 1:length(spot_struct)
         % create old version of field if appropriate
         if ~isfield(spot_struct,'APPosParticleNormOrig') && isfield(spot_struct,'APPosParticleNorm')
            spot_struct(s).APPosParticleNormOrig = spot_struct(s).APPosParticleNorm;
         elseif isfield(spot_struct,'APPosParticleNormOrig') && isfield(spot_struct,'APPosParticleNorm')
            if isempty(spot_struct(s).APPosParticleNormOrig)
              spot_struct(s).APPosParticleNormOrig = spot_struct(s).APPosParticleNorm;
            end
         end
         
         % get original AP coordinates         
         time_vec = spot_struct(s).time;
         setID = spot_struct(s).setID;
         if false
%          if false~isnan(spot_struct(s).particleID) && spot_struct(s).qcFlag
%             time_vec_interp = spot_struct(s).timeInterp;
%             ap_vec_interp = spot_struct(s).APPosParticleInterp;        
         
         else
            start_i = find(interp_grid>=time_vec(1),1);
            stop_i = find(interp_grid<=time_vec(end),1,'last');
            time_vec_interp = interp_grid(start_i:stop_i);
            xRaw = spot_struct(s).xPosNucleus;
            yRaw = spot_struct(s).yPosNucleus;
            ap_raw = xy_to_ap(setID,1) + xy_to_ap(setID,2)*double(xRaw) + xy_to_ap(setID,3)*double(yRaw);
            
            if length(ap_raw) > 1
                ap_vec_interp = interp1(time_vec,ap_raw,time_vec_interp);
            else                
                [~, mi] = min(abs(time_vec-time_vec_interp));
                ap_vec_interp = repelem(ap_raw,length(time_vec_interp));                
            end
            
         end
         if any(isnan(ap_vec_interp))
           error('wtf')
         end                  
         
         % extract relevant chunk of array
         setID = spot_struct(s).setID;
         
         % set-specific
         stripe_set_array = stripe_array_master(stripe_array_master(:,2)==setID,:);
         stripe_set_array = stripe_set_array(ismember(interp_grid,time_vec_interp),3:end);
         
         % average         
         stripe_ref_array = stripe_ref_array_cd(ismember(interp_grid,time_vec_interp),2:end);
         
         % remove designated stripes for each genotype to prevent us from
         % using them as references
         stripe_set_array(:,ismember(1:7,ref_struct(p).from_stripes)) = NaN;
         
         % find closest stripe for each time point
         [~,nearest_stripe_vec] = min(abs(ap_vec_interp' - stripe_set_array),[],2);
         lin_index_vec = sub2ind(size(stripe_ref_array),(1:size(stripe_ref_array,1))',nearest_stripe_vec);
         
         % calculate ap adjustment
         ap_vec_new = ap_vec_interp + stripe_ref_array(lin_index_vec)' - stripe_set_array(lin_index_vec)';
         
         % update
         spot_struct(s).APPosParticleNorm = 100*ap_vec_new;         
                  
      end
      
      % save spot struct
      save([liveProject.dataPath filesep 'spot_struct.mat'],'spot_struct')
    end
end    