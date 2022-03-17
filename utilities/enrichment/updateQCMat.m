function updateQCMat
%         qc_mat(numel(nc_x_vec)-j_pass+1).setID = currentSetID; 
%         qc_mat(numel(nc_x_vec)-j_pass+1).frame = currentFrame;
%         qc_mat(numel(nc_x_vec)-j_pass+1).nc_index = nc_lin_index_vec(j);
%         qc_mat(numel(nc_x_vec)-j_pass+1).nc_sub_index = nc_sub_index_vec(j);
%         qc_mat(numel(nc_x_vec)-j_pass+1).qcFlag = edge_qc_flag_vec(j);            
%         qc_mat(numel(nc_x_vec)-j_pass+1).xp = x_spot;
%         qc_mat(numel(nc_x_vec)-j_pass+1).yp = y_spot;  
%         qc_mat(numel(nc_x_vec)-j_pass+1).xp_sister = x_spot_sister;
%         qc_mat(numel(nc_x_vec)-j_pass+1).yp_sister = y_spot_sister;  
%         qc_mat(numel(nc_x_vec)-j_pass+1).xc_edge = edge_null_x_vec(j);
%         qc_mat(numel(nc_x_vec)-j_pass+1).yc_edge = edge_null_y_vec(j);
%         qc_mat(numel(nc_x_vec)-j_pass+1).xc_serial = serial_null_x_vec(j);
%         qc_mat(numel(nc_x_vec)-j_pass+1).yc_serial = serial_null_y_vec(j);       
%         qc_mat(numel(nc_x_vec)-j_pass+1).particleID = particle_id_vec(j);
%         qc_mat(numel(nc_x_vec)-j_pass+1).serial_qc_flag = serial_qc_flag_vec(j);
%         qc_mat(numel(nc_x_vec)-j_pass+1).edge_qc_flag = edge_qc_flag_vec(j);        
%         sz = nb_size;
%         edge_dist_mat = nc_dist_frame;
%         edge_dist_mat(~spot_nc_mask&~null_mask) = 0;
%         if edge_qc_flag_vec(j) == 2         
%             sz = max([nb_size,abs(x_nucleus - edge_null_x_vec(j)),abs(y_nucleus - edge_null_y_vec(j))...
%                 abs(x_nucleus - serial_null_x_vec(j)),abs(y_nucleus - serial_null_y_vec(j))]);            
%         end
%         y_range = max(1,y_nucleus-sz):min(yDim,y_nucleus+sz);
%         x_range = max(1,x_nucleus-sz):min(xDim,x_nucleus+sz);
%         qc_mat(numel(nc_x_vec)-j_pass+1).x_origin = x_range(1);
%         qc_mat(numel(nc_x_vec)-j_pass+1).y_origin = y_range(1);
%         qc_mat(numel(nc_x_vec)-j_pass+1).mcp_snip = mcp_slice(y_range,x_range);
%         qc_mat(numel(nc_x_vec)-j_pass+1).protein_snip = protein_slice(y_range,x_range);
%         qc_mat(numel(nc_x_vec)-j_pass+1).edge_dist_snip = edge_dist_mat(y_range,x_range);  


% disp('saving qc frames...')
% % save qc data
% tic
% particle_index = unique([spot_struct_protein.particleID]);
% particle_index = particle_index(~isnan(particle_index));
% qc_particles = randsample(particle_index,min([100,numel(particle_index)]),false);
% particle_index_full = [];
% particle_frames_full = [];
% for i = 1:numel(qc_structure)
%     qc_mat = qc_structure(i).qc_mat;
%     for  j = 1:numel(qc_mat)
%         qc_spot = qc_mat(j);
%         if ~isfield(qc_spot,'ParticleID')
%             continue
%         end
%         ParticleID = qc_spot.particleID;
%         if isempty(ParticleID) || ~ismember(ParticleID,qc_particles)
%             continue
%         end        
%         samplingInfo.Frame = qc_spot.frame;      
%         particle_index_full = [particle_index_full ParticleID];
%         particle_frames_full = [particle_frames_full samplingInfo.Frame];        
%         save_name = [snipPath 'pt' num2str(1e4*ParticleID) '_frame' sprintf('%03d',samplingInfo.Frame) '.mat'];
%         save(save_name,'qc_spot');
%     end
% end
% [particle_index_full, si] = sort(particle_index_full);
% particle_frames_full = particle_frames_full(si);