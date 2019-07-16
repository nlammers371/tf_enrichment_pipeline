% main03_check_control_selection(project)
%
% DESCRIPTION
% Generates figures overlaying segmentation results onto histone channel
% and interface in which user can accept or reject spot assignments
%
% ARGUMENTS
% project: master ID variable
%
% OPTIONS
% dropboxFolder: Path to data folder where you wish to save
%                pipeline-generated data sets and figures. If this
%                var is not specified, output will be saved one level
%                above git repo in the folder structure
%
% INTERFACE
%
% Movement:
%
% n: move back one sample
% m: move forward one sample
% j: jump to specified sample (entered in command line)
%
% x: Exit
%
% General instructions: For now this script is primarily intended as a
% means to spot-check the segementation and sample selection for different
% data sets and time points. If systematic issues are uncovered, we will
% need to adapt a different component of the pipleine either to resolve
% issue or to remove problematic observations

function main03_check_control_selection(project,varargin)
close all
% specify paths
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];
figPath = [dropboxFolder 'LocalEnrichmentFigures\' project '/'];
for i = 1:numel(varargin)    
    if strcmpi(varargin{i}, 'dropboxFolder')        
        dataPath = [varargin{i+1} '/ProcessedEnrichmentData/' project '/'];
        figPath = [varargin{i+1} '/LocalEnrichmentFigures/' project '/control_selection/'];
    end
end

snipPath = [dataPath 'qc_images/'];
mkdir(figPath);

% load data
load([dataPath '/qc_ref_struct.mat']);
particle_index_full = qc_ref_struct.particle_index_full;
particle_frames_full = qc_ref_struct.particle_frames_full;
% iterate through snip files
exit_flag = 0;
cm = jet(128);

% index = find(all_frames==outstanding_frames(1));
index = 1;
while ~exit_flag
    % create sister_struct(i) struct    
    frame = particle_frames_full(index);
    ParticleID = particle_index_full(index);
    setID = floor(ParticleID);
    % load snip data
    load([snipPath 'pt' num2str(1e4*ParticleID) '_frame' sprintf('%03d',frame) '.mat']);
    
    %%%%%%%%%%%%%%%%%%%%%%% load image stack %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    cc = '';
    while ~strcmp(cc,'0')&&~strcmp(cc,'1')      
        edge_dist_snip = qc_spot.edge_dist_snip;
        edge_dist_rescaled = 1 + ceil(edge_dist_snip/max(edge_dist_snip(:)) * 63);
        cm = [[1 1 1] ; cm(2:end,:)];
        edge_dist_rgb = ind2rgb(edge_dist_rescaled,cm);
               
        % get frame center
        x_origin = qc_spot.x_origin;
        y_origin = qc_spot.y_origin;
        
        qc_fig = figure(1);%('Position',[0 0 512 512]);                 
%         subplot(1,2,1)
        imshow(imadjust(mat2gray(qc_spot.mcp_snip)),'InitialMagnification','fit');                        
        hold on
        p = imshow(edge_dist_rgb);        
        p.AlphaData = .4;        
        s1 = scatter(qc_spot.xp-x_origin,qc_spot.yp-y_origin,40,'MarkerFaceColor',cm(30,:),'MarkerEdgeColor','black');
        s2 = scatter(qc_spot.xc_edge-x_origin,qc_spot.yc_edge-y_origin,30,'MarkerFaceColor',cm(90,:),'MarkerEdgeAlpha',0);        
        s3 = scatter(qc_spot.xp_sister-x_origin,qc_spot.yp_sister-y_origin,40,'d','MarkerFaceColor',cm(60,:),'MarkerEdgeColor','black');        
        s4 = scatter(qc_spot.xc_serial-x_origin,qc_spot.yc_serial-y_origin,30,'MarkerFaceColor',cm(120,:),'MarkerEdgeAlpha',0);
        legend([s1 s2 s3 s4], 'spot', 'edge control', 'sister spot', 'serialized')          
        title(['Particle: ' num2str(ParticleID) ' Frame: ' sprintf('%03d',frame)])
        
        
        set(gcf,'Name',['Particle ' num2str(ParticleID) ' Frame ' num2str(frame) ' (' num2str(index) ' of ' num2str(numel(particle_frames_full)) ')'])      
        ct=waitforbuttonpress;
        cc=get(qc_fig,'currentcharacter');
        if strcmp(cc,'1')||strcmp(cc,'0')                       
            nucleus_struct(qc_spot.nc_index).qc_review_vec(qc_spot.nc_sub_index) = eval(cc);
            index = min(numel(all_frames),index + 1);
        elseif strcmp(cc,'x')
            exit_flag = 1;
            break        
        elseif strcmp(cc,'n')
            index = max(1,index-1);
            break
        elseif strcmp(cc,'m')
            index = min(numel(particle_index_full),index+1);
            break
        elseif strcmp(cc,'j')
            index = input('enter desired index: ');
            break
        elseif strcmp(cc,'s')
            disp('saving figure')
            saveas(qc_fig,[figPath 'qc_pt' num2str(1e4*ParticleID) '_frame' sprintf('%03d',frame) '.png'])
            saveas(qc_fig,[figPath 'qc_pt' num2str(1e4*ParticleID) '_frame' sprintf('%03d',frame) '.pdf'])             
        end       
    end 
%     close all
    if exit_flag
        disp('Exiting')
        close all
        break
    end
end