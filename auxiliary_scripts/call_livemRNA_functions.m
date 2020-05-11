% clear
close all
% add path to mRNADynamics folder
livemRNAPath = 'E:\Nick\LivemRNA\';
addpath(genpath([livemRNAPath '\mRNADynamics']))
cd(livemRNAPath)
% addpath('D:\Nick\LivemRNA\ImportFromAWS.m')
% addpath('D:\Nick\LivemRNA\ComputerFolders.csv')
% set path to dropbox
sheet_path = 'S:\Nick\Dropbox\LocalEnrichmentResults\DataStatus_NLtemp.xlsx';
DropboxTab = '2xDl-Ven_snaBAC-mCh';
% specfiy functions to fun
% command_cell = {'ImportFromAWS(Prefix)'};
% command_cell = {'segmentSpots(Prefix,5000,"Weka","fit3DOnly")'};...,'TrackmRNADynamics(Prefix,5000,5000)',...
%     'AddParticlePosition(Prefix)','CompileParticles(Prefix,"SkipAll","ApproveAll")'};
command_cell = {'TrackNuclei(Prefix)','TrackmRNADynamics(Prefix,5000,5000)'};
% command_cell = {'AddParticlePosition(Prefix)','CompileParticles(Prefix)'};
% command_cell = {'CompileParticles(Prefix,"SkipAll","ApproveAll")'};

%%%%%%%%%%%%% get prefix list %%%%%%%%%%%%%%%%
% find sheet
[~ ,sheet_names] = xlsfinfo(sheet_path);
sheet_index = find(ismember(sheet_names,DropboxTab));
[~,~,sheet_cell] = xlsread(sheet_path,sheet_index);
name_col = sheet_cell(1:33,1); % hard coded for now
ready_ft = contains(name_col,'Ignore');
ready_cols = 1 + find([sheet_cell{ready_ft,2:end}]==0);
sheet_cell = sheet_cell(:,[1 ready_cols]);
% get list of project names
prefix_ft = contains(name_col,'Prefix');
prefix_cell_raw = sheet_cell(prefix_ft,2:end);
prefix_cell = {};
for i = 1:numel(prefix_cell_raw)
    if ~isempty(prefix_cell_raw{i})
        eval([prefix_cell_raw{i} ';'])
        prefix_cell = [prefix_cell{:} {Prefix}];
    end
end

error_array = false(numel(prefix_cell),numel(command_cell));
%%%%%%%%%% execute command list %%%%%%%%%%%%%%
for c = 1:numel(command_cell)
    cmd_string = command_cell{c};
    for p = 1:numel(prefix_cell)
        Prefix = prefix_cell{p};
        disp(['executing ' cmd_string '(' num2str(p) '/' num2str(numel(prefix_cell)) ')'])
        tic
        try            
            eval(cmd_string);
        catch
            error_array(p,c) = true;
        end
        toc
    end
end
        