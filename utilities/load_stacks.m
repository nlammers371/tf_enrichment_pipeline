function [mcp_stack, protein_stack] = load_stacks(RawPath, src, frame, mcp_channel)
    
    options = [1 2];
    protein_channel = options(options~=mcp_channel);
    % MCP
    mcp_files = dir([RawPath src '/*_' sprintf('%03d',frame) '*_ch0' num2str(mcp_channel) '.tif']);
    mcp_stack = [];
    % Protein
    protein_stack = [];%NaN(size(mcp_stack));
    protein_files = dir([RawPath src '/*_' sprintf('%03d',frame) '*_ch0' num2str(protein_channel) '.tif']);
    for im = 2:numel(mcp_files)-1
        mcp_stack(:,:,im-1) = double(imread([RawPath src '/' mcp_files(im).name]));
        protein_stack(:,:,im-1) = double(imread([RawPath src '/' protein_files(im).name]));
    end
    