function stack = load_stacks(RawPath, src, frame, channel)

    stack = [];%NaN(size(mcp_stack));
    files = dir([RawPath src '/*_' sprintf('%03d',frame) '*_ch0' num2str(channel) '.tif']);
    for im = 2:numel(files)-1      
        stack(:,:,im-1) = double(imread([RawPath src '/' files(im).name]));
    end
    