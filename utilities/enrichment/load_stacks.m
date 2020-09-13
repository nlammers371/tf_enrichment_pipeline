function stack = load_stacks(preFolder, frame, channel,xDim,yDim,zDim)

%   stack = [];%NaN(size(mcp_stack));
  file = dir([preFolder '/*_' sprintf('%03d',frame) '*_ch0' num2str(channel) '.tif']);
%     for im = 2:numel(files)-1      
%         stack(:,:,im-1) = double(imread([RawPath src '/' files(im).name]));
%     end
  stack = imreadStack2([file(1).folder '/' file(1).name], yDim, xDim, zDim+2);
  stack = stack(:,:,2:end-1);