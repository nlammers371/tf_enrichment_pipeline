function [imgout]=ci_loadLif(filename, getonlynumberofelements, number)
% ci_loadLif - Load Leica Image File Format
%
% Based on HKLoadlif by Hiroshi Kawaguchi
% http://de.mathworks.com/matlabcentral/fileexchange/36005-leica-image-format-file-loader
%
% [imgout]=ci_loadLif('filename', getonlynumberofelements, number) loads an image (series) from a Leica Image File
% imgout =
%               Image: {[1024x1024x101 uint8]}
%                Info: [1x1 struct]
%                Name: 'loc1'
%                Type: 'X-Y-Z'
%     NumberOfChannel: 1
%                Size: '1024  1024   101'
%
% diplibtensor=dip_image(imgout.Image); gives dipimage Tensor Image
% diplibimage=array2im(diplibtensor); gives nCh image
% channelimage=diplibimage(:,:,:,ch); gives ch image
%
% History
% Version 1.0
% Some LAS-AF 4.x files give errors, fixed in LAS-X
% (c) Ron Hoebe
% Cellular Imaging - Core Facility
% AMC - UvA - Amsterdam - The Netherlands
%% Main function
if nargin==0
    [filename, pathname]=uigetfile({'*.lif','Leica Image Format (*.lif)'});
    if filename==0; return; end
    filename=[pathname filename];
    getonlynumberofelements=false;
    number=1;
end
fp=fopen(filename,'r');
% Reading XML Part
[fp, xmlHdrStr] = ReadXMLPart(fp);
% xmlList is cell array (n x 5)
% rank(double) name(string) attributes(cell(n,2)) parant(double) children(double array)
% Changing XML to Cell
xmlList=XMLtxt2cell(xmlHdrStr);
lifVersion = GetLifVersion(xmlList(1,:)); % lifVersion is double scalar
% Reading Image Info
imgList = GetImageDescriptionList(xmlList);% imgList is struct vector
% memoryList is cell array (n x 4)
% ID(string), startPoint(uint64), sizeOfMemory, Index(double)
imgList  = ReadObjectMemoryBlocks(fp,lifVersion,imgList);
fclose(fp);
if exist('getonlynumberofelements', 'var')
    if getonlynumberofelements
        imgout=numel(imgList);
        return
    end
end
   
dat=cell(numel(imgList),5);
for n=1:numel(imgList)
    dat{n,1}=imgList(n).Name;           % Name of image
    dat{n,2}=numel(imgList(n).Channels);% number of channel
    [dimType, dimSize]=GetDimensionInfo(imgList(n).Dimensions);
    dat{n,3}=dimType;
    dat{n,4}=int2str(dimSize');
    dat{n,5}=false;
end
if exist('number', 'var')
    if number<=numel(imgList)
        n=number;
    elseif number<1
        error('Series number too low')
    else
        error('Series number too high')
    end
else
    n=1;
end
imgData  = ReadAnImageData(imgList(n), filename);
imgStruct= ReconstructImage(imgList(n), imgData);
imgStruct.Name            = dat{n,1};
imgStruct.Type            = dat{n,3};
imgStruct.NumberOfChannel = dat{n,2};
imgStruct.Size            = dat{n,4};
imgout=imgStruct;
return
function [dimType, dimSize]=GetDimensionInfo(dimensions)
ndims=numel(dimensions);
dimType=cell(2*ndims,1);
dimSize=zeros(ndims,1);
for n=1:ndims
    dimSize(n)=str2double(dimensions(n).NumberOfElements);
    switch dimensions(n).DimID
        case '0';dimType{2*n-1}='Not Valied';
        case '1';dimType{2*n-1}='X';
        case '2';dimType{2*n-1}='Y';
        case '3';dimType{2*n-1}='Z';
        case '4';dimType{2*n-1}='T';
        case '5';dimType{2*n-1}='Lambda';
        case '6';dimType{2*n-1}='Rotation';
        case '7';dimType{2*n-1}='XT';
        case '8';dimType{2*n-1}='TSlice';
        otherwise
    end
    dimType{2*n}='-';
end
dimType=strrep(sprintf('%s',char(dimType)'),' ','');
dimType=dimType(1:end-1);% Delete last '-'
function imgList=ReconstructImage(imgInfo, imgData)
% ===================================================================
% imgInfo description
% ===================================================================
% <ChannelDescription>
%%% DataType   [0, 1]               [Integer, Float]
%%% ChannelTag [0, 1, 2, 3]         [GrayValue, Red, Green, Blue]
% Resolution [Unsigned integer]   Bits per pixel if DataType is Float value can be 32 or 64 (float or double)
% NameOfMeasuredQuantity [String] Name
% Min        [Double] Physical Value of the lowest gray value (0). If DataType is Float the Minimal possible value (or 0).
% Max        [Double] Physical Value of the highest gray value (e.g. 255) If DataType is Float the Maximal possible value (or 0).
% Unit       [String] Physical Unit
% LUTName    [String] Name of the Look Up Table (Gray value to RGB value)
% IsLUTInverted [0, 1] Normal LUT Inverted Order
%%% BytesInc   [Unsigned long (64 Bit)] Distance from the first channel in Bytes
% BitInc     [Unsigned Integer]       Bit Distance for some RGB Formats (not used in LAS AF 1..0 ? 1.7)
% <DimensionDescription>
% DimID   [0, 1, 2, 3, 4, 5, 6, 7, 8] [Not valid, X, Y, Z, T, Lambda, Rotation, XT Slices, T Slices]
%%% NumberOfElements [Unsigned Integer] Number of elements in this dimension
% Origin           [Unsigned integer] Physical position of the first element (Left pixel side)
% Length   [String] Physical Length from the first left pixel side to the last left pixel side (Not the right. A Pixel has no width!)
% Unit     [String] Physical Unit
%%% BytesInc [Unsigned long (64 Bit)] Distance from one Element to the next in this dimension
% BitInc   [Unsigned Integer] Bit Distance for some RGB Formats (not used, i.e.: = 0 in LAS AF 1..0 ? 1.7)
% ===================================================================
% imgList info
% ===================================================================
% img
% Get Dimension info
dimension=ones(1,9);
for m=1:numel(imgInfo.Dimensions)
    dimension(str2double(imgInfo.Dimensions(m).DimID))=str2double(imgInfo.Dimensions(m).NumberOfElements);
end
% Separate to each channel image
nCh=numel(imgInfo.Channels);
imgList=struct('Image',[],'Info',[]);
imgList.Image = cell(nCh,1);
if nCh > 1;
    imgData=reshape(imgData,...
        str2double(imgInfo.Channels(2).BytesInc)-str2double(imgInfo.Channels(1).BytesInc),[]);
    for m=1:nCh
        tmp=imgData(:,m:nCh:end);
        imgList.Image{m}= reshape(typecast(tmp(:),GetType(imgInfo.Channels(m))),dimension);
    end
else
    imgList.Image{1}=reshape(typecast(imgData,GetType(imgInfo.Channels)),dimension);
end
imgList.Info = imgInfo;
function chType=GetType(Channels)
switch str2double(Channels.DataType)
    case 0 % int case
        switch str2double(Channels.Resolution);
            % currently, resolution is constant through the channels
            case 8;   chType='uint8';
            case 12;  chType='uint16';
            case 32;  chType='uint32';
            case 64;  chType='uint64';
            otherwise;error('Unsupported data bit. ')
        end
    case 1 % float case
        switch str2double(Channels.Resolution);
            % currently, resolution is constant through the channels
            case 32;  chType='single';
            case 64;  chType='double';
            otherwise;error('Unsupported data bit. ')
        end
end
function mems = GetImageDescriptionList(xmlList)
%  For the image data type the description of the memory layout is defined
%  in the image description XML node (<ImageDescription>).
% <ImageDescription>
imgIndex  =SearchTag(xmlList,'ImageDescription');
numImgs=numel(imgIndex);
% <Memory Size="21495808" MemoryBlockID="MemBlock_233"/>
memIndex  =SearchTag(xmlList,'Memory');
memSizes  =cellfun(@str2double,GetAttributeVal(xmlList,memIndex,'Size'));
memIndex=memIndex(memSizes~=0);
memSizes=memSizes(memSizes~=0);
if numImgs~=numel(memIndex);
    error('Number of ImageDescription and Memory did not match.')
end
% Matching ImageDescription with Memory
imgParentElmIndex =  zeros(numImgs,1);
for n=1:numImgs
    imgParentElmIndex(n) = SearchTagParent(xmlList,imgIndex(n),'Element');
end
memParentElmIndex =  zeros(numImgs,1);
for n=1:numImgs
    memParentElmIndex(n) = SearchTagParent(xmlList,memIndex(n),'Element');
end
[imgParentElmIndex, sortIndex]=sort(imgParentElmIndex); imgIndex=imgIndex(sortIndex);
[memParentElmIndex, sortIndex]=sort(memParentElmIndex); memIndex=memIndex(sortIndex);memSizes=memSizes(sortIndex);
if ~all(imgParentElmIndex==memParentElmIndex)
    error('Matching ImageDescriptions with Memories')
end
mems=struct('Name',[],'Channels',[],'Dimensions',[],'Memory',[]);
mems(numImgs).Memory.StartPosition=[];
for n=1:numImgs
    mems(n).Name = char(GetAttributeVal(xmlList, imgParentElmIndex(n),'Name'));
    [mems(n).Channels, mems(n).Dimensions]= MakeImageStruct(xmlList,imgIndex(n));
    mems(n).Memory.Size=memSizes(n);
    mems(n).Memory.MemoryBlockID=char(GetAttributeVal(xmlList,memIndex(n),'MemoryBlockID'));
end
return
function [C, D]=MakeImageStruct(xmlList,iid)
% ChannelDescription   DataType="0" ChannelTag="0" Resolution="8"
%                      NameOfMeasuredQuantity="" Min="0.000000e+000" Max="2.550000e+002"
%                      Unit="" LUTName="Red" IsLUTInverted="0" BytesInc="0"
%                      BitInc="0"
% DimensionDescription DimID="1" NumberOfElements="512" Origin="4.336809e-020"
%                      Length="4.558820e-004" Unit="m" BitInc="0"
%                      BytesInc="1"
% Memory ?@?@?@?@?@?@?@  Size="21495808" MemoryBlockID="MemBlock_233"
iidChildren=xmlList{iid,5};
for n=1:numel(iidChildren)
    if strcmp(xmlList{iidChildren(n),2},'Channels')
        id=xmlList{iidChildren(n),5};
        p=xmlList(id,3);
        nid=numel(id);
        tmp=cell(11,nid);
        for m=1:nid
            tmp(:,m)=p{m}(:,2);
        end
        C=cell2struct(tmp,p{1}(:,1),1);
    elseif strcmp(xmlList{iidChildren(n),2},'Dimensions')
        id=xmlList{iidChildren(n),5};
        p=xmlList(id,3);
        nid=numel(id);
        tmp=cell(7,nid);
        for m=1:nid
            tmp(:,m)=p{m}(:,2);
        end
        D=cell2struct(tmp,p{1}(:,1),1);
    else
        error('Undefined Tag')
    end
end
function lifVersion = GetLifVersion(xmlList)
% return version of header
index  =SearchTag(xmlList,'LMSDataContainerHeader');
value  =GetAttributeVal(xmlList,index,'Version');
lifVersion = str2double(cell2mat(value(1)));
return
function pindex=SearchTagParent(xmlList,index,tagName)
% return the row index of given tag name
pindex=xmlList{index,4};
while pindex~=0
    if strcmp(xmlList{pindex,2},tagName)
        return;
    else
        pindex=xmlList{pindex,4};
    end
end
error('Cannot Find the Parent Tag "%s"',tagName);
return;
function index=SearchTag(xmlList,tagName)
% return the row index of given tag name
listLen=size(xmlList,1);
index=[];
for n=1:listLen
    if strcmp(char(xmlList(n,2)),tagName)
        index=[index; n]; %#ok<AGROW>
    end
end
function value=GetAttributeVal(xmlList, index, attributeName)
% return cell array of attributes row index of given tag name
value={};
for n=1:length(index)
    currentCell=xmlList{index(n),3};
    for m=1:size(currentCell,1)
        if strcmp(char(currentCell(m,1)),attributeName)
            value=[value; currentCell(m,2)]; %#ok<AGROW>
        end
    end
end
function CheckTestValue(value,errorMsg)
switch class(value)
    case 'uint8';  trueVal=hex2dec('2A');
    case 'uint32'; trueVal=hex2dec('70');
    otherwise;
        error('Unsupported Error Number: %d',value)
end
if value~=trueVal
    error(errorMsg);
end
return;
function [fp, str, ketPos] = ReadXMLPart(fp)
% Size(bytes) Total(bytes) description (some LAS-AF 3.x version give errors
% here, fixed in LAS-X
CheckTestValue(fread(fp,1,'*uint32'),...        % 4  4 Test Value 0x70
    'Invalid test value at Part: XML.');
xmlChunk = fread(fp, 1, 'uint32');              % 4  8 Binary Chunk length NC*2 + 1 + 4
CheckTestValue(fread(fp,1,'*uint8'),...         % 1  9 Test Value 0x2A
    'Invalid test value at XML Content.');
nc = fread(fp,1,'uint32');                      % 4 13 Number of UTF-16 Characters (NC)
if (nc*2 + 1 + 4)~=xmlChunk;  %
    error('Chunk size mismatch at Part: XML.');
end
str= fread(fp,nc*2,'char');                     % 2*nc - XML Object Description
% UTF-16 -> UTF-8 (cut zeros)
str     = char(str(1:2:end)');
% Insert linefeed(char(10)) for facilitate visualization -----
% str=strrep(str,'><',['>' char(10) '<']);
ketPos =strfind(str,'>'); % find position of ">" for fast search of element
return;
function imgLists=ReadObjectMemoryBlocks(fp,lifVersion,imgLists)
% get end of file and return current point
cofp=    ftell(fp);
fseek(fp,0,'eof');
eofp=    ftell(fp);
fseek(fp,cofp,'bof');
nImgLists=length(imgLists);
memoryList=cell(nImgLists,4);
% ID(string), startPoint(uint64), sizeOfMemory, Index(double)
for n = 1:nImgLists
    memoryList{n,1}=imgLists(n).Memory.MemoryBlockID;
end
% read object memory blocks
while ftell(fp) < eofp;
    
    CheckTestValue(fread(fp,1,'*uint32'),...        % Test Value 0x70
        'Invalied test value at Object Memory Block');
    
    objMemBlkChunk = fread(fp, 1, '*uint32');%#ok<NASGU> % Size of Description
    
    CheckTestValue(fread(fp,1,'*uint8'),...         % Test Value 0x2A
        'Invalied test value at Object Memory Block');
    
    
    switch uint8(lifVersion)            % Size of Memory (version dependent)
        case 1; sizeOfMemory = double(fread(fp, 1, '*uint32'));
        case 2; sizeOfMemory = double(fread(fp, 1, '*uint64'));
        otherwise; error('Unsupported LIF version. Update this program');
    end
    
    CheckTestValue(fread(fp,1,'*uint8'),...         % Test Value 0x2A
        'Invalied test value at Object Memory Block');
    
    nc = fread(fp,1,'*uint32');                     % Number of MemoryID string
    
    str = fread(fp,nc*2,'*char')';                  % Number of MemoryID string (UTF-16)
    str = char(str(1:2:end));                       % convert UTF-16 to UTF-8
    
    if sizeOfMemory > 0
        for n=1:nImgLists
            if strcmp(char(memoryList{n,1}),str) % NEED CONSIDERATION !!!!!!
                if imgLists(n).Memory.Size ~= sizeOfMemory;
                    error('Memory Size Mismatch.');
                end
                imgLists(n).Memory.StartPosition=ftell(fp);
                fseek(fp,sizeOfMemory,'cof');
                break;
            end
        end
    end
end
return;
function imgData=ReadAnImageData(imgInfo,fileName)
fp = fopen(fileName,'rb');
if fp<0; errordlg('Cannot open file: \n\t%s', fileName); end
fseek(fp,imgInfo.Memory.StartPosition,'bof');
imgData = fread(fp,imgInfo.Memory.Size,'*uint8');
fclose(fp);
function tagList=XMLtxt2cell(c)
% rank(double) name(string) attributes(cell(n,2)) parant(double) children(double array)
tags  =regexp(c,'<("[^"]*"|''[^'']*''|[^''>])*>','match')';
nTags=numel(tags);
tagList=cell(nTags,5);
tagRank=0;
tagCount=0;
for n=1:nTags
    currentTag=tags{n}(2:end-1);
    if currentTag(1)=='/'
        tagRank=tagRank-1;
        continue;
    end
    tagRank=tagRank+1;
    tagCount=tagCount+1;
    [tagName, attributes]=ParseTagString(currentTag);
    tagList{tagCount,1}=tagRank;
    tagList{tagCount,2}=tagName;
    tagList{tagCount,3}=attributes;
    % search parant
    if tagRank~=1
        if tagRank~=tagList{tagCount-1,1};
            tagRankList=cell2mat(tagList(1:tagCount,1));
            parent=find(tagRankList==tagRank-1,1,'last');
            tagList{tagCount,4}=parent;
        else
            tagList{tagCount,4}=tagList{tagCount-1,4};
        end
    else
        tagList{tagCount,4} = 0;
    end
    if currentTag(end)=='/'
        tagRank=tagRank-1;
    end
end
tagList   =tagList(1:tagCount,:);
parentList=cell2mat(tagList(:,4));
% Make Children List
for n=1:tagCount
    tagList{n,5}=find(parentList==n);
end
return;
function [name, attributes]=ParseTagString(tag)
[name, tmpAttributes]=regexp(tag,'^\w+','match', 'split');
name=char(name);
attributesCell=regexp(char(tmpAttributes(end)),'\w+=".*?"','match');
if isempty(attributesCell)
    attributes={};
else
    nAttributes = numel(attributesCell);
    attributes=cell(nAttributes,2);
    for n=1:nAttributes
        currAttrib=char(attributesCell(n));
        dqpos=strfind(currAttrib,'"');
        attributes{n,1}=currAttrib(1:dqpos(1)-2);
        if dqpos(2)-dqpos(1)==1 % case attribute=""
            attributes{n,2}='';
        else
            attributes{n,2}=currAttrib(dqpos(1)+1:dqpos(2)-1);
        end
    end
end
