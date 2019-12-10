classdef MECTracker
    %MECHCYGROWTHMAIN  Main class for the mechanical effects on growth
    %
    %  This is the main class for running the tracking and analysis for the
    %  mechanical effects on cyanobacteria growth project.
    %
    %  For full documentation, see the project wiki: 
    %    https://biof-git.colorado.edu/cameron-lab/crowding-effects-on-cyanobacteria-growth/wikis/home
    %
    %  Copyright 2017 CU Boulder and the Cameron Lab
       
    
    properties
        
        FrameRange = Inf;        
        SeriesRange = Inf;
        OutputMovie = true;
        
        %Segmentation options
        ChannelToSegment = '!CStack';
        CellMarkerChannel = '';
        ThresholdLevel = 0.05;
                     
        %Track linking parameters
        LinkedBy = 'Centroid';
        LinkCalculation = 'euclidean';
        LinkingScoreRange = [-Inf, Inf];
        
        MaxTrackAge = 2;
        
        HighIntPrctile = 95;    %Prctile value for the high intensity value
        
        %Mitosis detection parameters
        TrackMitosis = true;
        MinAgeSinceMitosis = 2;
        MitosisParameter = 'PixelIdxList';          %What property is used for mitosis detection?
        MitosisCalculation = 'pxintersect';
        MitosisScoreRange = [-Inf, Inf];
        MitosisLinkToFrame = -1;                    %What frame to link to/ This should be 0 for centroid/nearest neighbor or -1 for overlap (e.g. check with mother cell)
        
        %LAP solver
        LAPSolver = 'lapjv';
        
        EnableParallel = false;
                
    end
    
    methods
        
        function trackLinker = processFiles(obj, varargin)
            %PROCESSFILE  Segment and track cells in a video file
            
            ip = inputParser;
            ip.addOptional('Filename','',@(x) ischar(x));
            ip.addParameter('OutputDir',0, @(x) exist(x,'dir'));
            ip.parse(varargin{:});
            
            %Validate the inputs
            if isempty(ip.Results.Filename)
                [fname, fpath] = uigetfile({'*.nd2','ND2 file (*.nd2)'},...
                    'Select a file');
                
                if fname == 0
                    %Stop running the script
                    return;
                end
                
                filename = fullfile(fpath,fname);
                
            else
                filename = ip.Results.Filename;
                
            end
            
            %Prompt user for a directory to save output files to
            if ip.Results.OutputDir == 0
                
                startPath = fileparts(filename);
                
                outputDir = uigetdir(startPath, 'Select output directory');
                
                if outputDir == 0
                    %Processing cancelled
                    return;
                end
                
            else
                
                outputDir = ip.Results.Filename;
                
            end
            
            %Get a reader object for the image
            bfReader = BioformatsImage(filename);
            
            %Get the frame range to process
            if isinf(obj.FrameRange)
                frameRange = 1:bfReader.sizeT;
            else
                %TODO: Check that requested range is actually within sizeT
                frameRange = obj.FrameRange;
            end
            
            if isinf(obj.SeriesRange)
                seriesRange = 1:bfReader.seriesCount;
            else
                seriesRange = obj.SeriesRange;                
            end
            
            for iSeries = seriesRange
                
                %Set the image series number
                bfReader.series = iSeries;
                
                %Segment cells
                for iT = frameRange
                    
                    %Get image to segment
                    if strcmpi(obj.ChannelToSegment,'!CStack')
                        
                        imgToSegment = zeros(bfReader.height, bfReader.width);
                        for iC = 1:bfReader.sizeC
                            
                            %Add the intensity of a channel, weighted by the
                            %number of channels
                            imgToSegment = imgToSegment + double(bfReader.getPlane(1, iC, iT));% ./ bfReader.sizeC;
                            
                        end
                        
                    else
                        imgToSegment = bfReader.getPlane(1, obj.ChannelToSegment, iT);
                    end
                    
                    if ~isempty(obj.CellMarkerChannel)
                        markerImg = bfReader.getPlane(1,obj.CellMarkerChannel,iT);
                    
                        %Segment cells and get cell label
                        cellLabels = MECTracker.getCellLabels(imgToSegment, markerImg);
                    else
                        cellLabels = MECTracker.getCellLabels(imgToSegment);
                    end
                    
                    %Get cell data
                    cellData = MECTracker.getCellData(cellLabels, bfReader, iT, obj.HighIntPrctile);
                    
                    if numel(cellData) == 0
                        warning('No cell data');
                        continue
                    end
                                       
                    %Link cells
                    if iT == frameRange(1)
                        %Set up the cell tracker
                        trackLinker = TrackLinker(iT, cellData);
                        trackLinker = trackLinker.setOptions(obj.propsToStruct);
                        trackLinker = trackLinker.setFilename(bfReader.filename);
                        trackLinker = trackLinker.setPxSizeInfo(bfReader.pxSize(1),bfReader.pxUnits);
                        trackLinker = trackLinker.setImgSize([bfReader.height, bfReader.width]);

                    else
                        try
                            trackLinker = trackLinker.assignToTrack(iT, cellData);
                        catch
                            
                            saveData = input('There was an error linking tracks. Would you like to save the tracked data generated so far?\n');
                            if saveData
                                trackArray = trackLinker.getTrackArray; %#ok<NASGU>
                                save(fullfile(outputDir, sprintf('%s_series%d.mat',fname, iSeries)), 'trackArray');
                            end
                            clear trackArray
                            
                            keyboard
                        end
                    end
                    
                    if obj.OutputMovie
                        
                        [~, fname] = fileparts(filename);
                        
                        if ~exist('vidObj','var')
                            vidObj = VideoWriter(fullfile(outputDir,sprintf('%s_series%d.avi',fname,iSeries)));
                            vidObj.FrameRate = 10;
                            vidObj.Quality = 100;
                            open(vidObj);
                        end
                        
                        imgOut = MECTracker.makeAnnotatedImage(iT, imgToSegment, cellLabels, trackLinker);
                        vidObj.writeVideo(imgOut);
                        
                    end
                    
                end
                
                if exist('vidObj','var')
                    close(vidObj);
                    clear vidObj
                end
                
                %Add timestamp information to the track array
                trackArray = trackLinker.getTrackArray; 
               
                %Add timestamp information
                [ts, tsunit] = bfReader.getTimestamps(1,1);
                ts = ts(frameRange);
                trackArray = trackArray.setTimestampInfo(ts,tsunit); %#ok<NASGU>
                
                %Save the data
                
                save(fullfile(outputDir, sprintf('%s_series%d.mat',fname, iSeries)), 'trackArray');
                clear trackArray
            end
                                   
            %Save the settings file
            obj.exportOptions(fullfile(outputDir,'settings.txt'));
            
        end
        
        function obj = setOptions(obj, varargin)
            %SETOPTIONS  Set options for the linker
            %
            %  linkerObj = linkerObj.SETOPTIONS(parameter, value) will set
            %  the parameter to value.
            %
            %  linkerObj = linkerObj.SETOPTIONS(O) where O is a data object
            %  with the same property names as the options will work.
            %
            %  linkerObj = linkerObj.SETOPTIONS(S) where S is a struct
            %  with the same fieldnames as the options will also work.
            %
            %  Non-matching parameter names will be ignored.
            
            if numel(varargin) == 1 && isstruct(varargin{1})
                %Parse a struct as input
                
                inputParameters = fieldnames(varargin{1});
                
                for iParam = 1:numel(inputParameters)
                    if ismember(inputParameters{iParam},fieldnames(obj.options))
                        obj.options.(inputParameters{iParam}) = ...
                            varargin{1}.(inputParameters{iParam});
                    else
                        %Just skip unmatched options
                    end
                    
                end
                
            elseif numel(varargin) == 1 && isobject(varargin{1})
                %Parse an object as input
                
                inputParameters = properties(varargin{1});
                
                for iParam = 1:numel(inputParameters)
                    if ismember(inputParameters{iParam},fieldnames(obj.options))
                        obj.options.(inputParameters{iParam}) = ...
                            varargin{1}.(inputParameters{iParam});
                    else
                        %Just skip unmatched options
                    end
                    
                end
                
            else
                if rem(numel(varargin),2) ~= 0
                    error('Input must be Property/Value pairs.');
                end
                inputArgs = reshape(varargin,2,[]);
                for iArg = 1:size(inputArgs,2)
                    if ismember(inputArgs{1,iArg},properties(obj))
                        
                        obj.(inputArgs{1,iArg}) = inputArgs{2,iArg};
                    else
                        %Just skip unmatched options
                    end
                end
            end
        end
        
        function obj = importOptions(obj, importFilename)
            %IMPORTOPTIONS  Import settings from file
            %
            %  S = L.IMPORTOPTIONS(filename) will import
            %  settings from the file specified.
            
            if ~exist('importFilename','var')
                
                [filename, pathname] = uigetfile({'*.txt','Text file (*.txt)';...
                    '*.*','All files (*.*)'},...
                    'Select output file location');
                
                importFilename = fullfile(pathname,filename);
                
            end
            
            fid = fopen(importFilename,'r');
            
            if fid == -1
                error('TrackLinker:importOptions:ErrorReadingFile',...
                    'Could not open file %s for reading.',filename);
            end
            
            while ~feof(fid)
                currLine = strtrim(fgetl(fid));
                
                if isempty(currLine)
                    %Empty lines should be skipped
                    
                elseif strcmpi(currLine(1),'%') || strcmpi(currLine(1),'#')
                    %Lines starting with '%' or '#' are comments, so ignore
                    %those
                    
                else
                    %Expect the input to be PARAM_NAME = VALUE
                    parsedLine = strsplit(currLine,'=');
                    
                    %Get parameter name (removing spaces)
                    parameterName = strtrim(parsedLine{1});
                    
                    %Get value name (removing spaces)
                    value = strtrim(parsedLine{2});
                    
                    if isempty(value)
                        %If value is empty, just use the default
                    else
                        obj = obj.setOptions(parameterName,eval(value));
                    end
                    
                end
                
            end
            
            fclose(fid);
        end
        
        function exportOptions(obj, exportFilename)
            %EXPORTOPTIONS  Export tracking options to a file
            %
            %  L.EXPORTOPTIONS(filename) will write the currently set
            %  options to the file specified. The options are written in
            %  plaintext, no matter what the extension of the file is.
            %
            %  L.EXPORTOPTIONS if the filename is not provided, a dialog
            %  box will pop-up asking the user to select a location to save
            %  the file.
            
            if ~exist('exportFilename','var')
                
                [filename, pathname] = uiputfile({'*.txt','Text file (*.txt)'},...
                    'Select output file location');
                
                exportFilename = fullfile(pathname,filename);
                
            end
            
            fid = fopen(exportFilename,'w');
            
            if fid == -1
                error('FRETtrackerOptions:exportSettings:CouldNotOpenFile',...
                    'Could not open file to write')
            end
            
            propertyList = properties(obj);
            
            %Write output data depending on the datatype of the value
            for ii = 1:numel(propertyList)
                
                if ischar(obj.(propertyList{ii}))
                    fprintf(fid,'%s = ''%s'' \r\n',propertyList{ii}, ...
                        obj.(propertyList{ii}));
                    
                elseif isnumeric(obj.(propertyList{ii}))
                    fprintf(fid,'%s = %s \r\n',propertyList{ii}, ...
                        mat2str(obj.(propertyList{ii})));
                    
                elseif islogical(obj.(propertyList{ii}))
                    
                    if obj.(propertyList{ii})
                        fprintf(fid,'%s = true \r\n',propertyList{ii});
                    else
                        fprintf(fid,'%s = false \r\n',propertyList{ii});
                    end
                    
                end
                
            end
            
            fclose(fid);
            
        end
        
        function sOut = propsToStruct(obj)
            
            propList = properties(obj);
            for iP = 1:numel(propList)
                sOut.(propList{iP}) = obj.(propList{iP});
            end
        end
        
    end
    
    methods (Static)
        
        function cellLabels = getCellLabels(cellImage, varargin)
            %GETCELLLABELS  Segment and label individual cells
            %
            %  L = MECTracker.GETCELLLABELS(I) will segment the cells in image
            %  I, returning a labelled image L. Each value in L should
            %  correspond to an individual cell.
            %
            %  L = MECTracker.GETCELLLABELS(I, M) will use image M to mark
            %  cells. M should be a fluroescent image (e.g. YFP, GFP) that
            %  fully fills the cells.
            
            if isempty(varargin)

                cellMarkerImg = [];
                
            else
                
                cellMarkerImg = varargin{1};
                
            end
            

            %Normalize the cellImage
            cellImage = MECTracker.normalizeimg(cellImage);
            
            if isinteger(cellImage)
                %Get an initial mask and clean it up
                mask = cellImage > 0.05 * 65535;
            else
                mask = cellImage > 0.05;
            end

            mask = imopen(mask,strel('disk',2));
            mask = imclearborder(mask);
            
            mask = activecontour(cellImage,mask);
            
            mask = bwareaopen(mask,100);
            mask = imopen(mask,strel('disk',2));
            
            mask = imfill(mask,'holes');
                     
%             MECTracker.showoverlay(MECTracker.normalizeimg(cellImage),bwperim(mask),[0 1 0]);
%             keyboard

%---
%---
            
            if isempty(cellMarkerImg)
                
                dd = -bwdist(~mask);
                dd(~mask) = -Inf;
                
                imgToWatershed = imhmin(dd,1);
                
            else
                %TODO
                
                %%Find extended maxima in the cell marker channel
                %fgmarker = ColonyTracker.getCellMarker(cellMarkerImg);
                
                %%Invert the cell marker image so that the centers are dark
                %img = imcomplement(cellMarkerImg);
                
                %Mark the negative regions (e.g. regions that are not the cell
                %and each cell region)
                %markedImage = imimposemin(img, ~colonyMask | fgmarker);
            
                cellMarkerImgFilt = medfilt2(cellMarkerImg,[10 10]);
                cellMarkerMask = imextendedmax(cellMarkerImgFilt, 50);
                
                %Invert the cell marker image so that the centers are dark
                img = imcomplement(cellMarkerImg);
                
                %Mark the negative regions (e.g. regions that are not the cell
                %and each cell region)
                imgToWatershed = imimposemin(img, ~mask | cellMarkerMask);
                
            end
            
            cellLabels = watershed(imgToWatershed);
            cellLabels = imclearborder(cellLabels);
            
        end
        
        function thLvl = getThreshold(imageIn)
            %GETTHRESHOLD  Get a threshold level to binarize the image
            %
            %  T = MECTracker.GETTHRESHOLD(I) will look for a suitable
            %  greyscale level T to binarize image I.
%             
            [nCnts, binEdges] = histcounts(imageIn(:),150);
            binCenters = diff(binEdges) + binEdges(1:end-1);
            
            %nCnts = smooth(nCnts,3);
%            nCnts(1) = 0;
            %[~,locs] = findpeaks(nCnts,'Npeaks',2,'sortStr','descend','MinPeakDistance',5);
                        
            %Determine the background intensity level
            [~,locs] = findpeaks(nCnts,'Npeaks',1,'sortStr','descend');
            
            thLvl = find(nCnts(locs:end) <= 0.2 * nCnts(locs));
            
%             %Find valley
%             [~,valleyLoc] = min(nCnts(locs(1):locs(2)));
%             
%             thLvl = binCenters(valleyLoc + locs(1));
            
            if isempty(thLvl)
                warning('Threshold level not found')
                keyboard
            end
            
        end
        
        function varargout = showoverlay(baseimage, mask, color, varargin)
            %SHOWOVERLAY    Plot an overlay mask on an image
            %
            %  SHOWOVERLAY(IMAGE,MASK,COLOR) will plot an overlay specified by a binary
            %  MASK on the IMAGE. The color of the overlay is specified using a three
            %  element vector COLOR.
            %
            %  Example:
            %
            %    mainImg = imread('cameraman')
            %
            %
            %  Downloaded from http://cellmicroscopy.wordpress.com
            
            if ~exist('color','var')
                color = [1 1 1]; %Default color of the overlay
            end
            
            if size(baseimage,3) == 3
                red = baseimage(:,:,1);
                green = baseimage(:,:,2);
                blue = baseimage(:,:,3);
                
            elseif size(baseimage,3) == 1
                red = baseimage;
                green = baseimage;
                blue = baseimage;
                
            else
                error('Image should be either NxNx1 (greyscale) or NxNx3 (rgb)')
            end
            
            %Make sure the mask is binary (anything non-zero becomes true)
            mask = (mask ~= 0);
            
            if isinteger(baseimage)
                maxInt = intmax(class(baseimage));
            else
                maxInt = 1;
            end
            
            red(mask) = color(1) .* maxInt;
            green(mask) = color(2) .* maxInt;
            blue(mask) = color(3) .* maxInt;
            
            %Concatenate the output
            outputImg = cat(3,red,green,blue);
            
            if nargout == 0
                %Get the current warning status
                warningStatus = warning;
                warningStatus = warningStatus(1).state;
                
                %Turn off warnings
                warning off
                
                %Show image
                imshow(outputImg,[])    
                
                %Restore warning state
                warning(warningStatus)
            else
                varargout{1} = outputImg;
            end
        end
        
        function imageOut = normalizeimg(imageIn,varargin)
            %NORMALIZEIMG   Linear dynamic range expansion for contrast enhancement
            %   N = NORMALIZEIMG(I) expands the dynamic range (or contrast) of image I
            %   linearly to maximize the range of values within the image.
            %
            %   This operation is useful when enhancing the contrast of an image. For
            %   example, if I is an image with uint8 format, with values ranging from
            %   30 to 100. Normalizing the image will expand the values so that they
            %   fill the full dynamic range of the format, i.e. from 0 to 255.
            %
            %   The format of the output image N depends on the format of the input
            %   image I. If I is a matrix with an integer classs (i.e. uint8, int16), N
            %   will returned in the same format. If I is a double, N will be
            %   normalized to the range [0 1] by default.
            %
            %   N = NORMALIZEIMG(I,[min max]) can also be used to specify a desired
            %   output range. For example, N = normalizeimg(I,[10,20]) will normalize
            %   image I to have values between 10 and 20. In this case, N will be
            %   returned in double format regardless of the format of I.
            %
            %   In situations where most of the interesting image features are
            %   contained within a narrower band of values, it could be useful to
            %   normalize the image to the 5 and 95 percentile values.
            %
            %   Example:
            %       I = imread('cameraman.tif');
            %
            %       %Calculate the values corresponding to the 5 and 95 percentile of
            %       %values within the image
            %       PRC5 = prctile(I(:),5);
            %       PRC95 = prctile(I(:),95);
            %
            %       %Threshold the image values to the 5 and 95 percentiles
            %       I(I<PRC5) = PRC5;
            %       I(I>PRC95) = PRC95;
            %
            %       %Normalize the image
            %       N = normalizeimg(I);%
            %
            %       %Display the normalized image
            %       imshow(N)
            
            %Define default output value range
            outputMin = 0;
            outputMax = 1;
            
            %Check if the desired output range is set. If it is, make sure it contains
            %the right number of values and format, then update the output minimum and
            %maximum values accordingly.
            if nargin >= 2
                if numel(varargin{1}) ~= 2
                    error('The input parameter should be [min max]')
                end
                
                outputMin = varargin{1}(1);
                outputMax = varargin{1}(2);
            else
                %If the desired output range is not set, then check if the image is an
                %integer class. If it is, then set the minimum and maximum values
                %to match the range of the class type.
                if isinteger(imageIn)
                    inputClass = class(imageIn);
                    
                    outputMin = 0;
                    outputMax = double(intmax(inputClass)); %Get the maximum value of the class
                    
                end
            end
            
            %Convert the image to double for the following operations
            imageIn = double(imageIn);
            
            %Calculate the output range
            outputRange = outputMax - outputMin;
            
            %Get the maximum and minimum input values from the image
            inputMin = min(imageIn(:));
            inputMax = max(imageIn(:));
            inputRange = inputMax - inputMin;
            
            %Normalize the image values to fit within the desired output range
            imageOut = (imageIn - inputMin) .* (outputRange/inputRange) + outputMin;
            
            %If the input was an integer before, make the output image the same class
            %type
            if exist('inputClass','var')
                eval(['imageOut = ',inputClass,'(imageOut);']);
            end
            
        end
        
        function cellData = getCellData(cellLabels, bfReader, iT, hiIntLevel)
            %GETCELLDATA  Get cell data
            
            %Let's get basic data for now
            
            %Get standard data
            cellData = regionprops(cellLabels, ...
                    'Area','Centroid','PixelIdxList','MajorAxisLength', 'MinorAxisLength');

            %Remove non-existing data
            cellData([cellData.Area] ==  0) = [];
            
            %Get intensity data. Names: PropertyChanName
            for iC = 1:bfReader.sizeC
                currImage = bfReader.getPlane(1, iC, iT);
                
                for iCell = 1:numel(cellData)
                    cellData(iCell).(['TotalInt',bfReader.channelNames{iC}]) = ...
                        sum(currImage(cellData(iCell).PixelIdxList));
                    
                    cellData(iCell).(['MaxInt',bfReader.channelNames{iC}]) = ...
                        max(currImage(cellData(iCell).PixelIdxList));
                    
                    cellData(iCell).(['MinInt',bfReader.channelNames{iC}]) = ...
                        min(currImage(cellData(iCell).PixelIdxList));
                    
                    cellData(iCell).(['HighInt',bfReader.channelNames{iC}]) = ...
                        prctile(currImage(cellData(iCell).PixelIdxList),hiIntLevel);
                    
                end
            end
            
        end
        
        function imgOut = makeAnnotatedImage(iT, baseImage, cellMasks, trackData)
            %MAKEANNOTATEDIMAGE  Make annotated images
            
            imgOut = MECTracker.showoverlay(MECTracker.normalizeimg(double(baseImage)),...
                bwperim(cellMasks),[0 1 0]);
            
            %Write frame number on top right
            imgOut = insertText(imgOut,[size(baseImage,2), 1],iT,...
                'BoxOpacity',0,'TextColor','white','AnchorPoint','RightTop');
                        
%             if iT >= 2
                
                for iTrack = 1:trackData.NumTracks
                    
                    currTrack = trackData.getTrack(iTrack);
                    
                    if iT > currTrack.FirstFrame && iT <= currTrack.LastFrame
                        
                        trackCentroid = cat(2,currTrack.Data.Centroid);
                        
                        imgOut = insertShape(imgOut, 'line', trackCentroid, 'color','white');
                        
                        imgOut = insertText(imgOut, currTrack.Data(end).Centroid, iTrack,...
                            'BoxOpacity', 0,'TextColor','yellow');
                        
                    end                  
                    
                end
                
%             end
        end
        
    end
    
    
end