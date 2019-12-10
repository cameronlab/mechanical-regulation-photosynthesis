classdef BioformatsImage
    % BIOFORMATSIMAGE  Class to read microscope images
    %
    %   R = BIOFORMATSIMAGE(filename) creates a new BioformatsImage object
    %   R which is linked to the filename specified.
    %
    %   To get a specific image plane, use ZCTS coordinates:
    %
    %     I = R.getPlane(iZ, iC, iT, iS);
    %
    %   where iZ = z-plane index, iC = channel index, iT = timepoint index,
    %   and iS = series. iS is optional; is not specified, the current
    %   series will be used.
    %
    %   **Bioformats toolbox users**: Be aware that ZCTS is one-based in
    %   this class.
    %
    %   You can also get an ROI:
    %
    %     I = R.getPlane(iZ, iC, iT, 'ROI',[top left, width, height]);
    %
    %   For large images, it can be helpful to get an image in tiles:
    %
    %     I = R.getTile(iZ, iC, iT, numTiles, tileIndex);
    %
    %   Note: getTile does not support an iS input.
    %
    %   where numTiles is the number of tiles to divide into [rows, cols],
    %   and the tileIndex is the desired tile number. The tileIndex follows
    %   MATLAB numbering convention and increases down the rows first, then
    %   along the columns.
    %
    %   NOTE: This class uses the Bioformats toolbox, a standalone Java
    %   library for reading and writing life sciences image formats
    %   developed by the Open Microscopy Environment team
    %   (http://http://www.openmicroscopy.org). Whenever an object is
    %   created, it will check if the toolbox exists. If it does not, it
    %   will download and install a compatible version (currently 5.4.0) to
    %   the MATLAB path.
    
    properties (AbortSet)   %Filename
        filename = '';
    end
    
    properties (Dependent)  %Image file attributes
        width
        height
        
        series
        seriesCount
        
        sizeZ   %Number of z-planes
        sizeC   %Number of channels
        sizeT   %Number of timepoints
        
        channelNames
        
        pxSize
        pxUnits
    end
    
    properties (Transient, Hidden, SetAccess = private)
        cleanup
        bfReader
        metadata
    end
    
    properties (Access = private)   %Bioformats toolbox download URL
        bfTbxURL = 'http://downloads.openmicroscopy.org/bio-formats/5.7.1/artifacts/bfmatlab.zip';
    end
    
    %--- Main methods ---%
    
    methods %Constructor, getters and setters
        
        function obj = BioformatsImage(varargin)
            % Create a new class object
            %
            % R = BIOFORMATSIMAGE(filename) returns the object in R.
            
            %Check if the toolbox exists
            if ~obj.bfTbxExist
                obj.installBFtbx;
            end
            
            if nargin > 0
                if nargin == 1
                    obj.filename = varargin{1};
                else
                    error('BioformatsImage:TooManyInputArguments',...
                        'Too many input arguments. Was expecting one.')
                end
            end
            
            obj.cleanup = onCleanup(@()delete(obj));
            
        end
        
        function obj = set.filename(obj,filename)
            
            %Verify that filename is a string
            if ~ischar(filename)
                error('BioformatsImage:FilenameNotString',...
                    'Input filename must be a string.');
            elseif isempty(filename)
                error('BioformatsImage:NoFilename',...
                    'Please specify an image file.');
            end
            
            %Check that the file exists
            if ~exist(filename,'file')
                error('BioformatsImage:CannotFindFile',...
                    'Could not find %s. Provide the full filename or make sure the image folder is on the MATLAB path.',filename);
            end
            
            %Make sure that the full path to the file is stored
            [fPath,fName,fExt] = fileparts(filename);
            
            if ~isempty(fPath)
                %Reconstruct the full path (this will correct any
                %system-dependent file separator strings)
                obj.filename = fullfile(fPath,[fName,fExt]);
                
            else
                %Since the file must exist on the MATLAB path (we checked
                %this above), we can use which to determine the full path
                %to the file
                obj.filename = which(filename);
                
            end
            
            %Get a reader for the image
            obj = obj.getReader;
            
        end
        
        function obj = set.series(obj,iS)
            %Set series number
            
            if isempty(obj.bfReader)
                obj = obj.getReader;
            end
            
            if rem(iS,1) ~= 0
                error('BioformatsImage:SeriesNumberMustBeInteger',...
                    'Expected series number to be an integer.')
            end
            
            if iS <= 0 || iS > obj.seriesCount
                error('BioformatsImage:RequestedSeriesOutOfBounds',...
                    'Requested series must be > 0 and <= than %d.',obj.seriesCount);
            end
            
            obj.bfReader.setSeries(iS - 1);
            
        end
        
        function seriesOut = get.series(obj)
            
            if ~obj.bfReaderExist
                seriesOut = NaN;
            else
                seriesOut = obj.bfReader.getSeries + 1;
            end
            
        end
        
        function height = get.height(obj)
            
            if ~isempty(obj.bfReader)
                height = obj.bfReader.getSizeY;
            end
            
        end
        
        function width = get.width(obj)
            
            if ~isempty(obj.bfReader)
                width = obj.bfReader.getSizeX;
            end
            
        end
        
        function sizeZ = get.sizeZ(obj)
            
            if ~isempty(obj.bfReader)
                sizeZ = obj.bfReader.getSizeZ;
            end
            
        end
        
        function sizeC = get.sizeC(obj)
            
            if ~isempty(obj.bfReader)
                sizeC = obj.bfReader.getSizeC;
            end
            
        end
        
        function sizeT = get.sizeT(obj)
            
            if ~isempty(obj.bfReader)
                sizeT = obj.bfReader.getSizeT;
            end
            
        end
        
        function seriesCount = get.seriesCount(obj)
            
            seriesCount = obj.bfReader.getSeriesCount;
            
        end
        
        function metadata = get.metadata(obj)
            
            metadata = obj.bfReader.getMetadataStore;
            
        end
        
        function chNames = get.channelNames(obj)
            chNames = cell(1,obj.sizeC);
            
            for iC = 1:obj.sizeC
                chNames{iC} = obj.metadata.getChannelName(0, iC - 1).char;
            end
        end
        
        function pxSize = get.pxSize(obj)
            
            pxSizeX = double(obj.metadata.getPixelsPhysicalSizeX(0).value);
            pxSizeY = double(obj.metadata.getPixelsPhysicalSizeY(0).value);
            
            pxSize = [pxSizeX, pxSizeY];
            
        end
        
        function pxUnit = get.pxUnits(obj)
            
            pxInfoStr = char(obj.metadata.getPixelsPhysicalSizeX(0).toString);
            
            pxUnit = BioformatsImage.unitStrToChar(pxInfoStr);
            
        end
        
        function versionOut = version(obj)
            
            verStr = ver(class(obj));
            versionOut = verStr.Version;
            
        end
        
        function delete(obj)
            
            obj.bfReader.close;
            
        end

    end
    
    methods %Base functions
        
        function obj = getReader(obj)
            %Get a Bioformats Reader object
            %
            % You should not need to call the function again, but I'm
            % leaving it here in case it is necessary to reboot it.
            
            if ~isempty(obj.filename)
                obj.bfReader = bfGetReader(obj.filename);
            else
                error('BioformatsImage:EmptyFilename',...
                    'Set filename first.');
            end
            
            %Set current series to 1
            obj.series = 1;
            
        end
        
        function [imgOut, timestamp, tsUnits] = getPlane(obj, iZ, iC, iT, varargin)
            %Get image at specified index
            %
            %  imgOut = getImage([iZ, iC, iT])
            %
            %  imgOut = getImage(iZ, iC, iT, iS);
            %
            %  imgOut = getImage(iZ, iC, iT, iS, ROI);
            %
            %  ROI = [XMIN YMIN WIDTH HEIGHT];
            
            %Parse the variable input argument
            ip = inputParser;
            ip.addOptional('iSeries',NaN,@(x) isnumeric(x) && isscalar(x));
            ip.addParameter('ROI',[],@(x) isempty(x) || (numel(x) == 4 && all(x > 0)));
            ip.parse(varargin{:});
            
            %Check that reader object already exists
            if ~bfReaderExist(obj)
                %If not, get a reader for the file
                obj = obj.getReader;
            end
            
            %Change the series (if set)
            if ~isnan(ip.Results.iSeries)
                obj.series = ip.Results.iSeries;
            end
            
            %Resolve the channel name
            iC = obj.channelname2ind(iC);
            
            %Get the image
            if ~isempty(ip.Results.ROI)
                xMin = ip.Results.ROI(1);
                yMin = ip.Results.ROI(2);
                roiWidth = ip.Results.ROI(3);
                roiHeight = ip.Results.ROI(4);
                imgOut = bfGetPlane(obj.bfReader,obj.getIndex(iZ, iC, iT),xMin,yMin,roiWidth,roiHeight);
            else
                %Get full image
                imgOut = bfGetPlane(obj.bfReader,obj.getIndex(iZ, iC, iT));
            end
            
            if nargout > 1
                %Get the timestamp
                [timestamp, tsUnits] = obj.getTimestamps(iZ,iC,iT);
            end
            
        end
        
        function [imgOut, timestamp, tsUnits] = getXYplane(obj, iC, iT, XYloc, varargin)
            %GETXYPLANE  Get plane from multi-XY images
            %
            %  I = O.GETXYPLANE(iZ, iC, iT, XYLocation)
            %
            %  Currently, the code assumes that iZ = 1.
            %
            %  When saving multi-XY images in Nikon Jobs, there seems to be
            %  an issue where the images are not stored in the correct
            %  order. The image sequence then becomes interleaved:
            %      {(XY1, T1, C1), (XY1, T1, C2), (XY2, T1, C1), (XY2, T1,
            %      C2),...}
            %
            %  The series count corresponds to the number of XY locations.
            %
            %  To solve this issue, this function recalculates the image
            %  index assuming this interleaved order.
            
            ip = inputParser;
            ip.addParameter('ROI',[],@(x) numel(x) == 4 && all(x > 0));
            ip.parse(varargin{:});
            
            %Check that reader object already exists
            if ~bfReaderExist(obj)
                %Get a reader for the file
                obj = obj.getReader;
            end
            
            %Resolve the channel index
            iC = obj.channelname2ind(iC);
            
            %NOTE: I have not verified that calculations involving iZ is
            %correct due to lack of data. Therefore, force this to = 1 to
            %prevent the user from getting incorrect data.
            iZ = 1;
            
            %Calculate the image index
            imgIndex = XYloc + (iT - 1) * obj.seriesCount - 1;
            
            %Calculate the actual series number and timepoint
            iS = floor(imgIndex/obj.sizeT);
            
            iT = imgIndex - (iS * obj.sizeT) + 1;
            
            %Get the image plane
            obj.series = iS + 1;
            
            %Get the image
            imgOut = obj.getPlane(iZ, iC, iT, 'ROI',ip.Results.ROI);
            
            %Get timestamp if output is assigned
            if nargout > 1
                %Get the timestamp
                [timestamp, tsUnits] = obj.getTimestamps(iZ, iC, iT);
            end
        end
        
        function [timestamps, tsunits] = getTimestamps(obj, iZ, iC, varargin)
            %GETTIMESTAMPS  Get timestamps from the specified channel
            %
            %  [tsData, tsUnits] = O.GETTIMESTAMPS(channelIndex) returns
            %  the timestamp information in a vector tsData. tsUnits is a
            %  character array which contains the unit of the timestamp
            %  data.
            
            ip = inputParser;
            ip.addOptional('TimeRange',Inf,@(x) all(isinf(x)) || isnumeric(x));
            ip.parse(varargin{:});
            
            %Resolve channel name
            iC = obj.channelname2ind(iC);
            
            if isinf(ip.Results.TimeRange)
                timeRange = 1:obj.sizeT;
            else
                timeRange = ip.Results.TimeRange;
            end
            
            %Initialize a vector for the timestamps
            timestamps = zeros(1, numel(timeRange));
            
            for iT = timeRange
                %Resolve the bioformats index
                bfIndex = obj.getIndex(iZ,iC,iT);
                currTS = obj.metadata.getPlaneDeltaT(obj.series - 1,bfIndex - 1);
                if ~isempty(currTS)
                    timestamps(iT) = double(currTS.value);
                else
                    warning('Timestamp empty');
                    tsunits = '';
                    return;
                end
            end
            
            %Get the unit string
            tsStr = char(obj.metadata.getPlaneDeltaT(obj.series - 1,bfIndex - 1).toString);
            
            %Convert the string into a char
            tsunits = BioformatsImage.unitStrToChar(tsStr);
            
        end
        
        function [tileDataOut, roiOut] = getTile(obj, planeSpec, numTiles, tileIndex)
            
            %Check that the planeSpec has 3 entries for iZ, iC, iT
            if ~(numel(planeSpec) == 3)
                error('The plane specification should be in [iZ, iC, iT] format');
            end
            
            ip = inputParser;
            ip.addRequired('NumTiles',@(x) numel(x) == 2 && all(x > 0));
            ip.addRequired('TileRange', @(x) all(x > 0));
            ip.parse(numTiles,tileIndex);
            
            %Check that the tile range fits within the total number of
            %tiles
            if any(ip.Results.TileRange > prod(ip.Results.NumTiles))
                error('BioformatsImage:IndexExceedsTileDimensions',...
                    'Index exceeds tile dimensions. Maximum index is %d.',...
                    prod(ip.Results.NumTiles))
            end
            
            %Create a storage cell (for multiple tile data, which might have different sizes)
            tileDataOut = cell(numel(ip.Results.TileRange),1);
            
            %Create a storage for the roi
            roiOut = zeros(numel(ip.Results.TileRange),4);
            
            for iT = 1:numel(ip.Results.TileRange)
                tileIndex = ip.Results.TileRange(iT);
                roi = getTileIndices(obj, ip.Results.NumTiles, tileIndex);
                
                %Get the ROI
                tileDataOut{iT} = obj.getPlane(planeSpec(1),planeSpec(2),planeSpec(3), 'ROI', roi);
                roiOut(iT,:) = roi;
            end
            
            %If it's only a single tile, output a matrix rather than a cell
            if numel(tileDataOut) == 1
                tileDataOut = tileDataOut{:};
            end
            
        end
        
    end
    
    %--- Aux methods ---%    
    
    methods (Access = private, Static)
        
        function tbxExists = bfTbxExist
            %Check if the BioFormats toolbox is installed
            
            if ~exist('bfGetReader','file')
                
                %Check if there is a bfmatlab folder in the toolbox folder
                
                tbxFolder = fileparts(which('BioformatsImage'));
                
                if exist(fullfile(tbxFolder,'bfmatlab'),'dir')
                    %Folder exists so add it to the path
                    addpath(fullfile(tbxFolder,'bfmatlab'));
                else
                    tbxExists = false;
                    return;
                end
            end
            
            tbxExists = true;
            
        end
        
        function unitStr = unitStrToChar(strIn)
            %Parses an input string for the unit
            
            %Find the 'unit[' pattern
            unitStartIdx = strfind(strIn,'unit[');
            %Offset by 5 characters to get start of unit string
            unitStartIdx = unitStartIdx + 5;
            
            %Find the closing brace, starting from the unit string
            unitEndIdx = strfind(strIn(unitStartIdx:end),']');
            %Reduce by 1 character to not include the closing brace
            unitEndIdx = unitEndIdx - 2;
            
            %Output the pixel unit
            unitStr = strIn(unitStartIdx: unitStartIdx + unitEndIdx);
            
        end
        
    end
    
    methods (Access = private)
        
        function installBFtbx(obj)
            %Downloads and installs the Bioformats Toolbox
            %
            % The default install location is the toolbox folder
            
            tbxFolder = fileparts(which('BioformatsImage'));
            
            %Create a temporary directory
            if ~exist(fullfile(tbxFolder,'temp'),'dir')
                mkdir(fullfile(tbxFolder,'temp'))
            end
            
            %Download the toolbox from the OME website
            fprintf('Downloading Bioformats toolbox...');
            
            fn = websave(fullfile(tbxFolder,'temp','bfmatlab.zip'),obj.bfTbxURL);
            fprintf('Done.\n');
            
            fprintf('Installing toolbox...')
            %Unzip the toolbox
            unzip(fn,tbxFolder);
            
            %Add toolbox to path
            addpath(fullfile(tbxFolder,'bfmatlab'))
            fprintf('Done.\n');
            
            %Remove temporary files and directory
            if exist(fullfile(tbxFolder,'temp'),'dir')
                rmdir(fullfile(tbxFolder,'temp'),'s')
            end
            
        end
        
        function bfrExists = bfReaderExist(obj)
            %Checks to see if the bfReader property contains a bfreader
            %object
            
            if ~isempty(obj.bfReader) && isa(obj.bfReader,'loci.formats.ChannelSeparator')
                bfrExists = true;
            else
                bfrExists = false;
            end
            
        end
        
        function chInd = channelname2ind(obj,channelIn)
            %Convert channel name to index
            
            %If it is already a number, then nothing to do
            if isnumeric(channelIn)
                chInd = channelIn;
                return
            end
            
            %If it is a cell (happens if iterating over channelNames), then
            %make sure that the cell has a single string entry
            if iscell(channelIn)
                if numel(channelIn) > 1
                    error('BioformatsImage:channelname2ind:TooManyChannels',...
                        'Too many channel names specified at once. Only one channel expected.');                    
                end
                
                channelIn = channelIn{1};
                
            end
            
            %Find the matching name
            chTF = strcmpi(channelIn, obj.channelNames);
            
            if ~any(chTF)
                error('BioformatsImage:ChannelNameNotFound',...
                    'Could not find channel %s.',channelIn);
            end
            
            chInd = find(chTF);
            
        end
        
        function roiOut = getTileIndices(obj, numTiles, tileIndex)
            %Given a tile size and index, calculate the resulting ROI
            %vector
            %ROI = [XMIN YMIN WIDTH HEIGHT];
            
            if numel(tileIndex) == 1
                %If it's a single index, resolve it into row,col
                %coordinates
                [indRow, indCol] = ind2sub(numTiles,tileIndex);
            else
                indRow = tileIndex(1);
                indCol = tileIndex(2);
            end
            
            %Calculate the tile height/width start and end indices
            %
            % Tile indices are:
            %     start floor(height/numrows) * (tile row number - 1) + 1
            %     end floor(height/numrows) * (tile row number)
            
            %Row:
            tileHeight = floor(obj.height/numTiles(1));
            
            rowStart = tileHeight * (indRow - 1) + 1;
            
            if indRow ~= numTiles(1)
                rowEnd = tileHeight * indRow;
            else
                %If it's the last row, then just use the image height
                rowEnd = obj.height;
            end
            
            %Convert to ROI height
            roiHeight = rowEnd - rowStart + 1;
            
            %Col:
            tileWidth = floor(obj.width/numTiles(2));
            
            colStart = tileWidth * (indCol - 1) + 1;
            
            if indCol ~= numTiles(2)
                colEnd = tileWidth * indCol;
            else
                %If it's the last row, then just use the image height
                colEnd = obj.width;
            end
            
            %Convert to ROI width
            roiWidth = colEnd - colStart + 1;
            
            roiOut = [colStart, rowStart, roiWidth, roiHeight];
            
        end
        
        function imgIndex = getIndex(obj, iZ, iC, iT)
            %Get image index from ZCT coordinates.
            imgIndex = obj.bfReader.getIndex(iZ - 1, iC - 1, iT - 1) + 1;
        end
        
    end
        
end