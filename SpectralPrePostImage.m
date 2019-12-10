classdef SpectralPrePostImage
    %SPECTRALPREPOSTIMAGE  Class for pre- and post- illumination images
    %
    %  imgObj = SPECTRALPREPOSTIMAGE(filename) returns an object
    %  that will contain images and image data based on the filename
    %  specified. The filename can be either of the post- or the
    %  pre-illumination image (the code will look for the other).
    %
    %  imgObj = SPECTRALPREPOSTIMAGE will provide a dialog box to select
    %  the file of interest.
    %
    %  See also: SpectralPrePostImage.processData
    
    % TODO:
    %  * Plotting function?
    %
    % This class represents a pair of pre/post images.
    %
    % The images should be labelled "post <filename>" and "pre <filename>".
    % The images are expected to be stills (sizeT == 1). Cells cannot
    % change much between the two

    properties (Dependent)
        datasetName        
    end
    
    properties (SetAccess = private, Hidden)
        preIllumImageObj    %Bioformats Image Object for pre-illumination
        postIllumImageObj   %Bioformats Image Object for post-illumination
    end
    
    properties (SetAccess = private)
        NumCells            %Number of cells in image
        CellData            %Structure to hold cell data
        WavelengthData      %Holds detector wavelength calibration data
        Masks               %Holds the mask data
    end
    
    methods
        
        function obj = SpectralPrePostImage(varargin)
            %SPECTRALPREPOSTIMAGE  Constructor function
            %
            %  Upon construction, the object will call the processData
            %  function to get cell data.
            %
            %  See also: SpectralPrePostImage.processData

            if nargin == 0
                [fname, fpath] = uigetfile('*.nd2','Select a file');
                
                if fname == 0
                    error('SpectralPrePostImage:NoFileSelected','Select or specify a file.')
                end
            else
                [fpath,fname,fext] = fileparts(varargin{1});
                fname = [fname, fext];
            end
            
            %Look for a wavelength calibration file
            csvFileList = dir(fullfile(fpath,'*.csv'));
            if numel(csvFileList) == 1
                obj = obj.loadWavelengthData(fullfile(fpath,csvFileList.name));
            else
                obj = obj.loadWavelengthData('default');
            end
                       
            %Create the Bioformats Image objects            
            obj = obj.loadImages(fullfile(fpath,fname));
            
            %Get masks
            obj = obj.getMasks;
            
            %Get cell data
            obj = obj.ProcessImages;
            
        end
        
        function numCells = numel(obj)
            
            numCells = obj.NumCells;
            
        end
 
        function structOut = getCellData(obj, iCell)
            %GETCELLDATA  Get cell data
            %
            %  S = SpectralPrePostImage.GETCELLDATA(N) gets the data from
            %  cell N, returning it as a structure S.
            
            structOut = obj.CellData(iCell);
            
        end
        
        function imgOut = getAnnotatedImage(obj, channel, type)
            %GETANNOTATEDIMAGE  Makes an annotated image
            %
            %  I = SpectralPrePostImage.GETANNOTATEDIMAGE(channel, type)
            %  will generate an image with each cell labelled, returned as
            %  a matrix I. Channel should be the channel number. 
            %
            %  The following image types can be specified:
            %       'pre'  - Pre-illumination image
            %      'post' - Post-illumination image
            %     'merge' - Both images side by side (pre-illumination
            %               image on the left)
            %
            %  See also: SpectralPrePostImage.getImage            
            
            %Set defaults
            if ~exist('channel','var')
                channel = 1;
            end
            
            if ~exist('type','var')
                type = 'merge';
            end
            
            %Options for the text color and box opacity (see insertText)
            labelOptions = {'TextColor','blue','BoxOpacity',0.5,'AnchorPoint','CenterTop'};
            
            switch lower(type)
                
                case {'pre', 'post'}
                    
                    frameImg = obj.getImage(channel, type);
                    imgOut = showoverlay(normalizeimg(frameImg),bwperim(obj.Masks.Center),[0 1 0]);
                    imgOut = showoverlay(imgOut,bwperim(obj.Masks.Edge),[1 0 0]);
                    
                    %Label the cells
                    for iCell = 1:obj.NumCells
                        imgOut = insertText(imgOut, obj.CellData(iCell).Centroid, iCell,...
                            labelOptions{:});                        
                    end
                    
                case 'merge'
                    
                    preImg = obj.getImage(channel, 'pre');
                    imgPre = showoverlay(normalizeimg(preImg),bwperim(obj.Masks.Center),[0 1 0]);
                    imgPre = showoverlay(imgPre,bwperim(obj.Masks.Edge),[1 0 0]);
                    
                    postImg = obj.getImage(channel, 'post');
                    imgPost = showoverlay(normalizeimg(postImg),bwperim(obj.Masks.Center),[0 1 0]);
                    imgPost = showoverlay(imgPost,bwperim(obj.Masks.Edge),[1 0 0]);
                    
                    %Label the cells
                    for iCell = 1:obj.NumCells
                        imgPre = insertText(imgPre, obj.CellData(iCell).Centroid, iCell,...
                           labelOptions{:});
                        
                        imgPost = insertText(imgPost, obj.CellData(iCell).Centroid, iCell,...
                           labelOptions{:});
                    end
                    
                    %Label the images
                    imgPre = insertText(imgPre, [1, 1], 'Pre','FontSize',20);
                    imgPost = insertText(imgPost, [1, 1], 'Post','FontSize',20);
                    
                    imgOut = [imgPre,...
                        ones(obj.preIllumImageObj.height,2,3) * 65535,...
                        imgPost];
                    
            end
            
        end
        
        function frameImage = getImage(obj, channel, type)
            %GETIMAGE  Get image from ND2 file
            %
            %  I = SpectralPrePostImage.GETIMAGE(channel, type) returns the
            %  image frame from the specified channel to matrix I.
            %
            %  The following image types can be specified:
            %       'pre'  - Pre-illumination image
            %      'post' - Post-illumination image
            %     'merge' - Both images side by side (pre-illumination
            %               image on the left)
            %
            %  Note that these images are not contrast adjusted. If the
            %  image does not appear when using imshow, try setting the
            %  limits to automatic, e.g. imshow(I,[]).
            
            switch lower(type)
                
                case 'pre'
                    frameImage = obj.preIllumImageObj.getPlane(1,channel,1);
                    
                case 'post'
                    frameImage = obj.postIllumImageObj.getPlane(1,channel,1);
                    
                case 'merge'
                    
                    frameImage = [obj.preIllumImageObj.getPlane(1,channel,1),...
                        zeros(obj.preIllumImageObj.height,2),...     %Draw a thin line between the images
                        obj.postIllumImageObj.getPlane(1,channel,1)];
                
                otherwise
                    error('SpectralPrePostImage:InvalidImageType',...
                        'Invalid type %s. Expected ''pre'', ''post'', or ''merge''', type);
            end
            
        end

        function exportData(obj, outputFN)
            %EXPORTDATA  Exports data to CSV
            %
            %  SpectralPrePostImage.exportData(filename) will export the
            %  data to the file specified. The exported format will be
            %  CSV.
            
            if obj.NumCells == 0
                error('SpectralPrePostImage:NoData','No data currently exists.');
            end
            
            [fpath, fname] = fileparts(outputFN);
            if isempty(fname)
                fname = 'exportedData';
            end
            
            fname = sprintf('%s_%s_%s.csv', fname, obj.datasetName, datestr(today));
            
            fid = fopen(fullfile(fpath,fname),'w');
            
            %Write headers
            fprintf(fid,'CellID, Wavelength, Area Center, Area Edge, Total Center (pre), Mean Center (pre), Total Edge (pre),  Mean Edge (pre), Total Center (post), Mean Center (post), Total Edge (Center), Mean Edge (post)\n');
            
            for ii = 1:obj.NumCells
                
                for iChan = 1:numel(obj.CellData(ii).Channel)
                    
                    if iChan == 1
                        fprintf(fid,'%d, %d, %d, %d, %d, %.2f, %d, %.2f,  %d, %.2f,  %d, %.2f\n',...
                            ii, ...                                                                         %ID
                            0.5 * (obj.WavelengthData(iChan).End + obj.WavelengthData(iChan).Start),...     %WAVELENGTH
                            obj.CellData(ii).AreaCenter,...                                                 %Area Center (PRE)
                            obj.CellData(ii).AreaEdge,...                                                   %Area Edge (PRE)
                            obj.CellData(ii).Channel(iChan).PreIllumCenterTotalInt, ...                     %TOTAL CENTER (PRE)
                            obj.CellData(ii).Channel(iChan).PreIllumCenterMeanInt, ...                      %Mean Center (Pre)
                            obj.CellData(ii).Channel(iChan).PreIllumEdgeTotalInt,...                         %Total Center (pre)
                            obj.CellData(ii).Channel(iChan).PreIllumEdgeMeanInt,...                         %Mean Edge (pre)
                            obj.CellData(ii).Channel(iChan).PostIllumCenterTotalInt,...                        %Total Center (post)
                            obj.CellData(ii).Channel(iChan).PostIllumCenterMeanInt,...                      %Mean Center (post)
                            obj.CellData(ii).Channel(iChan).PostIllumEdgeTotalInt,...                         %Total Edge (post)
                            obj.CellData(ii).Channel(iChan).PostIllumEdgeMeanInt);                          %Mean Edge (post)
                            
                    else
                        fprintf(fid,'  ,%d, %d, %d, %d, %.2f, %d, %.2f,  %d, %.2f,  %d, %.2f\n',...
                            0.5 * (obj.WavelengthData(iChan).End + obj.WavelengthData(iChan).Start),...     %WAVELENGTH
                            obj.CellData(ii).AreaCenter,...                                  %Area Center (PRE)
                            obj.CellData(ii).AreaEdge,...                                    %Area Edge (PRE)
                            obj.CellData(ii).Channel(iChan).PreIllumCenterTotalInt, ...                     %TOTAL CENTER (PRE)
                            obj.CellData(ii).Channel(iChan).PreIllumCenterMeanInt, ...                      %Mean Center (Pre)
                            obj.CellData(ii).Channel(iChan).PreIllumEdgeTotalInt,...                      %Total Center (pre)
                            obj.CellData(ii).Channel(iChan).PreIllumEdgeMeanInt,...                         %Mean Edge (pre)
                            obj.CellData(ii).Channel(iChan).PostIllumCenterTotalInt,...                        %Total Center (post)
                            obj.CellData(ii).Channel(iChan).PostIllumCenterMeanInt,...                      %Mean Center (post)
                            obj.CellData(ii).Channel(iChan).PostIllumEdgeTotalInt,...                         %Total Edge (post)
                            obj.CellData(ii).Channel(iChan).PostIllumEdgeMeanInt);                          %Mean Edge (post)
                    end
                    
                end
                
            end
            
            fclose(fid);
            
        end
        
        function fname = get.datasetName(obj)
            
            %Get the current filename
            [~,fname] = fileparts(obj.postIllumImageObj.filename);
            
            %Look for either pre or post
            fname = fname(6:end);
            
        end
        
    end
       
    methods (Hidden) %Functions used to process images in call order
       
       function obj = loadWavelengthData(obj,varargin)
           %LOADWAVELENGTHDATA  Load wavelength data from calibration file
           %
           %  SpectralPrePostImage.LOADWAVELENGTHDATA(filename) will load
           %  the file specified and attempt to read in wavelength data.
           %  The expected data format should be as a comma-delimited,
           %  with the following columns: Channel, Start (nm), End (nm),
           %  Diff.
           %
           %  SpectralPrePostImage.LOADWAVELENGTHDATA will open a dialog
           %  box and ask the user to select the file containing
           %  wavelength data.
           %
           %  SpectralPrePostImage.LOADWAVELENGTHDATA('default') will
           %  attempt to read the channel names from the ND2 files and use
           %  those to determine the wavelength range.
           
           if isempty(varargin)
               [fname,fpath] = uigetfile({'*.txt; *.csv','Wavelength data (*.txt, *.csv)'},'Select wavelength data');
               filename = fullfile(fpath, fname);
               
               if filename == 0
                   error('SpectralPrePostImage:NoWavelengthFileSelected','Wavelength data not selected.');
               end
               
           elseif strcmpi(varargin{1},'default')
               %Otherwise use the default values from the Bioformats Image
               %objects
               for iChannel = 1:obj.preIllumImageObj.sizeC - 1
                   obj.WavelengthData(iChannel).Start = str2double(obj.preIllumImageObj.channelNames{iChannel});
                   obj.WavelengthData(iChannel).End = str2double(obj.preIllumImageObj.channelNames{iChannel + 1}) - 0.01;
               end
               
               %Special handling for last channel
               obj.WavelengthData(obj.preIllumImageObj.sizeC).Start = str2double(obj.preIllumImageObj.channelNames{obj.preIllumImageObj.sizeC});
               obj.WavelengthData(obj.preIllumImageObj.sizeC).End = str2double(obj.preIllumImageObj.channelNames{obj.preIllumImageObj.sizeC}) - 0.01;
               
           else
               filename = varargin{1};
           end
           
           %Delete old data if it exists
           obj.WavelengthData = [];
           
           %Open the file
           fid = fopen(filename,'r');
           
           while ~feof(fid)
               
               currLine = strsplit(fgetl(fid),',');
               if isempty(regexp(currLine{1},'\d*','ONCE'))
                   %Skip
               else
                   %Data should be Channel, Start (nm), End (nm), Diff
                   obj.WavelengthData(end + 1).Start = str2double(currLine{2});
                   obj.WavelengthData(end).End = str2double(currLine{3});
               end
               
           end
           
           %Close the file
           fclose(fid);
           
       end
        
       function obj = loadImages(obj,inputFilename)
           %LOADIMAGES  Create BioformatsImage objects
           %
           %  O = O.LOADIMAGES(filename) creates the required
           %  BioformatsImage objects from the input file.
           
           %Get the current filename
           [fpath,fname,fext] = fileparts(inputFilename);
           
           %Look for either pre or post
           preStrInd = strfind(fname,'pre ');
           postStrInd = strfind(fname,'post ');
           
           if ~isempty(preStrInd)
               %This filename provided is the pre-illumination file
               preFname = fname;
               postFname = ['post',preFname(4:end)];
               
           elseif ~isempty(postStrInd)
               %The filename provided is the post-illumination file
               postFname = fname;
               preFname = ['pre', postFname(5:end)];
               
           else
               error('SpectralPrePostImage:InvalidFilename',...
                   'Expected the filename to start with ''post'' or ''pre''.');
               
           end
           
           obj.preIllumImageObj = BioformatsImage(fullfile(fpath,[preFname, fext]));
           obj.postIllumImageObj = BioformatsImage(fullfile(fpath,[postFname, fext]));
           
       end
        
       function cellMask = segmentCells(obj, thLvl, varargin)
            %SEGMENTCELLS  Segments and displays cells
            
            ip = inputParser;
            ip.addParameter('Watershed', true, @(x) islogical(x));
            ip.addParameter('ROI',[]);
            ip.parse(varargin{:});
            
            %             %Sum the channels to get a clearer image of the cells
            %             imgData = obj.getCStack(iT,'ROI',ip.Results.ROI);
            
            %Sum the channels to get a clearer image of the cells
            imgData = obj.postIllumImageObj.getPlane(1,1,1,'ROI',ip.Results.ROI);
            
            thImg = imgData > thLvl;
            thImg = imopen(thImg,strel('disk',3));
            thImg = imfill(thImg,'holes');
            thImg = bwareaopen(thImg,100);
            %             thImg = imclearborder(thImg);
            
            if ip.Results.Watershed
                
                dd = -bwdist(~thImg);
                dd = imhmin(dd,0.8);
                dd(~thImg) = -Inf;
                
                cellLabels = watershed(dd);
                cellLabels = imclearborder(cellLabels);
%                 imshow(label2rgb(LL));
%                 keyboard
            else
                
                cellLabels = bwlabel(thImg);

            end           
            
            if nargout > 0
                cellMask = cellLabels;
            else
                showoverlay(normalizeimg(imgData),bwperim(cellLabels),[0 1 0]);
            end
            
       end
       
       function obj = getMasks(obj)
            %GETMASKS  Segment cells and make the center and edge masks
            
            obj.Masks.Full = obj.segmentCells(100);
            
            %Erode the mask to get the center and edge masks
            obj.Masks.Center = bwmorph(obj.Masks.Full,'shrink',3);
            obj.Masks.Edge = obj.Masks.Full & ~obj.Masks.Center;
            
        end
        
       function obj = ProcessImages(obj)
           %PROCESSIMAGES  Get cell data
           %
           %  PROCESSIMAGES(Obj) populates the object data structure with
           %  cell data. This method calls the function segmentCells (with
           %  a threshold greyscale of 100 by default), then generates a
           %  center mask and an edge mask by shrinking the cell labels by
           %  3 pixels.
           %
           %  The intensity of the images pre- and post-illumination is
           %  calculated for both the center and the edge masks, and for
           %  each channel.
           
           %Check that BioformatsImage objects are available
           if isempty(obj.preIllumImageObj) || isempty(obj.preIllumImageObj)
               error('SpectralPrePostImage:LoadImagesFirst','Images are not loaded.')
           end
           
           for iChannel = 1:obj.preIllumImageObj.sizeC
               %Get data
               preIllumImage = obj.preIllumImageObj.getPlane(1,iChannel,1);
               preIllumData_center = regionprops(obj.Masks.Center,preIllumImage,'MeanIntensity','Centroid','Area','PixelValues');
               preIllumData_edge = regionprops(obj.Masks.Edge,preIllumImage,'MeanIntensity','Area','PixelValues');
               
               postIllumImage = obj.postIllumImageObj.getPlane(1,iChannel,1);
               postIllumData_center = regionprops(obj.Masks.Center,postIllumImage,'MeanIntensity','PixelValues');
               postIllumData_edge = regionprops(obj.Masks.Edge,postIllumImage,'MeanIntensity','PixelValues');
               
               %Add cells to the data property
               for iCell = 1:numel(preIllumData_center)
                   obj.CellData(iCell).Centroid = preIllumData_center(iCell).Centroid;
                   obj.CellData(iCell).AreaCenter = preIllumData_center(iCell).Area;
                   obj.CellData(iCell).AreaEdge = preIllumData_edge(iCell).Area;
                   
                   %--- Pre-illumination data ---%
                   obj.CellData(iCell).Channel(iChannel).PreIllumCenterMeanInt = preIllumData_center(iCell).MeanIntensity;
                   obj.CellData(iCell).Channel(iChannel).PreIllumEdgeMeanInt = preIllumData_edge(iCell).MeanIntensity;
                   
                   obj.CellData(iCell).Channel(iChannel).PreIllumCenterTotalInt = sum(preIllumData_center(iCell).PixelValues);
                   obj.CellData(iCell).Channel(iChannel).PreIllumEdgeTotalInt = sum(preIllumData_edge(iCell).PixelValues);
                   
                   %--- Post-illumination data ---%
                   obj.CellData(iCell).Channel(iChannel).PostIllumCenterMeanInt = postIllumData_center(iCell).MeanIntensity;
                   obj.CellData(iCell).Channel(iChannel).PostIllumEdgeMeanInt = postIllumData_edge(iCell).MeanIntensity;
                   
                   obj.CellData(iCell).Channel(iChannel).PostIllumCenterTotalInt = sum(postIllumData_center(iCell).PixelValues);
                   obj.CellData(iCell).Channel(iChannel).PostIllumEdgeTotalInt = sum(postIllumData_edge(iCell).PixelValues);
                   
               end
               
           end
           
           obj.NumCells = numel(obj.CellData);
           
       end
       
    end
    
end