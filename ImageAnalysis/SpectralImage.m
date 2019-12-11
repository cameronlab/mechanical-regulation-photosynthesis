classdef SpectralImage < BioformatsImage
    %SPECTRALIMAGE  Spectral Image class
    %  
    % This class represents a spectral image. It is a subclass of
    % BioformatsImage.

    properties
        
        maskList = {};
        
    end
    
    methods
        
        function obj = SpectralImage(varargin)
            %SPECTRALIMAGE  Constructor function
            %
            %  imgObj = SPECTRALIMAGE(filename) will return a SpectralImage
            %  object with the filename specified.

            %Call the superclass constructor
            obj = obj@BioformatsImage(varargin{:});
            
        end
        
        function imshow(obj, zctCoords, varargin)
            %IMSHOW  Show specified image plane
            %
            %  Syntax: imshow(obj, zctCoords, varargin)
            %
            %  varargin can either by input to imshow or ROI.
            
            ip = inputParser;
            ip.addOptional('imshowRange',NaN);
            ip.addParameter('ROI',[]);
            ip.parse(varargin{:});
            
            if isempty(obj.filename)
                error('Set the filename first')
            end  
            
            %Get image
            imgData = obj.getPlane(zctCoords(1),zctCoords(2),zctCoords(3),'ROI',ip.Results.ROI);

            if isnan(ip.Results.imshowRange)
                imshow(imgData);
            elseif isempty(ip.Results.imshowRange)
                imshow(imgData,[]);                
            else
                imshow(imgData,ip.Results.imshowRange);
            end
            
        end
        
        function plotSpectrum(obj, iT, varargin)
            %PLOTSPECTRUM  Plots the spectrum of the image
            
            spectrum = obj.getSpectrum(iT,varargin{:});
            
            plot(str2double(obj.channelNames),spectrum)
%             xticks(str2double(obj.channelNames))
            xlabel('Wavelength [nm]')
            ylabel('Intensity')
            
        end
        
        function [spectrumOut, ROI] = getSpectrum(obj, iT, varargin)
            %GETSPECTRUM  Get spectrum of the image
            %
            %  Works in three flavors: (1) full image, (2) from a
            %  selectable ROI, and (3) from a given ROI.
            %
            %  [S, ROI] = GETSPECTRUM(obj, iT) shows the first image, then
            %  allows the user to select a rectangular ROI. To plot the
            %  ROI:
            %
            %      imshow(obj.getPlane(1,1,iT))
            %      hold on
            %      rectangle('Position',ROI,'EdgeColor','green');
            %      hold off
           
            %TODO:
            % * Tidy up code
            % * Allow ROI from segmentCells
            
            if isempty(varargin)
                %Display the image and use rect to select an ROI
                
                hf = figure;
                imshow(obj, [1, 1, iT],[]);
                rectData = imrect(gca);
                ROI = round(rectData.getPosition);
                
                close(hf)
                
            elseif ischar(varargin{1})
                
                if strcmpi(varargin{1},'full')
                    
                    ROI = [1, 1, obj.width, obj.height];
                    
                else
                    error('SpectralImage:UnknownROImode','Unknown argument %s.',varargin{1});
                end
            else
                if size(varargin{1},1) == obj.height && size(varargin{1},2) == obj.width
                    
                    %For efficiency, find the smallest ROI to load the image.
                    tVert = any(varargin{1},2);
                    top = find(tVert,1,'first');
                    height = find(tVert,1,'last') - top;
                    
                    tHorz = any(varargin{1},1);
                    left = find(tHorz,1,'first');
                    width = find(tHorz,1,'last') - left;
                    
                    ROI = [left, top, width, height];
                      
                    mask = varargin{1}(top:top + height - 1, left: left+width - 1);
                    mask = mask > 0;
                    
                elseif numel(varargin{1}) == 4
                    ROI = varargin{1};
                    
                else
                    
                    error('SpectralImage:MaskNotSameSize','The image is %d-by-%d but the input mask is %d-by-%d',...
                        obj.height, obj.width,...
                        size(varargin{1},1), size(varargin{1},2));
                end
                
            end
            
            storeSpectra = zeros(1,obj.sizeC);
            
            %Get the spectrum
            for iC = 1:obj.sizeC
                
                currSpecImg = obj.getPlane(1,iC,iT,'ROI',ROI);
                
                if exist('mask','var')
                    storeSpectra(iC) = sum(currSpecImg(mask));
                else
                    storeSpectra(iC) = sum(currSpecImg(:));
                end
            end

            spectrumOut = storeSpectra;
            
           
        end
        
        function imgData = getCStack(obj,iT,varargin)
            %GETCSTACK  Get a combined image of all the wavelength channels
            
            ip = inputParser;
            ip.addRequired('iT',@(x) isscalar(x) && isnumeric(x));
            ip.addParameter('ROI',[]);
            ip.parse(iT,varargin{:});
            
            if isempty(ip.Results.ROI)
                ROI = [1, 1, obj.width, obj.height];
            else
                ROI = ip.Results.ROI;
            end
            
            %Get the zstack
            imgData = zeros(ROI(4),ROI(3),obj.sizeC);
            for iC = 1:obj.sizeC
                imgData(:,:,iC) = obj.getPlane(1,iC,ip.Results.iT,'ROI',ROI);
            end
            imgData = sum(imgData,3);
            
        end
        
        function cellMask = segmentCells(obj, iT, thLvl, varargin)
            %SEGMENTCELLS  Segments and displays cells
            
            ip = inputParser;
            ip.addRequired('iT',@(x) isscalar(x) && isnumeric(x));
            ip.addParameter('Watershed', false, @(x) islogical(x));
            ip.addParameter('ROI',[]);
            ip.parse(iT,varargin{:});
            
            %Sum the channels to get a clearer image of the cells
            imgData = obj.getCStack(iT,'ROI',ip.Results.ROI);
            
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
%                 cellLabels = imclearborder(cellLabels);
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
        
        function ratioImg = getRatioImage(obj, iT, iC1, iC2)
            %GETRATIOIMAGE  Get ratiometric image
            
            img1 = double(obj.getPlane(1, iC1, iT)) + 1;
            img2 = double(obj.getPlane(1, iC2, iT)) + 1;
            
            ratioImg = img1 ./ img2;
            
            %Find the maximum value that is not Inf
%             tmp = ratioImg(~isnan(ratioImg) & ~isinf(ratioImg));

            ratioImg(ratioImg > 100) = nan;
            ratioImg(ratioImg > 100) = nan;

%             maxVal = max(tmp(:));
%             
%             ratioImg(isinf(ratioImg)) = maxVal;
%             ratioImg(isnan(ratioImg)) = 0;
        end
        
        function exportImages(obj,iC, frameRange, outputDir,varargin)
            %EXPORTIMAGES  Export images as a Tiff stack
            
            ip = inputParser;
            ip.addParameter('ExportAsStack',true);
            ip.parse(varargin{:});
            
            if ~exist(outputDir,'dir')
                mkdir(outputDir)
            end
            
            [~,fname] = fileparts(obj.filename);
            
            for iT = frameRange
                currImg = normalizeimg(obj.getPlane(1, iC, iT));
                
                if ip.Results.ExportAsStack
                    if iT == frameRange(1)
                        imwrite(currImg,fullfile(outputDir,sprintf('%s (series %d).tif',fname, obj.series)));
                    else
                        imwrite(currImg,fullfile(outputDir,sprintf('%s (series %d).tif',fname, obj.series)),'WriteMode','append');
                    end
                else
                    imwrite(currImg,fullfile(outputDir,sprintf('%s (series %d)_T%02d.tif',fname, obj.series,iT)));
                end
                                
            end
            
            
            
        end        
        
        function obj = loadMasks(obj, maskDir)
            %LOADMASKS  Load PNG files of masks, assuming the GIMP process
            %
            %  Masks should be ordered in a folder, starting with 001.png
            %  being the LAST frame. The masks should be interlaced with
            %  the frame, e.g. 001.png would the mask of the last frame.
            
            %Get filelist
            maskList = dir(fullfile(maskDir,'*.png')); %#ok<*PROPLC>
            maskList = {maskList.name};
            
            %Make sure list is sorted
            maskList = sort(maskList);
            
            %Skip every second, then flip the maskList
            maskList(2:2:end) = [];
            maskList = fliplr(maskList);
            
            obj.maskList = strcat([maskDir,'/'],maskList);
            
%             %Test
%             for ii = 1:numel(maskList)
%                 imgData = imread(fullfile(maskDir,maskList{ii}));
%                 imshow(imgData,[])
%             end
%             
%             keyboard
            
        end
        
        function trackData = trackCells(obj)
            
            if isempty(obj.maskList)
                error('Need masks');
            end
            
            %Get number of masks (e.g. frame range)
            frameRange = 1:numel(obj.maskList);
            
            for iF = frameRange
                
                currMask = imread(obj.maskList{iF});
                currMask = currMask > 0;
                
                currData = regionprops(currMask,'Area','Centroid','PixelIdxList');
                
%                 %Reformat for tracker
%                 cellAreas = cat(1,currData.Area);
%                 cellPositions = cat(1,currData.Centroid);
%                 
%                 avgIntensities = zeros(numel(currData),obj.sizeC);
                
                for iC = 1:obj.sizeC
                    currImg = double(obj.getPlane(1,iC,iF));
                    
                    for iDetObj = 1:numel(currData)
                        %Load all the channels of interest, calculate the
                        %average intensity
                        currData(iDetObj).AvgIntensities(iC) = mean(currImg(currData(iDetObj).PixelIdxList));
                    end
                end
                
                if iF == frameRange(1)
                    %Initialize the tracker
                    trackerObj = LAPtracker(iF,currData);
                    
                else
                    %Track cells
                    trackerObj = trackerObj.assignToTrack(iF,currData);
                    
                end
                
            end
            trackData = trackerObj;
%             trackData = trackerObj.tracks;
            
        end
        
    end
        
end