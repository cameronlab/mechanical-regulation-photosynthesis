classdef MECDataAnalyzer < handle
    %MECDATANALYZER  Analysis class for the mechanical experiment
    %
    %  A = MECDATANALYZER will create a new analyzer object.
    %
    %
      
    properties (Access = private)
        
        trackArray
        cellData
        ColonyAncestorCells
        
    end
    
    properties (Dependent)
        
        NumTracks  %Number of tracks in the track array object
        TrackDataCols %Fieldnames of the tracked data
        
        NumCellRecords  %Number of records in the single cell data table
        CellDataCols   %Names of the cell data table columns
        
        Channels
        
    end
    
    properties (SetAccess = private)
        
        isSpectralData
        
        MeanDeltaT  %Mean time between frames
        TimeUnit    %Time unit
        
        PxSize      %Pixel size in physical units
        PxSizeUnit  %Size unit
        
    end
        
    methods
        
        function obj = MECDataAnalyzer(varargin)
            %Constructor function
            
            if isempty(varargin)
                
                [fname, fdir] = uigetfile({'*.mat','*.mat (MAT-file)'},...
                    'Select tracked data file');
                
                if ~ischar(fname)
                    %Cancelled
                    return;
                else
                    filename = fullfile(fdir, fname);                    
                end
                                
            else
                filename = varargin{1};
                
            end
            
            obj.importData(filename);
            
        end  %Constructor
        
        %--- Get/Set functions
        
        function numTracks = get.NumTracks(obj)
            
            numTracks = numel(obj.trackArray);
            
        end
        
        function trackCols = get.TrackDataCols(obj)
            
            trackCols = obj.trackArray.TrackedDataFields;
            
        end
        
        function numCellRecords = get.NumCellRecords(obj)
            
            numCellRecords = numel(obj.cellData);
                        
        end
        
        function cellCols = get.CellDataCols(obj)
            
            cellCols = fieldnames(obj.cellData);
            
        end
        
        function chNames = get.Channels(obj)
            
            %Return the channel names
            idx = find(startsWith(obj.TrackDataCols,'TotalInt'));
            chNames = cell(1,numel(idx));
            
            for ii = 1:numel(idx)
                chNames{ii} = strtok(obj.TrackDataCols{idx(ii)},'TotalInt');
            end
            
        end
        
        %--- Track functions
        
        function tracksOut = getTrack(obj, trackIdxs)
            %GETTRACK  Get specific tracks
            %
            %  T = A.GETTRACK(I) returns the track(s) specified by I. The
            %  variable T will be an array of TrackData objects.            
            
            tracksOut = obj.trackArray.getTrack(trackIdxs);
            
        end
        
        function trackDataOut = getTrackData(obj, trackIdx, reqProp)
            %GETTRACKDATA  Get track data
            %
            %  D = GETTRACKDATA(A, I, Property) will return the data in the
            %  property specified in track I.
            
            trackDataOut = obj.trackArray.getTrack(trackIdx).getData(reqProp);
            
        end
        
        function arrayOut = getTrackArray(obj)
            %GETTRACK  Get the track array object
            %
            %  T = A.GETTTRACKARRAY returns the track array object T.         
            
            arrayOut = obj.trackArray;
            
        end
        
        %--- Data table functions
        
        function importData(obj,inputData)
            %IMPORTDATA  Import data from file
            %
            %  A.IMPORTDATA(filename) will import the data from the file
            %  specified. The file should contain a TrackDataArray object.
            
            %Parse the input
            if ischar(inputData)
                %Input is a filename
                
                if ~exist(inputData,'file')
                    error('MECDataAnalyzer:importData:FileNotFound',...
                        'File %s not found.', inputData);
                end
                
                data = load(inputData,'trackArray');
                
                %Validate that the variable was loaded as a TrackDataArray
                %object
                if ~isa(data.trackArray,'TrackDataArray')
                    error('MECDataAnalyzer:importData:ErrorImportingClass',...
                        'Expected data to be a ''TrackDataArray'' object, instead it is a %s.',...
                        class(data.trackArray));
                end
                
            elseif isa(inputData, 'TrackDataArray')
                
                data.trackArray = inputData;
                
            else
                error('MECDataAnalyzer:importData:UnknownInputType',...
                        'Unknown input type. Expected a filename or TrackDataArray object.');
            end
                            
            %If the timestamp information was not populated, prompt the
            %user to enter a time between frames in mins
            if isempty(data.trackArray.FileMetadata.Timestamps)
                s = input('Enter time between frames in minutes: ');
                
                if ~isempty(s)
                    %Update the track array (for saving purposes)
                    data.trackArray = data.trackArray.setTimestampInfo(s,'mins');
                else
                    error('MECDataAnalyzer:importData:NoTimeData',...
                        'Time data must be supplied.');
                end
                
            end
            
            %Load the data into the trackArray property
            obj.trackArray = data.trackArray;
            
            %Set the time and pixel size properties
            obj.MeanDeltaT = obj.trackArray.MeanDeltaT;
            obj.TimeUnit = obj.trackArray.FileMetadata.TimestampUnit;
            
            [obj.PxSize, obj.PxSizeUnit] = obj.trackArray.getPxSizeInfo;
            
            %Check if data is spectral data. The spectral data channels are
            %numbers (e.g. 650), while the other movie has letters (e.g.
            %RFP)
            obj.isSpectralData = ~any(ismember(obj.trackArray.TrackedDataFields, 'TotalIntRFP'));
            
            if obj.isSpectralData
                %Analyze the track array
                obj.analyzeSpectralTrackArray;
            else
                %Analyze the track array
                obj.analyzeTrackArray;
            end
        end
        
        function recordOut = getRecord(obj, tableName, recIdxs)
            %GETRECORD  Get rows from the specified data table
            %
            %  R = A.GETRECORD(T, I) will return rows I from the data table
            %  T. T can either be 'cell' for the single cell data table or
            %  'colony' for the colony data table.
            %
            %  R will be a struct containing the rows requested.
            
            
            switch lower(tableName)
                
                case 'cell'
                    
                    recordOut = obj.singleCellDT.getRecord(recIdxs);
                    
                case 'colony'
                    
                    
                otherwise
                    
                    error('MECDataAnalyzer:getRecord:UnknownTableName',...
                        'Table name should be either ''cell'' or ''colony''.');
            
            end
            
        end
        
        function tableOut = getDataTable(obj)
            %GETDATATABLE  Get the data table
            
            tableOut = obj.cellData;
            
        end
        
        function [dataOut, isValid] = queryData(obj, reqData, varargin)
            %QUERYDATA  Returns requested data from the data column
            %
            %  D = DT.QUERYDATA('NuclDiff') returns all the nuclear difference
            %  D = DT.QUERYDATA('TrackIdx', 'NuclDiff > 300') returns all the
            %  track indices for nuclear differences > 300
            
            %Make sure that the requested data is a field in the Data
            %property
            if ~isfield(obj.cellData,reqData)
                error('MECCellDataTable:get:InvalidRequest',...
                    '%s is not a data column. Available data are: %s', ...
                    reqData, strjoin(obj.DataCols,', '));
            end
            
            if strcmp(reqData,'DaughterIdxs')
                %These columns have uneven number of columns
                dataOut = {obj.cellData.(reqData)};
            else
                dataOut = cat(1,obj.cellData.(reqData));
            end
            
            %Parse varargin if it is populated
            if ~isempty(varargin)
                
                %Split string into operations at AND and OR operator
                [expList, opList] = strsplit(varargin{1}, {'&','|','AND','OR'});
                
                %Parse each expression
                for iExp = 1:numel(expList)
                    
                    %Tokens will be {data column} {Operator} {value}
                    strToParse = regexprep(expList{iExp},'\s','');
                    
                    [propName, matches] = strsplit(strToParse,{'>','<','=','~'});
                    
                    tokens = regexp(strToParse,'(\w*)([><=~]{1,2})([0-9]*)','tokens');
                    
                    if isempty(propName)
                        error('Invalid string %s',strToParse);
                    end
                    
                    comparisonStr = [matches{1},propName{2}];
                    
                    %Concatenate the data
                    filtCol = cat(1, obj.cellData.(propName{1}));
                    
                    %Check the relation
                    if ~exist('isValid','var')
                        isValid = eval(sprintf('filtCol %s',comparisonStr));
                    else
                        %Run the comparison depending on the operator
                        switch opList{iExp - 1}
                            
                            case {'&', 'AND'}
                                isValid = isValid & ...
                                    eval(sprintf('filtCol %s',comparisonStr));
                                
                            case {'|', 'OR'}
                                isValid = isValid | ...
                                    eval(sprintf('filtCol %s',comparisonStr));
                                
                        end
                    end
                end
                
                %Remove rows which are not valid
                dataOut(~isValid) = [];
            end
            
        end
        
        %--- Colony analysis functions
        function [timeVec, cellCnts, meanCy5] = getColonyData(obj, colonyID)
            %GETCOLONYDATA  Get colony intensity and cell counts
            %
            %  [T, C, I] = GETCOLONYDATA(A, CID) will return the time
            %  series data containing the cell counts C and the mean Cy5
            %  intensity I of the colony labelled CID. T is the timestamp
            %  vector.
            
            if colonyID <= 0 || colonyID > numel(obj.ColonyAncestorCells)
                error('MECDataAnalyzer:getColonyCellCount:InvalidColonyID',...
                    'Expected colony ID to be between 1 and %d', numel(obj.ColonyAncestorCells))
            end
            
            %Get track IDs for the cells in the specified colony
            cellIDs = [obj.cellData([obj.cellData.ColonyIdx] == colonyID).TrackIdx];
            
            cellCnts = zeros(obj.trackArray.NumFrames,1);
            sumCy5 = zeros(obj.trackArray.NumFrames,1);
            totalArea = zeros(obj.trackArray.NumFrames,1);
            
            timeVec = (1:obj.trackArray.NumFrames) .* obj.MeanDeltaT;
            
            if strcmpi(obj.TimeUnit,'s')
                timeVec = timeVec ./ 60;
            end
            
            for iT = cellIDs
                currTrack = obj.trackArray.getTrack(iT);
                
                
                currCy5Data = currTrack.getData('TotalIntCy5');
                currAreaData = currTrack.getData('Area');
                
                missingIdx = find(isnan(currCy5Data));
                
                %Look for missing values (NaNs)
                tt = currTrack.FirstFrame:currTrack.LastFrame;
                
                tt(missingIdx) = [];
                currCy5Data(missingIdx) = [];
                currAreaData(missingIdx) = [];
                
                                
                cellCnts(tt) = ...
                    cellCnts(tt) + 1;
                
                sumCy5(tt) = ...
                    sumCy5(tt) + currCy5Data;
                
                totalArea(tt) = ...
                    totalArea(tt) + currAreaData;
                
            end
            
            meanCy5 = sumCy5./totalArea;
            
        end
        
        
        %--- Plotting functions
        
        function plotTrack(obj, trackIdxList, trackProp, varargin)
            %PLOTTRACK  Plot data from specified track(s)
            %
            %  PLOTTRACK(A, I, P) will plot data property P from track(s)
            %  I. I can either be a single number of a vector containing
            %  the indices of the tracks to plot.
            %
            %  PLOTTRACK(A, I, 'MajorAxisLength','plotfit') will plot the
            %  cell length and its fitted growth rate curve.
                        
            %Check that the requested property is in the track data
            if ~ismember(trackProp, obj.getTrack(1).TrackDataProps)
                error('MECDataAnalyzer:plotTrack:InvalidProperty',...
                    'Invalid property. Available properties are: %s',...
                    strjoin(obj.getTrack(1).TrackDataProps,', '));
            end
            
            %Plotting options (default)
            plotFittedCurve = false;
            tUnits = 'Frames';
            
            %Parse the variable inputs
            while ~isempty(varargin)
                
                switch lower(varargin{1})
                    
                    case 'plotfit'
                        plotFittedCurve = true;
                end
                
                varargin(1) = [];
            end
            
                                   
            %Gather track data            
            tracks = obj.getTrack(trackIdxList);
            channelNames = cell(1, numel(trackIdxList));    
            
            for iT = 1:numel(tracks)
                
                currData = tracks(iT).getData(trackProp);
                
                %Make the time vector for plotting
                tt = tracks(iT).FirstFrame:tracks(iT).LastFrame;
                
                if ~isempty(obj.trackArray.MeanDeltaT)
                    tt = tt .* obj.trackArray.MeanDeltaT;
                end
                
                %Convert time to minutes if it is in seconds
                if strcmpi(obj.TimeUnit,'s')
                    tt = tt/60;
                    tUnits = 'mins';
                else
                    tUnits = obj.TimeUnit;
                end
                
                if plotFittedCurve && strcmpi(trackProp,'MajorAxisLength')
                    %Plot the growth curve fit
                    
                    fitTime = (1:numel(currData)) .* obj.MeanDeltaT;
                    
                    if strcmpi(obj.TimeUnit,'s')
                       fitTime = fitTime ./ 60; 
                    end
                    
                    %Get fitting parameters
                    gr = obj.cellData(trackIdxList(iT)).GrowthRate;
                    yint = obj.cellData(trackIdxList(iT)).FitYIntercept;
                    
                    %Plot
                    plot(tt,currData .* obj.PxSize,'o',tt,exp(gr .* (fitTime) + yint))
                    
                    legend('Data','Fit','Location','best')
                else
                    %Plot everything else
                    plot(tt, currData);
                    channelNames{iT} = sprintf('Track %d',trackIdxList(iT));
                end

                hold on
                
            end
            hold off
            
            xlabel(tUnits)
            ylabel(sprintf('%s (%s)',trackProp, obj.PxSizeUnit))
            
            if ~(plotFittedCurve && strcmpi(trackProp,'MajorAxisLength'))
                legend(channelNames,'Location','best')
            end
            
        end
        
        function plotRatio(obj, trackIdxList, chNumer, chDenom)
            %PLOTRATIO  Plots ratio of two channels
            %
            %  A.PLOTRATIO(I, numerator, denominator)
            %
            %  Examples:
            %    A.PLOTRATIO(1, 'RFP', 'Cy5')
            %    A.PLOTRATIO(1, 651, 684)
            
            channelNames = cell(1, numel(trackIdxList));            
            
            for iT = 1:numel(trackIdxList)
                currTrack = obj.trackArray.getTrack(trackIdxList(iT));
                
                %Make the time vector for plotting
                tt = currTrack.FirstFrame:currTrack.LastFrame;
                
                if ~isempty(obj.trackArray.MeanDeltaT)
                    tt = tt .* obj.trackArray.MeanDeltaT;
                end
                
                %Convert time to minutes if it is in seconds
                if strcmpi(obj.TimeUnit,'s')
                    tt = tt/60;
                    tUnits = 'mins';
                else
                    tUnits = obj.TimeUnit;
                end
                
                if ischar(chNumer)
                
                    dataNumer = currTrack.getData(sprintf('TotalInt%s',chNumer));
                    dataDenom = currTrack.getData(sprintf('TotalInt%s',chDenom));
                    ylblTxt = sprintf('Ratio %s nm/%s nm',chNumer, chDenom);
                    
                else
                    
                    dataNumer = currTrack.getData(sprintf('TotalInt%d',chNumer));
                    dataDenom = currTrack.getData(sprintf('TotalInt%d',chDenom));
                    ylblTxt = sprintf('Ratio %d nm/%d nm',chNumer, chDenom);
                end
                
                channelNames{iT} = sprintf('Track %d',trackIdxList(iT));
                
                plot(tt,dataNumer./dataDenom)
                hold on
            end
            hold off
            
            legend(channelNames,'Location','best')
            xlabel(tUnits)
            ylabel(ylblTxt)
            
        end
        
        function analyzeTrackArray(obj)
            %ANALYZETRACKARRAY  Analyze the track array data            
            
            obj.ColonyAncestorCells = obj.labelColonies(obj.trackArray);
            
            %Initialize the data array
            obj.cellData = struct('TrackIdx', NaN,...
                'TrackLength', NaN,...
                'DidDivide', false, ...
                'MotherIdx', NaN,...
                'DaughterIdxs', NaN, ...
                'Generation', NaN,...
                'ColonyIdx', NaN,...
                'GrowthRate', NaN, ...
                'FitYIntercept', NaN, ...
                'GrowthRateFitRes', NaN, ...
                'SumMeanIntCy5', NaN, ...
                'MaxMeanIntCy5', NaN, ...
                'MaxHighIntCy5', NaN, ...
                'SumMeanIntRFP', NaN, ...
                'MaxMeanIntRFP', NaN, ...
                'MaxHighIntRFP', NaN);
            
            obj.cellData(obj.NumTracks).TrackIdx = NaN;
            
            %Get data for each track
            for iTrack = 1:obj.trackArray.NumTracks
 
                currTrack = obj.trackArray.getTrack(iTrack);
                
                obj.cellData(iTrack).TrackIdx = iTrack;
                obj.cellData(iTrack).TrackLength = currTrack.NumFrames;
                
                obj.cellData(iTrack).DidDivide = ~any(isempty(currTrack.DaughterIdxs));
                obj.cellData(iTrack).MotherIdx = currTrack.MotherIdx;
                obj.cellData(iTrack).DaughterIdxs = currTrack.DaughterIdxs;
                
                tVec = (1:currTrack.NumFrames) .* obj.MeanDeltaT;
                
                if strcmpi(obj.TimeUnit,'s')
                    tVec = tVec/60;
                end
                
                [obj.cellData(iTrack).GrowthRate, ...
                    obj.cellData(iTrack).FitYIntercept, ...
                    obj.cellData(iTrack).GrowthRateFitRes] = obj.fitGrowthRate(tVec,currTrack.getData('MajorAxisLength') .* obj.PxSize);
                
                obj.cellData(iTrack).SumMeanIntCy5 = sum(currTrack.getData('TotalIntCy5')./currTrack.getData('Area')) ./ currTrack.NumFrames;
                obj.cellData(iTrack).MaxMeanIntCy5 = max(currTrack.getData('TotalIntCy5')./currTrack.getData('Area'));
                obj.cellData(iTrack).MaxHighIntCy5 = max(currTrack.getData('HighIntCy5'));
                
                obj.cellData(iTrack).SumMeanIntRFP = sum(currTrack.getData('TotalIntRFP')./currTrack.getData('Area')) ./ currTrack.NumFrames;
                obj.cellData(iTrack).MaxMeanIntRFP = max(currTrack.getData('TotalIntRFP')./currTrack.getData('Area'));
                obj.cellData(iTrack).MaxHighIntRFP = max(currTrack.getData('HighIntRFP'));

                [obj.cellData(iTrack).Generation, ...
                    obj.cellData(iTrack).ColonyIdx] = obj.getGenAndCol(iTrack);
                
                %If the track exists in the first frame, then classify it
                %as an existing cell
                
            end
            
        end
                
        function analyzeSpectralTrackArray(obj)
            %ANALYZETRACKARRAY  Analyze the track array data            
            
            obj.ColonyAncestorCells = obj.labelColonies(obj.trackArray);
            
            %Initialize the data array
            obj.cellData = struct('TrackIdx', NaN,...
                'TrackLength', NaN,...
                'DidDivide', false, ...
                'MotherIdx', NaN,...
                'DaughterIdxs', NaN, ...
                'Generation', NaN,...
                'ColonyIdx', NaN,...
                'GrowthRate', NaN, ...
                'FitYIntercept', NaN, ...
                'GrowthRateFitRes', NaN);
            
            obj.cellData(obj.NumTracks).TrackIdx = NaN;
            
            %Get data for each track
            for iTrack = 1:obj.trackArray.NumTracks
 
                currTrack = obj.trackArray.getTrack(iTrack);
                
                obj.cellData(iTrack).TrackIdx = iTrack;
                obj.cellData(iTrack).TrackLength = currTrack.NumFrames;
                
                obj.cellData(iTrack).DidDivide = ~any(isempty(currTrack.DaughterIdxs));
                obj.cellData(iTrack).MotherIdx = currTrack.MotherIdx;
                obj.cellData(iTrack).DaughterIdxs = currTrack.DaughterIdxs;
                
                tVec = (1:currTrack.NumFrames) .* obj.MeanDeltaT;
                
                if strcmpi(obj.TimeUnit,'s')
                    tVec = tVec/60;
                end
                
                [obj.cellData(iTrack).GrowthRate, ...
                    obj.cellData(iTrack).FitYIntercept, ...
                    obj.cellData(iTrack).GrowthRateFitRes] = obj.fitGrowthRate(tVec,currTrack.getData('MajorAxisLength') .* obj.PxSize);

                [obj.cellData(iTrack).Generation, ...
                    obj.cellData(iTrack).ColonyIdx] = obj.getGenAndCol(iTrack);
                
                %If the track exists in the first frame, then classify it
                %as an existing cell
                
            end
            
        end
    end

    methods (Access = private)
        
        function [genOut, colonyIdx] = getGenAndCol(obj, trackIdx)
            %GETGENANDCOL  Set the generation number
            %
            %  This algorithm works by searching for the specified track
            %  index to see if it can find a matching entry in the daughter
            %  idx fields.
            %
            %  To speed up the search, we only want to look at preceeding
            %  records, which is why the recIdx must be specified.
            
            isMatch = false;
            colonyIdx = NaN;
            for iRec = 1:trackIdx
                if ismember(trackIdx, obj.cellData(iRec).DaughterIdxs)
                    isMatch = true;
                    break;
                end
            end
            
            if isMatch
                genOut = obj.cellData(iRec).Generation + 1;
                colonyIdx = obj.cellData(iRec).ColonyIdx;
            else
                genOut = 1;
                
                %Check if it is one of the ancestor cells
                for ii = 1:numel(obj.ColonyAncestorCells)
                    if ismember(trackIdx, obj.ColonyAncestorCells{ii})
                        colonyIdx = ii;
                        break;
                    end
                end

            end
            
            
        end
        
    end
    
    methods (Static)
        
        function ancestorCells = labelColonies(trackArrayIn)
            %LABELCOLONIES  Label the colonies from the first frame
            %
            %  L = LABELCOLONIES(A) will return a cell of pixel indices for
            %  cells belonging to the same colony, as determined by the
            %  information from the first frame in track array A.
            %
            %  In particular, cells which touch each other are considered
            %  as belonging to the same colony.
            
            %Create a labelled image of the first frame. Tracks are created
            %based on time, so stop searching after the first frame.
            iT = 1;
            LL = nan(trackArrayIn.FileMetadata.ImgSize);
            while trackArrayIn.getTrack(iT).FirstFrame == 1
                LL(trackArrayIn.getTrack(iT).getData('PixelIdxList','first')) = iT;
                iT = iT + 1;
            end
            
            %Mask the label, then dilate to join neighboring cells
            MM = ~isnan(LL);
            MM = imdilate(MM, strel('disk',2));
            
            cc = bwconncomp(MM,8);
            
            ancestorCells = cell(1,cc.NumObjects);
            for iColony = 1:cc.NumObjects
                ancestorCells{iColony} = unique(LL(cc.PixelIdxList{iColony}));
                ancestorCells{iColony}(isnan(ancestorCells{iColony})) = [];
            end
            
        end
        
        function [growthRate, fitYintercept, fitRes] = fitGrowthRate(tt, cellLengthvTime)
            %FITGROWTHRATE  Calculate the growth rate
            %
            %  The estimated growth rate is fitted to the (natural) log of
            %  the cell length. The fitting uses polyval.
            %
            %  If the cell length data contains NaN values, they will be
            %  removed before fitting. This is because polyval will only
            %  return NaN values if the data is NaN.
            %
            %  Reference:
            %    http://www.sciencedirect.com/science/article/pii/S0092867414014998
            
            logCL = log(cellLengthvTime);
            
            %Remove points where the data is nan - this causes polyfit to
            %always return nans
            delInd = ~isfinite(logCL);
            logCL(delInd) = [];
            tt(delInd) = [];
            
            
            %If there is only 1 point, the code cannot fit a line to it
            if numel(cellLengthvTime) < 2
                growthRate = NaN;
                fitYintercept = NaN;
                fitRes = Inf;
            else
                [fitParams, S] = polyfit(tt',logCL,1);
                
                %The growth rate is the slope
                growthRate = fitParams(1);
                fitYintercept = fitParams(2);
                fitRes = S.normr;
            end
        end
    end
        
end