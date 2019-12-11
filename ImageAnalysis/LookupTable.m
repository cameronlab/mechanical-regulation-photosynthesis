classdef LookupTable
    %LookupTable information
    
    properties
       
        LUT8Bit
        LUT16Bit
        
    end
    
    properties (Dependent)
        
        bitRange
        
    end
    
    methods
        
        function obj = LookupTable(varargin)
            
            if ~isempty(varargin)
                
                if isa(varargin{1},'loci.formats.ChannelSeparator')
                    obj.LUT16Bit = varargin{1}.get16BitLookupTable;
                    obj.LUT8Bit = varargin{1}.get8BitLookupTable;
                elseif isa(varargin{1},'BioformatsImage')
                    obj.LUT16Bit = varargin{1}.bfReader.get16BitLookupTable;
                    obj.LUT8Bit = varargin{1}.bfReader.get8BitLookupTable;
                    
                end
                
            end
            
        end
       
        function imgOut = adjustImage(obj, imgIn,varargin)
            
            switch obj.bitRange
                
                case 8
                    
                    lut = obj.LUT8Bit;
                    
                case 16
                    lut = obj.LUT16Bit;
                otherwise
                
                
            end
            
            %Apply LUT to each channel
            %Format of LUT is: N x bitDepth
            imgOut = zeros(size(imgIn),class(imgIn));
            for ii = 1:size(imgIn,3)
                currLUT = lut(ii,:);
                imgOut(:,:,ii) = currLUT(1 + imgIn);
            end
            
        end
        
        function bitRange = get.bitRange(obj)
            
            if isempty(obj.LUT8Bit) && ~isempty(obj.LUT16Bit)
                bitRange = 16;
            elseif ~isempty(obj.LUT8Bit) && isempty(obj.LUT16Bit)
                bitRange = 8;
            else
                bitRange = nan;   
            end
            
        end
        
    end
    
end