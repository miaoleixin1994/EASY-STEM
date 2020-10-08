%serReader() function
%	Description:
%		File reader for .ser files from FEI's TIA or Emispec acquisition programs for electron microscopes. Reads in data from .SER files and the data is divided into "elements" such as line scans, spectra or images
%	Parameters:
%		fname - the filename of the SER file to be read. If not provided the program opens a dialog box to choose a file
%	Limitations:
%		Complex datatypes are currenty unsupported and will cause the program to fail
%	Output:
% 		Ser - a matlab structure that contains the data and the pixel sizes for each data element
%	Author:
%		Peter Ercius 2009/08/06
%
%---------------------- NO WARRANTY ------------------
%THIS PROGRAM IS PROVIDED AS-IS WITH ABSOLUTELY NO WARRANTY
%OR GUARANTEE OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
%INCLUDING BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
%MERCHANABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
%IN NO EVENT SHALL THE AUTHOR BE LIABLE
%FOR DAMAGES RESULTING FROM THE USE OR INABILITY TO USE THIS
%PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA
%BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR
%THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH
%ANY OTHER PROGRAM). 
%------------------------------------------------------------------------

function Ser = serReader_forapp(fName)

	fPath = '';

    if nargin == 0
        [fName, fPath] = uigetfile('*.ser');
        fullFileName = fullfile(fPath, fName);
        [~,name,~] = fileparts(fullFileName);
        Ser.fullFileName = fullFileName;
        Ser.name = name;
    end
	
	if fName == 0
		error('No file opened.');
	end
	
    [FID, FIDmessage] = fopen([fPath fName],'rb');
	if FID == -1
		disp(fName)
		error(['Issue opening file: ' FIDmessage])
    end
    
	%%
    fseek(FID,0,'eof');
    totalBytes = ftell(FID);
    fseek(FID,0,-1);
    %%
	%Read in the header information
	
	byteOrder = fread(FID,1,'int16'); %dec2hex(byteOrder) should be 4949
	seriesID = fread(FID,1,'int16');
	seriesVersion = fread(FID,1,'int16');
	dataTypeID = fread(FID,1,'int32'); %0x4120 = 1D arrays; 0x4122 = 2D arrays
	tagTypeID = fread(FID,1,'int32'); %0x4152 = Tag is time only; 0x4142 = Tag is 2D with time
	
	totalNumberElements = fread(FID,1,'int32'); %the number of data elements in the original data set
	validNumberElements = fread(FID,1,'int32'); %the number of data elements written to the file
	offsetArrayOffset = fread(FID,1,'int32'); %the offset (in bytes) to the beginning of the data offset array
	
	numberDimensions = fread(FID,1,'int32'); %the number of dimensions of the indices (not the data)
	
	%%
	%Dimension arrays (starts at byte offset 30)
    for ij = 1:numberDimensions
        dimensionSize = fread(FID,1,'int32'); %the number of elements in this dimension
        calibrationOffset = fread(FID,1,'float64'); %calibration value at element calibrationElement
        calibrationDelta = fread(FID,1,'float64'); %calibration delta between elements of the series
        calibrationElement = fread(FID,1,'int32'); %the element in the series with a calibraion value of calibrationOffset
        descriptionLength = fread(FID,1,'int32'); %length of the description string
        description = fread(FID,descriptionLength,'*char')'; %describes the dimension
        unitsLength = fread(FID,1,'int32'); %length of the units string
        units = fread(FID,unitsLength,'*char'); %name of the units in this dimension
    end
    
	%%
	%Get arrays containing byte offsets for data and tags
    fseek(FID,offsetArrayOffset,-1); %seek in the file to the offset arrays
	%Data offset array (the byte offsets of the individual data elements)
	dataOffsetArray = fread(FID,totalNumberElements,'int32');
	%Tag Offset Array (byte offsets of the individual data tags)
	tagOffsetArray = fread(FID,totalNumberElements,'int32');
	
	%%
	%Get data from 1D elements
	if dataTypeID == hex2dec('4120')
		%Ser.calibration = zeros(1, validNumberElements);
		for ii=1:validNumberElements
			fseek(FID,dataOffsetArray(ii),-1);
			calibrationOffset = fread(FID,1,'float64'); %calibration value at element calibrationElement
			calibrationDelta = fread(FID,1,'float64'); %calibration delta between elements of the array
			calibrationElement = fread(FID,1,'int32'); %element in the array with calibration value of calibrationOffset
			dataType = fread(FID,1,'int16'); %
            
			Type = getType(dataType);
			arrayLength = fread(FID,1,'int32'); %number of elements in the array
            
			Ser.data{ii} = fread(FID,arrayLength,Type);
			energy = linspace(calibrationOffset,calibrationOffset+(arrayLength-1)*calibrationDelta,arrayLength);
			Ser.calibration{ii} = energy;
			%Ser.calibration(ii) = length(Ser.data{ii})*calibrationDelta + calibrationOffset;
		end
	%Get data from 2D elements	
	elseif dataTypeID == hex2dec('4122')
		%Ser.calibration = zeros(2, validNumberElements);
		for ii=1:validNumberElements
			fseek(FID,dataOffsetArray(ii),-1);
			calibrationOffsetX = fread(FID,1,'float64'); %calibration at element calibrationElement along x
			calibrationDeltaX = fread(FID,1,'float64');
			calibrationElementX = fread(FID,1,'int32'); %element in the array along x with calibration value of calibrationOffset
			calibrationOffsetY = fread(FID,1,'float64');
			calibrationDeltaY = fread(FID,1,'float64');
			calibrationElementY = fread(FID,1,'int32');
			dataType = fread(FID,1,'int16');
		
			Type = getType(dataType);
			arraySizeX = fread(FID,1,'int32');
			arraySizeY = fread(FID,1,'int32');
            
			Ser.data{ii} = fread(FID,[arraySizeX arraySizeY],Type);
			Ser.calibration(:,ii) = [calibrationDeltaX calibrationDeltaY]';
		end
		
		%Check to see if only one object was valid
		%Remove cell 
		if validNumberElements == 1
			Ser.image = Ser.data{1};
			Ser = rmfield(Ser,'data');
		end

		if validNumberElements == 0
			Ser.image = []; %return empty matrix
			Ser.calibration = []; %return empty matrix
		end
	end
	%%
	%Read in the tags
	%Get data from "Time-only" tags
	if tagTypeID == hex2dec('4152')
		for jj=1:validNumberElements
			fseek(FID,tagOffsetArray(jj),-1);
			tagTypeID2 = fread(FID,1,'int16'); %type of the tag (should be 0x4152 or 16722)
			timeTag = fread(FID,1,'int32'); %number of seconds since Jan. 1, 1970
	end
	%Get data from Time-and-Position tags
	elseif tagTypeID == hex2dec('4142')
        newTags = zeros(3,validNumberElements);
		for jj=1:validNumberElements
			fseek(FID,tagOffsetArray(jj),-1);
			tagTypeID2 = fread(FID,1,'int16');
			newTags(1,jj) = fread(FID,1,'int32'); %time Tag
            newTags(2,jj) = fread(FID,1,'float64'); %xTag
            newTags(3,jj) = fread(FID,1,'float64'); %yTag
            %timeTag = fread(FID,1,'int32');
			%xTag = fread(FID,1,'float');
			%yTag = fread(FID,1,'float');
            %disp(tagTypeID2)
            %disp([num2str(xTag) ' ' num2str(yTag)]);
        end
	else disp(['Unknown time tag type ' num2str(tagTypeID) ' in ' fName])	
	end

	%%
	%Clean up and write out the data
% 	disp(['Total elements read: ' num2str(validNumberElements)])
	fclose(FID);
        
end

function Type = getType(dataType)
	switch dataType
		case 1
			Type = 'uint8';
		case 2
			Type = 'uint16';
		case 3
			Type = 'uint32';
		case 4
			Type = 'int8';
		case 5
			Type = 'int16';
		case 6
			Type = 'int32';
		case 7
			Type = 'float32';
		case 8
			Type = 'float64';
		otherwise
			error('Unsupported data type') %complex data types are unsupported currently
	end

end