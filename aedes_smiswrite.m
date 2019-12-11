function aedes_smiswrite(data,filename,varargin)
% AEDES_SMISWRITE - Write S.M.I.S. Data Files (*.sur)
%   
%
% Synopsis: 
%       aedes_smiswrite(data,filename,param1,value1,param2,value2,...)
%
% Description:
%       Write data to S.M.I.S. SUR-Format files.
%
% Examples:
%
% See also:
%       AEDES_SMISREAD, AEDES_READFID, AEDES_READ_NIFTI

% This function is a part of Aedes - A graphical tool for analyzing 
% medical images
%
% Copyright (C) 2006 Juha-Pekka Niskanen <Juha-Pekka.Niskanen@uku.fi>
% 
% Department of Physics, Department of Neurobiology
% University of Kuopio, FINLAND
%
% This program may be used under the terms of the GNU General Public
% License version 2.0 as published by the Free Software Foundation
% and appearing in the file LICENSE.TXT included in the packaging of
% this program.
%
% This program is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
% WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

% Defaults
DataTypeCode = 3;
BitsPerPixel = 12;
TextDescription = 'Exported from Aedes';

if nargin<1
  error('Too few input arguments!')
elseif nargin==1
  % Prompt for a file
  [fname,fpath]=uiputfile({'*.sur','S.M.I.S. SUR-Files (*.sur)';...
	'*.*','All Files (*.*)'},'Save Data in SUR Format',...
	'surdata.sur');
  if isequal(fname,0)
	% Canceled
	return
  end
  filename = fullfile(fpath,fname);
end
[fp,fn,fe]=fileparts(filename);
if isempty(fe) || ~any(strcmpi(fe,{'.sur','.mri'}))
  fe = '.sur';
end
filename = fullfile(fp,[fn,fe]);

if isstruct(data) && isfield(data,'FTDATA')
  data = data.FTDATA;
end

% Handle varargin
for ii=1:2:length(varargin)
  prop = lower(varargin{ii});
  value = varargin{ii+1};
  switch prop
	case 'datatype'
	  if strcmpi(value,'uint8')
		DataTypeCode = 0;
		BitsPerPixel = 8;
	  elseif strcmpi(value,'int8')
		DataTypeCode = 1;
		BitsPerPixel = 8;
	  elseif any(strcmpi(value, {'uint16','int16'}))
		DataTypeCode = 3;
		BitsPerPixel = 16;
	  elseif strcmpi(value,'int32')
		DataTypeCode = 4;
		BitsPerPixel = 32;
	  elseif any(strcmpi(value,{'float','single'}))
		DataTypeCode = 5;
		BitsPerPixel = 32;
	  elseif any(strcmpi(value,{'float64','double'}))
		DataTypeCode = 6;
		BitsPerPixel = 64;
	  else
		DataTypeCode = 6;
		BitsPerPixel = 64;
	  end
	case 'description'
	  TextDescription = value;
	otherwise
	  error('Unknown property %s',upper(prop))
  end
end

% Get data dimensions
x_size = size(data,2);
y_size = size(data,1);
z_size = size(data,3);
v_size = size(data,4);

% Open file for writing
fid = fopen(filename,'w');
if fid < 0
  error('Could not open file "%s" for writing!',filename)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write header information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allocate 512 bytes from the header
fwrite(fid,zeros(1,512),'char');

% Write size information
fseek(fid,0,-1);
fwrite(fid,x_size,'int32');
fwrite(fid,y_size,'int32');
fwrite(fid,z_size,'int32');
fwrite(fid,v_size,'int32');

% Write DataType
fseek(fid,18,-1);
fwrite(fid,DataTypeCode,'int16');

% Write Bits Per Pixel
fseek(fid,109,-1);
fwrite(fid,BitsPerPixel,'uchar');

% Write Text Description
fseek(fid,256,-1);
fwrite(fid,TextDescription(1:min(256,end)),'char');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Image Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Crop image data to BitsPerPixel if 12-bit
if BitsPerPixel==12
  data(data<0) = 0;
  data(data>4096) = 4096;
end

% Seek to the beginning of the data
fseek(fid,512,-1);

% Determine precision
dataTypes = {...
  'uint8',...
  'int8',...
  'int16',...
  'int16',...
  'int32',...
  'single',...
  'double'};
precision = dataTypes{DataTypeCode+1};

% Write image data
count=fwrite(fid,permute(data,[2 1 3 4]),precision);

% Close file
fclose(fid);

% Check if all elements were written
if count~=prod(size(data))
  error('The file "%s" was not written properly!',filename)
end


% - EOF -
