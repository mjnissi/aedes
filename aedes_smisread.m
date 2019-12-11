function DATA = aedes_smisread(filename,opt)
% AEDES_SMISREAD - Read S.M.I.S. Data Files (*.sur)
%   
%
% Synopsis: 
%       DATA = aedes_smisread(filename)
%
% Description:
%       Read image data from S.M.I.S. SUR-Files.
%
% Examples:
%
% See also:
%       AEDES_READFID, AEDES_READ_NIFTI

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

% *************************************************************************
% SMIS Data format documentation:
%
% SMIS reconstructed image files are in a format very similar to the first
% part of that used for raw data (.MRD).  They consist of:
% 
% 1)      a 256 byte header
% 2)      a 256 byte text description
% 3)      the data proper
% 4)      the pule program parameters (PPR file)
% 5)      other parameters
% 
% 256 byte header
% Bytes (hex)     'C' data type           Usage
% 0-3             Long (4 byte integer)   Number of pixels in X direction
% 4-7             Long                    Number of pixels in Y direction
% 8-B             Long                    Number of pixels in Z direction (currently 1)
% C-F             Long                    Dimension 4 (currently 1)
% ...                                     Unspecified
% 12-13           Int (2 byte integer)    Data type code
% ...                                     Unspecified
% 30-33           Float   (4 bytes)       Scaling factor from highest pixel to 4095
% ...                                     Unspecified
% 6D              Unsigned char (1 byte)  Number of bits per pixel (currently 12)
% ...                                     Unspecified    
% 
% The data type code indicates the format of individual data elements, as
% follows:
% 
%         0x00    unsigned char           1 byte
%         0x01    signed char             1 byte
%         0x02    short                   2 bytes
%         0x03    int                     2 bytes
%         0x04    long                    4 bytes
%         0x05    float                   4 bytes
%         0x06    double                  8 bytes
% 
%         0x10    bit mask indicating complex data as real, imaginary pairs of
% base data type
% 
% All SMIS .SUR files are at present integer, indicated by data type 0x03.
% 
% 256 text description
% At present this is unused and contains the ASCII text "There is no script
% !".
% 
% Data
% Data follows.  The pixel positions increments in the X dimension until the
% next increment in the Y dimension.
% *************************************************************************

DATA = [];
ReadHeader = true;
ReadData = true;

if nargin == 0 || isempty(filename)
  % Prompt for a file
  [fname,fpath,findex] = uigetfile({'*.sur',...
	'S.M.I.S. Data Files (*.sur)';...
	'*.*','All Files (*.*)'},'Open S.M.I.S. Data File');
  if isequal(fname,0)
	return
  end
  filename = fullfile(fpath,fname);
end
[fp,fn,fe]=fileparts(filename);

if nargin == 2
  % Read only header
  if strcmpi(opt,'header')
	ReadData = false;
  end
end

% Open the file for reading
fid = fopen(filename,'r');
if fid < 0
  error('Could not open file "%s" for reading!',filename)
end

% Read data header
if ReadHeader
  hdr=l_ReadHdrInfo(fid);
end

% Read image data
if ReadData
  [data,msg]=l_ReadData(hdr,fid);
  if ~isempty(msg)
	fclose(fid);
	error(msg);
  end
end

% Close file
fclose(fid);

% Construct DATA structure
DATA.DataFormat = 'SUR';
DATA.FTDATA = data;
DATA.KSPACE = [];
DATA.HDR.FileHeader = hdr;
DATA.HDR.fname = [fn,fe];
DATA.HDR.fpath = [fp,filesep];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Header Information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hdr=l_ReadHdrInfo(file_fid)

% 256 byte header
% Bytes (hex)     'C' data type           Usage
% 0-3             Long (4 byte integer)   Number of pixels in X direction
% 4-7             Long                    Number of pixels in Y direction
% 8-B             Long                    Number of pixels in Z direction (currently 1)
% C-F             Long                    Dimension 4 (currently 1)
% ...                                     Unspecified
% 12-13           Int (2 byte integer)    Data type code
% ...                                     Unspecified
% 30-33           Float   (4 bytes)       Scaling factor from highest pixel to 4095
% ...                                     Unspecified
% 6D              Unsigned char (1 byte)  Number of bits per pixel (currently 12)
% ...                                     Unspecified    
% 
% The data type code indicates the format of individual data elements, as
% follows:
% 
%         0x00    unsigned char           1 byte
%         0x01    signed char             1 byte
%         0x02    short                   2 bytes
%         0x03    int                     2 bytes
%         0x04    long                    4 bytes
%         0x05    float                   4 bytes
%         0x06    double                  8 bytes
% 
%         0x10    bit mask indicating complex data as real, imaginary pairs of
% base data type

hdr = [];

% Seek to the beginning of the file
fseek(file_fid,0,-1);

hdr.nPixelsX=fread(file_fid,1,'int32');	 % Number of pixels in X direction
hdr.nPixelsY=fread(file_fid,1,'int32');	     % Number of pixels in Y direction
hdr.nPixelsZ=fread(file_fid,1,'int32');	     % Number of pixels in Z direction
hdr.nPixels4D=fread(file_fid,1,'int32');	     % Number of pixels in 4 dimension

% Seek over the unspecified bytes
fseek(file_fid,18,-1);

hdr.DataTypeCode = fread(file_fid,1,'int16'); % Data type code

% Seek over the unspecified bytes
fseek(file_fid,48,-1);

hdr.ScalingFactor = fread(file_fid,1,'float'); % Scaling factor from highest 
                                          % pixel to 4095

% Seek over the unspecified bytes
fseek(file_fid,109,-1);

hdr.BitsPerPixel = fread(file_fid,1,'uchar');  % Number of bits per pixel 
                                          % (currently 12)

% Read the text description
fseek(file_fid,256,-1);
hdr.TextDescription = char(fread(file_fid,256,'char')).';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Image Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,msg]=l_ReadData(hdr,file_fid)

data = [];
msg = '';

% Determine data type
%         0x00    unsigned char           1 byte
%         0x01    signed char             1 byte
%         0x02    short                   2 bytes
%         0x03    int                     2 bytes
%         0x04    long                    4 bytes
%         0x05    float                   4 bytes
%         0x06    double                  8 bytes
dataTypes = {...
  '*uint8',1;...
  '*int8',1;...
  '*int16',2;...
  '*int16',2;...
  '*int32',4;...
  '*single',4;...
  '*double',6};
precision = dataTypes{hdr.DataTypeCode+1,1};

% Seek to the beginning of the data
fseek(file_fid,512,-1);

% Size of the data in bytes
data_size = hdr.nPixelsX*hdr.nPixelsY*hdr.nPixelsZ*hdr.nPixels4D;

try
  % Read image data
  data = fread(file_fid,data_size,precision);
  
  % Reshape and permute to correct size
  data = reshape(data,[hdr.nPixelsX hdr.nPixelsY hdr.nPixelsZ hdr.nPixels4D]);
  data = permute(data,[2 1 3 4]);
catch
  msg = ['Could not read data: ',lasterr];
  data = [];
end

