function [DATA,msg] = aedes_readswiftsgl(filename,varargin)
% AEDES_READSWIFTSGL - Read SWIFT reconstruction (*.sgl) files
%   
%
% Synopsis: 
%       [DATA,msg]=aedes_readswiftsgl(filename,varargin)
%
% Description:
%       The function reads SWIFT reconstruction SGL files. If not
%       specified, [256 256 256] data size is assumed. The DATASIZE
%       property can be used to specify the data size, and the INTERPSIZE
%       property can be used to upsample data. The DATAINTERPSIZEPROMPT
%       property can be used to prompt for these data sizes.
%
%       If called without input arguments or if the FILENAME arguments is
%       an empty string, the user is prompted to input file.
%
%       property/value pairs
%
%       'DataSize'   - Size of the data matrix (default: [256 256 256])
%       
%       'InterpSize' - Interpolated size (default: no interpolation, empty)
%
%       'DataInterpSizePrompt' - Prompt for Data/Interp sizes (default: 0)
%
%       'Wbar' - Show/hide waitbar (default: false)
%       
% Examples:
%       [DATA,msg]=aedes_readswiftsgl(filename,'prop1',value1,'prop2',value2,...);
%
% See also:
%       AEDES_READFID, AEDES_READ_NIFTI, AEDES

% This function is a part of Aedes - A graphical tool for analyzing 
% medical images
%
% Copyright (C) 2006 Juha-Pekka Niskanen <Juha-Pekka.Niskanen@uku.fi>
% 
% Department of Physics and Mathematics, Department of Neurobiology
% University of Eastern Finland, FINLAND
%
% This program may be used under the terms of the GNU General Public
% License version 2.0 as published by the Free Software Foundation
% and appearing in the file LICENSE.TXT included in the packaging of
% this program.
%
% This program is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
% WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

% Defaults
DATA = [];
msg = '';
InterpSize = [];  % If empty, don't interpolate
DataSize = [256 256 256];
DataInterpSizePrompt = false;
ShowWbar = false;

% Check if we need to be silent about errors
if nargout>1
  silent = true;
else
  silent = false;
end

%% Parse varargin
for ii=1:2:length(varargin)
  switch lower(varargin{ii})
   case {'interpsize'}
     InterpSize = varargin{ii+1};
     
    case 'datasize'
      DataSize = varargin{ii+1};
      
    case 'datainterpsizeprompt'
      DataInterpSizePrompt = varargin{ii+1};
      
    case 'wbar'
      ShowWbar = varargin{ii+1};
      
   otherwise
    done=false;
    msg=sprintf('Unknown parameter "%s"!',varargin{ii});
    return
  end
end

% Prompt for file if not given as an input argument
if nargin==0 || isempty(filename)
  [fname, fpath, findex] = uigetfile( ...
    {'*.sgl','SWIFT SGL-files (*.sgl)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Select SWIFT SGL file','');
  if isequal(fname,0)
    msg = 'Canceled';
    return
  end
  filename = [fpath,fname];
end

% Get file name and path
[fp,fn,fe]=fileparts(filename);


% Prompt for data/interpolation size
if DataInterpSizePrompt
  done=false;
  while ~done
    answer = inputdlg({'Data size','Interpolated size'},...
      'Input data/interp size',1,...
      {'[256 256 256]','[256 256 256]'});
    if isempty(answer)
      msg = 'Canceled';
      return
    end
    DataSize = str2num(answer{1});
    InterpSize = str2num(answer{2});
    if length(DataSize)~=3 || length(InterpSize)~=3
      h=errordlg('Invalid Data/Interpolation size!',...
        'Invalid size','modal');
      uiwait(h)
    else
      done = true;
    end
  end
  if all(DataSize==InterpSize)
    InterpSize = [];
  end
end

% Check interpolation size
if ~isempty(InterpSize)
  if length(InterpSize)~=3
    msg = 'Invalid interpolation size!';
    if ~silent
      error(msg)
    end
    return
  end
end

% Show waitbar
if ShowWbar
  [wbh,txh] = aedes_calc_wait('Reading data in SWIFT SGL format...');
end

% Open file for reading
fid=fopen(filename,'r');
if fid<0
  msg = sprintf('Could not open file %s for reading',filename);
  if ~silent
    error(msg)
  end
  return
end

% Read data
I=fread(fid,inf,'float32','ieee-be');

% Close file
fclose(fid);

% Reshape to correct size
try
  I=reshape(I,DataSize);
  I=permute(I,[2 1 3]);
catch
  msg = sprintf('Could not reshape data to size [%d %d %d]',...
    DataSize);
  if ~silent
    error(msg)
  end
  return
end

% Interpolate data if requested
if ~isempty(InterpSize)
  
  % Update waitbar
  set(txh,'string',sprintf('Interpolating data to [%d %d %d]',...
    InterpSize))
  drawnow
  
  % Fix the InterpSize to be greater or equal to the DataSize
  for ii=1:3
    InterpSize(ii) = max(InterpSize(ii),DataSize(ii));
  end
  
  if ~all(InterpSize==DataSize)
  
    % Scale to positive
    min_val = min(I(:));
    I = I+min_val;
    
    % Interpolate to InterpSize using Fourier interpolation
    I = fftshift(fftn(I));
    I(InterpSize(1),InterpSize(2),InterpSize(3))=single(0);
    I=abs(ifftn(ifftshift(I)));
    
    % Scale back to original minimum
    I=I-min_val;
  end
end

% Construct DATA structure
DATA.DataFormat = 'swift_sgl';
DATA.HDR.fname = [fn,fe];
DATA.HDR.fpath = [fp,filesep];
DATA.FTDATA = I;

if ShowWbar
  delete(wbh)
end



