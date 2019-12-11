function DATA = aedes_readfdf(filename,varargin)
% AEDES_READFDF - Read Varian Flexible Data Format (FDF) Files
%   
%
% Synopsis: 
%
% Description:
%
% Examples:
%
% See also:
%

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


if nargin==0 || isempty(filename)
  %% Ask for a file name
  [f_name, f_path, f_index] = uigetfile( ...
    {'*.fdf;*.FDF','Varian FDF-files (*.fdf)'; ...
     '*.*', 'All Files (*.*)'}, ...
    'Select Varian FDF-file');
  if ( all(f_name==0) | all(f_path==0) ) % Cancel is pressed
    DATA=[];
    return
  else
    filename = [f_path,f_name];
  end
end

% Try to open file for reading
file_fid = fopen(filename,'r','ieee-le');
if file_fid<0
  DATA=[];
  error('Could not open file "%s" for reading!',filename);
end

[fp,fn,fe]=fileparts(filename);

% Read Header
[FileHeader,msg]=l_ReadHdr(file_fid);
if isempty(FileHeader)
  DATA=[];
  error(msg);
end

% Check endianness
if isfield(FileHeader,'bigendian') && FileHeader.bigendian==1
  pos = ftell(file_fid);
  fclose(file_fid);
  file_fid = fopen(filename,'r','ieee-be');
  if file_fid<0
    DATA=[];
    error('Could not open file "%s" for reading!',filename);
  end
  
  % Seek to data start
  fseek(file_fid,pos,-1);
end

% Read Data
[data,msg]=l_ReadData(FileHeader,file_fid);

% Close file
fclose(file_fid);

% Construct DATA structure
DATA.DataFormat = 'FDF';
DATA.HDR.FileHeader = FileHeader;
DATA.HDR.fname = [fn,fe];
DATA.HDR.fpath = [fp,filesep];
DATA.FTDATA = data;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Data Header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FileHeader,msg]=l_ReadHdr(file_fid)

FileHeader=[];
msg = '';


%% Check if file is a valid FDF-file
magic_str=char(fread(file_fid,25,'char'))';
if length(magic_str)<25 || ~strcmpi(magic_str(1:end-1),'#!/usr/local/fdf/startup')
  msg='File is not a valid FDF-file!';
  return
end

% Seek the file back to beginning

try
%% Read ASCII header to a cell array
done=false;
count=0;
header={};
while ~done
  count=count+1;
  
  % Check first and second char on the line. 
  % In SWIFT reconstruction fdf:s there can be an additional
  % form feed character (12) before NULL...
  tmp=fread(file_fid,2,'uint8');
  if any(tmp==0)
    if tmp(1)==0
      fseek(file_fid,-1,0);
    end
    done=true;
    continue;
  else
    % Move file indicator back
    fseek(file_fid,-2,0);
  end
  
  header{count}=fgetl(file_fid);
end



%% Parse the header cell array
for ii=1:length(header)
  if length(header{ii})<=2 || all(header{ii}==32) || ...
        all(header{ii}==0)
    continue;
  end
  
  % Remove * and [] chars
  header{ii}=strrep(header{ii},'[','');
  header{ii}=strrep(header{ii},']','');
  header{ii}=strrep(header{ii},'*','');
  
  % Replace " chars with '
  header{ii}=strrep(header{ii},'"','''');
  
  %% Find the first space char
  ind=find(header{ii}==32);
  if length(ind)>1 && diff([ind(1:2)])==1
    eval(['FileHeader.',header{ii}(ind(2)+1:end)])
  else
    eval(['FileHeader.',header{ii}(ind(1)+1:end)])
  end
end

%% Check that numeric fields are numeric
fldnames=fieldnames(FileHeader);
for ii=1:length(fldnames)
  if iscell(FileHeader.(fldnames{ii}))
    if all(cellfun(@isnumeric,FileHeader.(fldnames{ii})))
      tmp=FileHeader.(fldnames{ii});
      FileHeader.(fldnames{ii}) = [tmp{:}];
    end
  end
end

catch
  FileHeader=[];
  msg='Error occurred while reading FDF file header!';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,msg]=l_ReadData(FileHeader,file_fid)

data=[];
msg='';

%% Read data parameters
data_sz = FileHeader.matrix;
datatype = FileHeader.storage;
bits = FileHeader.bits;
dataformat=FileHeader.type;

if strcmpi(datatype,'integer')
  data_str = ['*int',num2str(bits)];
elseif strcmpi(datatype,'float')
  data_str = ['*float',num2str(bits)];
end

%% Read data
[data,count]=fread(file_fid,inf,data_str);
data=reshape(data,data_sz);

% Permute to correct orientation
data=permute(data,[2 1 3 4]);

%for ii=1:size(data,3)
%  data(:,:,ii)=rot90(data(:,:,ii));
%end

