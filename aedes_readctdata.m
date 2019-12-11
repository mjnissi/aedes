function [DATA,msg] = aedes_readctdata(filename,varargin)
% AEDES_READCTDATA - Read Gamma Medica CT/SPECT data format
%   
% 
% Synopsis: 
%        [DATA,msg] = aedes_readctdata(filename,varargin)
%
% Description:
%
% Examples:
%
% See also:
%        AEDES

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


DATA=[];
msg='';
ShowWbar = true;
rotation = 0;
flipping = '';

if nargin==0 || isempty(filename)
  
  % Ask for a file
  [fname,fpath,findex]=uigetfile({'*.hdr;*.HDR','CT/SPECT header files (*.hdr)';...
                      '*.xxm;*.XXM','Reconstruction parameter file (*.xxm)';...
                      '*.*','All Files (*.*)'},...
                                 'Select a file');
  if isequal(fname,0) || isequal(fpath,0)
    if nargout==2
      msg = 'Cancelled';
    end
    return
  end
else
  [fpath,fname,fext]=fileparts(filename);
  fpath = [fpath,filesep];
  fname = [fname,fext];
  if strcmpi(fext,'.hdr')
    findex = 1;
  elseif strcmpi(fext,'.xxm')
    findex = 2;
  end
end

% Go through varargin
for ii=1:2:length(varargin)
  switch lower(varargin{ii})
   case 'rotate'
    rotation = varargin{ii+1};
   case 'flip'
    flipping = varargin{ii+1};
  end
end

if findex==1
  
  % Read header file
  [hdr,msg] = l_ReadHDRFile([fpath,fname]);
  if ~isempty(msg)
    DATA=[];
    return
  end
  
  if ShowWbar
    % Show aedes_calc_wait
    [h,txh]=aedes_calc_wait('Reading SPECT/CT data...');
    drawnow
  end
  
  % Read data
  [data,msg] = l_ReadImgFile([fpath,hdr.DataFileName],hdr);
  if ~isempty(msg)
    DATA=[];
    return
  end
  
  if iscell(data)
    
    DATA={};
    
    for ii=1:length(data)
      DATA{ii}=struct('DataFormat','spect/ct',...
                      'FTDATA',data{ii},...
                      'HDR',hdr);
    end
  else
    DATA.FTDATA = data;
    DATA.HDR = hdr;
    
    % Include file information
    DATA.HDR.fname = fname;
    DATA.HDR.fpath = fpath;
  end
  
  if ShowWbar
    delete(h)
  end
  
elseif findex==2
  
  % Look for data files in the same directory
  %D=dir(fpath);
  %f_names={D(~[D(:).isdir]).name};
  %slice_ind=strncmpi(fnames,'slice.',6);
  
  % Read parameter file
  [DATA.HDR,msg]=l_ReadXXMFile([fpath,fname]);
  if ~isempty(msg)
    DATA=[];
    return
  end
  
  % Number of slice files
  nFiles = DATA.HDR.PARTAG_CUBESIZEZ;
  
  slice_files = {};
  for ii=1:nFiles
    slice_files{ii}=[fpath,'slice.',sprintf('%04i',ii-1)];
  end
  
% $$$   if all(slice_ind==0)
% $$$     if nargout==2
% $$$       msg = 'Could not find slice data files!'; 
% $$$     else
% $$$       error('Could not find slice data files!');
% $$$     end
% $$$     return
% $$$   end
% $$$   slice_files = {f_names{slice_ind}};
  
  % Read data files
  [DATA.FTDATA,msg]=l_ReadSliceData(slice_files,DATA.HDR,ShowWbar);
  if ~isempty(msg)
    DATA=[];
    return
  end
  
  % Include file information
  DATA.HDR.fname = fname;
  DATA.HDR.fpath = fpath;
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read *.HDR header file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [HDR,msg]=l_ReadHDRFile(filename)

HDR=[];
msg='';

% Try to open file for reading
fid = fopen(filename,'r');
if fid < 0
  msg = sprintf('Could not open file "%s"',filename);
  return
end

% Read parameter file
C=textscan(fid,'%s','delimiter','\n');
fclose(fid); % Close file
C=C{1};
if ~strncmp(C{1},'!INTERFILE',10)
  msg = sprintf('File "%s" is not valid SPECT/CT HDR-file',filename);
  return
end

try
  for ii=1:length(C)
    if strncmpi(C{ii},'!name of data file',18)
      if strcmpi(C{ii}(end),'=')
        HDR.DataFileName = '';
      else
        ind=findstr(C{ii},'=');
        HDR.DataFileName = C{ii}(ind+2:end);
      end
    elseif strncmpi(C{ii},'!patient ID',11)
      if strcmpi(C{ii}(end),'=')
        HDR.PatientID = '';
      else
        ind=findstr(C{ii},'=');
        HDR.PatientID = C{ii}(ind+2:end);
      end
    elseif strncmpi(C{ii},'!study ID',9)
      if strcmpi(C{ii}(end),'=')
        HDR.StudyID = '';
      else
        ind=findstr(C{ii},'=');
        HDR.StudyID = C{ii}(ind+2:end);
      end
    elseif strncmpi(C{ii},'!matrix size',12)
      if strncmpi(C{ii},'!matrix size [1]',16)
        ind=findstr(C{ii},'=');
        HDR.MatrixSize_1 = str2num(C{ii}(ind+2:end));
      elseif strncmpi(C{ii},'!matrix size [2]',16)
        ind=findstr(C{ii},'=');
        HDR.MatrixSize_2 = str2num(C{ii}(ind+2:end));
      elseif strncmpi(C{ii},'!matrix size [3]',16)
        ind=findstr(C{ii},'=');
        HDR.MatrixSize_3 = str2num(C{ii}(ind+2:end));
      end
    elseif strncmpi(C{ii},'!number format',14)
      if strcmpi(C{ii}(end),'=')
        HDR.NumberFormat = '';
      else
        ind=findstr(C{ii},'=');
        HDR.NumberFormat = C{ii}(ind+2:end);
      end
    elseif strncmpi(C{ii},'!number of bytes per pixel',26)
      ind=findstr(C{ii},'=');
      HDR.BytesPerPixel = str2num(C{ii}(ind+2:end));
    elseif strncmpi(C{ii},'!number of images this frame group',34)
      ind=findstr(C{ii},'=');
      HDR.ImagesInFrameGroup = str2num(C{ii}(ind+2:end));
    elseif strncmpi(C{ii},'!number of frame groups',23)
      ind=findstr(C{ii},'=');
      HDR.NbrOfFrameGroups = str2num(C{ii}(ind+2:end));
      
      %% Parse frame groups
      done = false;
      count = ii+1;
      framecount = 0;
      while ~done
        if strncmpi(C{count},'!frame group number',19)
          framecount = framecount+1;
          ind=findstr(C{count},'=');
          HDR.FrameGroups(framecount).FrameGroupNbr = ...
              str2num(C{count}(ind+2:end));
        elseif strncmpi(C{count},'!matrix size [1]',16)
          ind=findstr(C{count},'=');
          HDR.FrameGroups(framecount).MatrixSize_1 = str2num(C{count}(ind+2:end));
        elseif strncmpi(C{count},'!matrix size [2]',16)
          ind=findstr(C{count},'=');
          HDR.FrameGroups(framecount).MatrixSize_2 = ...
              str2num(C{count}(ind+2:end));
        elseif strncmpi(C{count},'!number format',14)
          if strcmpi(C{count}(end),'=')
            HDR.FrameGroups(framecount).NumberFormat = '';
          else
            ind=findstr(C{count},'=');
            HDR.FrameGroups(framecount).NumberFormat = C{count}(ind+2:end);
          end
        elseif strncmpi(C{count},'!number of bytes per pixel',26)
          ind=findstr(C{count},'=');
          HDR.FrameGroups(framecount).BytesPerPixel = ...
              str2num(C{count}(ind+2:end));
        elseif strncmpi(C{count},'!number of images this frame group',34)
          ind=findstr(C{count},'=');
          HDR.FrameGroups(framecount).ImagesInFrameGroup = ...
              str2num(C{count}(ind+2:end));
        elseif strncmpi(C{count},'image duration',14)
          ind=findstr(C{count},'=');
          HDR.FrameGroups(framecount).ImageDuration = ...
              str2num(C{count}(ind+2:end));
        elseif strncmpi(C{count},'!maximum pixel count in group',29)
          ind=findstr(C{count},'=');
          HDR.FrameGroups(framecount).MaxPixelCount = ...
              str2num(C{count}(ind+2:end));
          
          %% Check if we can bail out from the while loop
          if HDR.NbrOfFrameGroups==framecount
            done=true;
          end
        elseif strncmpi(C{count},'!END OF INTERFILE',17)
          %% End of file -> break from while loop
          done=true;
        end
        count =count+1;
      end
      
      %% Break away from the FOR-loop
      break;
    else
      continue
    end
  end
catch
  msg = sprintf('Error occurred while parsing file "%s"',filename)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read IMG file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,msg]=l_ReadImgFile(filename,HDR);

msg='';
data = [];

% Try to open img file for reading
fid = fopen(filename,'r');
if fid < 0
  msg = sprintf('Could not open file "%s" for reading',filename);
  return
end

% Data size for Frame Group Data
if isfield(HDR,'FrameGroups')
  
  %% Check if the slices in frame groups are of different sizes
  if all([HDR.FrameGroups(:).MatrixSize_1]==...
         HDR.FrameGroups(1).MatrixSize_1) && ...
        all([HDR.FrameGroups(:).MatrixSize_2]==...
            HDR.FrameGroups(1).MatrixSize_2)
    equalSize = true;
    
    % Allocate variable 'data'
    data=[];
  else
    equalSize = false;
    
    % Allocate variable 'data'
    data={};
  end
  
  %% Read frame groups separately
  for ii=1:HDR.NbrOfFrameGroups
    
    xSize = HDR.FrameGroups(ii).MatrixSize_1;
    ySize = HDR.FrameGroups(ii).MatrixSize_2;
    zSize = HDR.FrameGroups(ii).ImagesInFrameGroup;
    
    % Get data type
    if strcmpi(HDR.FrameGroups(ii).NumberFormat,'unsigned integer')
      if HDR.FrameGroups(ii).BytesPerPixel==2
        DataStr = '*uint16';
      elseif HDR.FrameGroups(ii).BytesPerPixel==1
        DataStr = '*uint8';
      end
    elseif strcmpi(HDR.FrameGroups(ii).NumberFormat,'signed integer')
      if HDR.FrameGroups(ii).BytesPerPixel==2
        DataStr = '*int16';
      elseif HDR.FrameGroups(ii).BytesPerPixel==1
        DataStr = '*int8';
      end
    end
    
    % Read data
    try
      data_tmp = fread(fid,xSize*ySize*zSize,DataStr);
      data_tmp = reshape(data_tmp,[xSize ySize zSize]);
      
      % Orient data
      data_tmp = permute(data_tmp,[2 1 3]);
      data_tmp = data_tmp(:,size(data_tmp,2):-1:1,:);
    catch
      msg = sprintf('Error while reading data from "%s"',filename);
      
      % Close file
      fclose(fid);
      return
    end
    
    if equalSize
      if ii==1
        data(:,:,end:end+(size(data_tmp,3)-1))=data_tmp;
      else
        data(:,:,end+1:end+(size(data_tmp,3)))=data_tmp;
      end
    else  
      for kk=1:size(data_tmp,3)
        data{end+1}=data_tmp(:,:,kk);
      end
    end
  end
else
  
  % Get data size from header
  xSize = HDR.MatrixSize_1;
  ySize = HDR.MatrixSize_2;
  if isfield(HDR,'MatrixSize_3')
    zSize = HDR.MatrixSize_3;
  elseif isfield(HDR,'ImagesInFrameGroup')
    zSize = HDR.ImagesInFrameGroup;
    if isfield(HDR,'NbrOfFrameGroups')
      zSize = zSize*HDR.NbrOfFrameGroups;
    end
  else
    zSize = 1;
  end

  % Get data type
  if strcmpi(HDR.NumberFormat,'unsigned integer')
    if HDR.BytesPerPixel==2
      DataStr = '*uint16';
    elseif HDR.BytesPerPixel==1
      DataStr = '*uint8';
  end
  elseif strcmpi(HDR.NumberFormat,'signed integer')
    if HDR.BytesPerPixel==2
      DataStr = '*int16';
    elseif HDR.BytesPerPixel==1
      DataStr = '*int8';
    end
  end
  
  % Read data
  try
    data = fread(fid,xSize*ySize*zSize,DataStr);
    data = reshape(data,[xSize ySize zSize]);
    
    % Orient data
    data = permute(data,[2 1 3]);
    data = data(:,size(data,2):-1:1,:);
  catch
    msg = sprintf('Error while reading data from "%s"',filename);
    
    % Close file
    fclose(fid);
    return
  end

end

% Close file
fclose(fid);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read XXM reconstruction parameter file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [HDR,msg]=l_ReadXXMFile(filename)


HDR=[];
msg='';

% Try to open file for reading
fid = fopen(filename,'r');
if fid < 0
  msg = sprintf('Could not open file "%s"',filename);
  return
end

% Read parameter file
C=textscan(fid,'%s','delimiter','\n');
fclose(fid); % Close file
C=C{1};

HDR.Tags={};

% Loop over lines in the file
for ii=1:length(C)
  
  if isempty(C{ii})
    continue
  end
  
  if strncmpi(C{ii},'BPMODETAG',9)
    HDR.Tags{end+1} = C{ii};
    continue;
  end
  
  if strncmpi(C{ii},'****',4)
    if strncmpi(C{ii},'**** R',6)
      HDR.Date = C{ii};
    elseif strncmpi(C{ii},'**** T',6)
      HDR.TotalReconstructionTime = C{ii};
    end
    continue
  end
  
  ind=findstr(C{ii},'=');
  if isempty(ind)
    continue
  end
  
  % Extract parameter and value
  param=C{ii}(1:ind-1);
  value=C{ii}(ind+1:end);
  
  % Test if the parameter is a number
  val_num=str2num(value);
  if isempty(val_num);
    val_num=value;
  end
  
  HDR.(param)=val_num;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Slice Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,msg]=l_ReadSliceData(slice_files,HDR,ShowWbar)


msg='';
% Get data size from header
xSize = HDR.PARTAG_CUBESIZEX;
ySize = HDR.PARTAG_CUBESIZEY;
nFiles = HDR.PARTAG_CUBESIZEZ;

% Allocate space for data
data = zeros([xSize,ySize,nFiles],'int16');

% Show waitbar
if ShowWbar
  hWbar=aedes_wbar(0,['Reading file 1/',num2str(nFiles)]);
  tmp_h=findall(hWbar,'string',['Reading file 1/',num2str(nFiles)]);
  set(tmp_h,'interpreter','none')
end

% Loop over slice files
for ii=1:length(slice_files)
  
  % Update waitbar
  if ShowWbar
    aedes_wbar(ii/nFiles,hWbar,...
			{['Reading file ',num2str(ii),'/',num2str(nFiles)],...
			slice_files{ii}});
  end
  
  fid = fopen(slice_files{ii},'r');
  if fid < 0
    data = [];
    msg = sprintf('Could not open file "%s" for reading!', ...
			slice_files{ii});
		if ShowWbar
			close(hWbar)
		end
    return
  end
  
  % Read data
  data(:,:,ii)=fread(fid,[xSize ySize],'*int16');
  fclose(fid); % Close file
end

% Permute data to correct orientation
data = permute(data,[2 1 3]);
data = data(:,size(data,2):-1:1,:);

if ShowWbar
  close(hWbar);
end
