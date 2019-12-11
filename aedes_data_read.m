function DATA = aedes_data_read(filename,varargin)
% AEDES_DATA_READ - Read various image data formats to data structure
%   
%
% Synopsis: 
%       DATA=aedes_data_read(filename,file_format,default_dir,varargin)
%       
%       or
%
%       DATA=aedes_data_read;  % (interactive mode, opens a file dialog)
%
% Description:
%       The function reads image data from various different file formats
%       into a DATA-structure. The first input argument "filename" is the
%       full path to the data file. If the first input argument is given as
%       an empty string, the open file dialog is shown. The second input
%       argument is format string that defines the data format; valid
%       format strings are: 'vnmr', 'nifti', 'sur', 'mri', 'dcm',
%       'spect/ct', 'mat'. If the format string is not given as an input
%       argument the file extension is used to determine the data format.
%
%       The "default_dir" input argument is a path string defining the
%       default directory for the open file dialog. If "default_dir" is
%       omitted, the current directory (pwd) is used to open the file
%       dialog.
%
% Examples:
%       DATA=aedes_data_read;   % Read image data
%       aedes(DATA)    % Open data in Aedes
%
% See also:
%       AEDES_READFID, AEDES_READCTDATA, AEDES_READ_NIFTI, AEDES

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


showWbar = true; % Show wbar by default
ddir = [pwd,filesep];
if ~isempty(varargin)
  if any(strcmpi(varargin,'default_dir'))
    ind=find((strcmpi(varargin,'default_dir')));
    try
      ddir = varargin{ind+1};
    catch
      ddir = [pwd,filesep];
    end
  end
end


%% Parse input arguments
if nargin<1 || isempty(filename)
  
  % Check if default directory is given
  %if nargin==3
  %  ddir = default_dir;
  %else
  %  ddir = [pwd,filesep];
  %end
  [filefilt,dataformats] = aedes_getfilefilter;
  [f_name, f_path, f_index] = uigetfile(...
	filefilt, ...
	'Select data file',ddir,...
	'MultiSelect', 'off');
  if isequal(f_name,0) % Cancel is pressed
	DATA=[];
	return
	end
	
	filename = fullfile(f_path,f_name);
	
	% There is a bug in uigetfile in OSX version of Matlab and the
	% fid-directory may be returned instead of the fid-file.
	if ismac
		if length(f_name)>3 && strcmpi(f_name(end-3:end),'.fid')
			f_path = [filename,filesep];
			f_name = 'fid';
			filename = [f_path,f_name];
		end
	end
	dataformat = aedes_getdataformat(filename);
	
  
elseif nargin>=1
  if ischar(filename)
    [f_path,f_name,f_ext] = fileparts(filename);
    f_path=[f_path,filesep];
	f_name = [f_name,f_ext];
	dataformat = aedes_getdataformat(filename);
  else
    error('First input argument has to be of class STRING!')
  end
end

% Parse varargin
for ii=1:2:length(varargin)
  switch lower(varargin{ii})
   case 'wbar'
    if strcmpi(varargin{ii+1},'on')
      showWbar = true;
    else
      showWbar = false;
    end
    
    case 'dataformat'
	  if ~isempty(varargin{ii+1})
		dataformat = varargin{ii+1};
	  end
      
    case 'default_dir'
      % This is just a dummy case to prevent from accidently going to
      % the "otherwise" case...
      
    otherwise
      error('Unknown parameter "%s"',varargin{ii})
    
  end
end


%% Read data
switch dataformat
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Read Matlab MAT-File
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
 case 'mat'

  % Show aedes_calc_wait
  if showWbar
    [h,txh]=aedes_calc_wait('Loading data from Matlab MAT-File...');
    drawnow
  end
  
  %% Set data format string
  DATA.DataFormat = 'mat';
  
  %% Check variables in the MAT-file
  tmp=who('-file',filename);
  if isempty(tmp) || ~(any(strncmpi(tmp,'data',4)) || ...
      any(strncmpi(tmp,'images',6)))
    if showWbar
      delete(h)
    end
    DATA=[];
    error('The MAT-file doesn''t contain the required variable "Data" or "images"!')
  end
  
  % Use the "data" variable by default
  ind=find(strncmpi(tmp,'data',4));
  if isempty(ind)
    % If "data" variable is not found, use the "images" variable
    ind=find(strncmpi(tmp,'images',6));
  end
  dataFieldName = tmp{ind(1)};
  
  if showWbar
    set(txh,'string',sprintf('%s\n%s',...
	  'Loading data from Matlab MAT-File...',...
	  ['using variable "',dataFieldName,'"']));
    drawnow
  end
  
  % Load MAT-file
  try
    img=load(filename,'-mat');
  catch
    if showWbar
      delete(h)
	end
	DATA=[];
	error('Could not read MAT-file "%s"',filename)
  end
  
  % Check if data is structure or matrix
  data = img.(dataFieldName);
  if isstruct(data) || iscell(data)
    DATA = data;
    if isfield(img,'DataRotation') && ...
        isfield(img,'DataFlip') && iscell(DATA)
      DATA{1}.DataRotation = img.DataRotation;
      DATA{1}.DataFlip = img.DataFlip;
    end
    if isfield(img,'SliceClim') && iscell(DATA)
      DATA{1}.SliceClim = img.SliceClim;
    end
  elseif isnumeric(data) || islogical(data)
     DATA.FTDATA = img.(dataFieldName);
     DATA.HDR.fname = f_name;
     DATA.HDR.fpath = f_path;
     DATA.HDR.DataFormat = dataformat;
  else
    if showWbar
      delete(h)
    end
    clear img data;
	DATA=[];
	error('The variable "DATA" is invalid!')
  end
	
  if showWbar
	pause(0.3)
    delete(h)
  end
  
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Read S.M.I.S. SUR-files
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 case {'sur','mri'}
  
  %% Set data format string
  try
	DATA=aedes_smisread(filename);
  catch
	DATA=[];
	error('Could not read file "%s"!',filename)
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Read S.M.I.S. MRD-files
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'mrd'
	
	% To be written ...
	error('Reading of S.M.I.S. MRD-Files has not been implemented!')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Read Analyze 7.5 and NIfTI files
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 case {'nifti'}
  
  % Show aedes_calc_wait
  if showWbar
    [h,txh]=aedes_calc_wait('Reading data in NIfTI/Analyze75 format...');
  end
  
  % Read NIfTI and Analyze 7.5 header
  [fp,fn,fe] = fileparts(filename);
  if ~isempty(fe) && strcmpi(fe,'.img')
    [DATA,msg]=aedes_read_nifti(filename);
    if isempty(DATA)
      if showWbar
        delete(h)
      end
      error('Could not read data from file "%s"!',filename)
    end
    delete(h);
  else
    [DATA,msg]=aedes_read_nifti(filename,'header');
    if isempty(DATA)
      if showWbar
        delete(h)
      end
      error('Could not read header from file "%s"!',filename)
    end
    
    % Read NIfTI and Analyze 7.5 format data
    [DATA,msg]=aedes_read_nifti(DATA.HDR);
    if isempty(DATA)
      if showWbar
        delete(h)
	  end
	  error('Could not read data from file "%s"!',filename)
    end
    if showWbar
      delete(h);
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Read Varian VNMR files (FID)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 case 'vnmr'
  
  if showWbar
    showWbar = 'on';
  else
    showWbar = 'off';
  end
  
  %% Read parameters from procpar file
  [procpar,msg]=aedes_readprocpar([f_path,'procpar']);%,'wbar',showWbar);
  if isempty(procpar)
    DATA=[];
    error(msg);
    return
  end
  
  %% Fallback defaults for reading VNMR files
  ReadfidReturn = 1;
  ReadfidDCcorrection = 'off';
  ReadfidZeropadding = 'auto';
  ReadfidSorting = 'on';
  ReadfidFastRead = 'on';
  ReadfidPrecision = 'single';
  OrientImages = 'on';
  RemoveEPIphaseIm = 'off';
  VnmrUseOldReadFcn = true;
  
  %% Get defaults for Return
  if ispref('Aedes','ReadfidReturn')
    ReadfidReturn = getpref('Aedes','ReadfidReturn');
  end

  %% Get defaults for DC correction
  if ispref('Aedes','ReadfidDCcorrection')
    if getpref('Aedes','ReadfidDCcorrection')
      ReadfidDCcorrection = 'on';
    else
      ReadfidDCcorrection = 'off';
    end
  end

  %% Get defaults for Zeropadding
  if ispref('Aedes','ReadfidZeropadding')
    if getpref('Aedes','ReadfidZeropadding')==0
      ReadfidZeropadding = 'off';
    elseif getpref('Aedes','ReadfidZeropadding')==1
      ReadfidZeropadding = 'on';
    elseif getpref('Aedes','ReadfidZeropadding')==2
      ReadfidZeropadding = 'auto';
    else
      ReadfidZeropadding = 'auto';
    end
  end

  %% Get defaults for Sorting
  if ispref('Aedes','ReadfidSorting')
    if getpref('Aedes','ReadfidSorting')
      ReadfidSorting = 'on';
    else
      ReadfidSorting = 'off';
    end
  end
  
  %% Get defaults for FastRead
  if ispref('Aedes','ReadfidFastRead')
    if getpref('Aedes','ReadfidFastRead')
      ReadfidFastRead = 'on';
    else
      ReadfidFastRead = 'off';
    end
  end
  
  %% Get defaults for Precision
  if ispref('Aedes','ReadfidPrecision')
    if strcmpi(getpref('Aedes','ReadfidPrecision'),'single')
      ReadfidPrecision = 'single';
    else
      ReadfidPrecision = 'double';
    end
  end
  
%   %% Get defaults for Reorienting EPI data
%   if ispref('Aedes','ReadfidReorientEPI')
%     ReorientEPI = getpref('Aedes','ReadfidReorientEPI');
%   else
%     ReorientEPI = 'off';
%   end
  
  %% Get defaults for Reorienting EPI data
  if ispref('Aedes','ReadfidOrientImages')
    OrientImages = getpref('Aedes','ReadfidOrientImages');
  else
    OrientImages = 'on';
  end
  
  %% Get defaults for removing phase image from EPI
  if ispref('Aedes','ReadfidRemoveEPIphaseIm')
    RemoveEPIphaseIm = getpref('Aedes','ReadfidRemoveEPIphaseIm');
  else
    RemoveEPIphaseIm = 'off';
  end
  
  %% Get default read function
  if ispref('Aedes','VnmrUseOldReadFcn')
    VnmrUseOldReadFcn = getpref('Aedes','VnmrUseOldReadFcn');
  else
    VnmrUseOldReadFcn = true;
  end
  
  
  
  %% Read data from fid file
  try
    if VnmrUseOldReadFcn
      DATA=aedes_readfid([f_path,'fid'],...
        'procpar',procpar,...
        'wbar',showWbar,...
        'Return',ReadfidReturn,...
        'DCcorrection',ReadfidDCcorrection,...
        'Zeropadding',ReadfidZeropadding,...
        'sorting',ReadfidSorting,...
        'FastRead',ReadfidFastRead,...
        'Precision',ReadfidPrecision,...
        'OrientImages',OrientImages,...
        'RemoveEPIphaseIm',RemoveEPIphaseIm);
    else
      DATA=aedes_readvnmr([f_path,'fid'],...
        'procpar',procpar,...
        'wbar',showWbar,...
        'Return',ReadfidReturn,...
        'DCcorrection',ReadfidDCcorrection,...
        'Zeropadding',ReadfidZeropadding,...
        'sorting',ReadfidSorting,...
        'FastRead',ReadfidFastRead,...
        'Precision',ReadfidPrecision,...
        'OrientImages',OrientImages,...
        'RemoveEPIphaseIm',RemoveEPIphaseIm);
    end
    if isempty(DATA)
      DATA=[];
      error('Unknown error while reading "%s".',[f_path,'fid'])
      return
    end
  catch
    DATA=[];
    error(lasterr)
	end
  
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Read Bruker file formats
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	case {'bruker_raw','bruker_reco'}

		% Read Bruker data file
		DATA = aedes_readbruker(filename);
		
		
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Read DICOM image files (DCM)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 case 'dcm'
  
  % Show aedes_calc_wait
  if showWbar
    [h,txh]=aedes_calc_wait('Reading data from DICOM file...');
  end
   
  %% Set data format string
  DATA.DataFormat = 'dcm';
  
  %% Read DICOM header
  try
    hdr = dicominfo(filename);
    DATA.HDR.FileHeader=hdr;
    DATA.HDR.fname = f_name;
    DATA.HDR.fpath = f_path;
  catch
	DATA=[];
	if showWbar
	  delete(h)
	end
	error('Could not read header from DICOM file "%s"',...
	  filename);
	
  end
  
  %% Read DICOM data
  try
    DATA.FTDATA=dicomread(hdr);
  catch
	DATA=[];
	if showWbar
	  delete(h)
	end
	error('Could not read image data from DICOM file "%s"',...
	  filename);
	end
	
	% Scale data if requested in DICOM headers
	if ~isempty(DATA) && isfield(hdr,'RescaleSlope') && isfield(hdr,'RescaleIntercept')
		DATA.FTDATA = double(DATA.FTDATA);
		DATA.FTDATA = DATA.FTDATA*hdr.RescaleSlope+hdr.RescaleIntercept;
	end
	
	
  if showWbar
	delete(h)
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Read SPECT/CT Files
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 case 'spect/ct'
  
  %% Read Data
  try
    [DATA,msg]=aedes_readctdata(filename);
    if isempty(DATA)
      return
    end
  catch
	DATA=[];
	error('An error occurred while reading SPECT/CT data from "%s"!',...
	  filename)
  end
  
  %% Set data format string
  DATA.DataFormat = 'spect/ct';
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Read Varian FDF Files
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'fdf'
    
    %% Read data
    DATA=aedes_readfdf(filename);
    if isempty(DATA)
      DATA=[];
      error('An error occurred while reading FDF file "%s"!',...
        filename)
    end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Read SWIFT SGL-Files
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'swift_sgl'
    
    % Read data
    [DATA,msg]=aedes_readswiftsgl(filename,...
      'datainterpsizeprompt',true,...
      'wbar',true);
    if isempty(DATA)
      error(msg)
    end
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Read Aedes ROI-Files
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'roi'
	
	%% Try to load the ROI-file
    try
      roifile = load(filename,'-mat');
	catch
	  error('Error while reading ROI-file "%s".',...
		filename)
    end
    
    %% See if the ROI-file includes file information
    if ~isfield(roifile,'FileInfo') || ...
        any(cellfun(@isempty,roifile.FileInfo.DataFileName))
	  error('The ROI-file does not contain file information.')
    end
    
    filenames=roifile.FileInfo.DataFileName;
    pathnames=roifile.FileInfo.DataPathName;
    
    DATA={};
	if showWbar
	  [calc_h,txh]=aedes_calc_wait({['Reading file 1/' num2str(length(filenames))], ...
		'""'});
	end

    %% Read the files
    for ii=1:length(filenames)
      
      %% Check if the file exists
      if exist([pathnames{ii},filenames{ii}],'file')==0
        delete(calc_h)
		DATA=[];
		error('Cannot find the data file "%s".',...
		  [pathnames{ii},filenames{ii}])
      end
      
      set(txh,'String',...
              sprintf('%s\n%s',...
                      ['Reading file ' num2str(ii) ...
                       '/' num2str(length(filenames))],...
                      ['"',pathnames{ii},filenames{ii},'"']))
      drawnow
	  
	  try
		DATA{ii}=aedes_data_read([pathnames{ii},filenames{ii}],'wbar','off');
	  catch
		error(['Unknown error occurred while reading ROI-file "%s". ',...
		  lasterr],...
		  [pathnames{ii},filenames{ii}])
	  end
    end
    delete(calc_h)
    
    %Dat.LoadRoiAtStartUp = true;
    
    %% Use rotation information
    if isfield(roifile,'RotateFlip')
	  if iscell(roifile.RotateFlip)
		RotateFlip3d = roifile.RotateFlip;
		DataRotation = false(1,length(DATA));
		DataFlip = false(1,length(DATA));
	  else
		DataRotation = roifile.RotateFlip.Rotate;
		DataFlip = roifile.RotateFlip.Flip;
		RotateFlip3d = {};
	  end
	else
	  DataRotation = false(1,length(DATA));
	  DataFlip = false(1,length(DATA));
	  RotateFlip3d = {};
	end
	
	%% Use SliceClim information
	if isfield(roifile,'SliceClim')
	  DATA{1}.SliceClim = SliceClim;
	end
    
    %% Assign the ROI-structure to the first structure
    DATA{1}.ROI = roifile.ROI;
    clear('roifile')
    
    %% Rotate images if necessary
	if ~isempty(RotateFlip3d)
	  
	  for ii=1:length(RotateFlip3d)
		if strcmpi(RotateFlip3d{ii}{1},'rotate')
		  k = RotateFlip3d{ii}{2};
		  dim = RotateFlip3d{ii}{3};
		  DATA{1}.FTDATA = aedes_rot3d(DATA{1}.FTDATA,k,dim);
		elseif strcmpi(RotateFlip3d{ii}{1},'flip')
		  dim = RotateFlip3d{ii}{2};
		  DATA{1}.FTDATA = flipdim(DATA{1}.FTDATA,dim);
		end
	  end
	  DATA{1}.RotateFlip3d = RotateFlip3d;
	elseif ( ~all(DataRotation==0) || ~all(DataFlip==0) ) 
      for ii=1:length(DATA)
        if DataRotation(ii)~=0
          DATA{ii}.FTDATA = rot90(DATA{ii}.FTDATA,DataRotation(ii));
        end
      
        if DataFlip(ii)~=0
          if DataFlip(ii)==1
            DATA{ii}.FTDATA = flipud(DATA{ii}.FTDATA);
          elseif DataFlip(ii)==2
            DATA{ii}.FTDATA = fliplr(DATA{ii}.FTDATA);
          end
        end
	  end
	  DATA{1}.DataRotation = DataRotation;
	  DATA{1}.DataFlip = DataFlip;
	end
	
    
  
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Unknown File Format
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  otherwise
	DATA=[];
	error('Unknown file format')
end

