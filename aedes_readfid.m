function [DATA,msg_out] = aedes_readfid(filename,varargin)
% AEDES_READFID - Read VNMR (Varian) FID-files 
%   
%
% Synopsis: 
%        [DATA,msg]=aedes_readfid(filename,'PropertyName1',value1,'PropertyName2',value2,...)
%        [DATA,msg]=aedes_readfid(filename,'header')
%        [DATA,msg]=aedes_readfid(filename,output_filename)
%
% Description:
%        The function reads VNMR (Varian) FID-file and return a structure
%        DATA with fields DataFormat, HDR, FTDATA, KSPACE, PROCPAR, and
%        PHASETABLE. The fields of the DATA structure are constructed as
%        follows: 
%        
%        DATA
%          |-> DataFormat         (string identifier for data format 'vnmr')
%          |-> HDR
%                |-> FileHeader   (data file header)
%                |-> BlockHeaders (data block headers, not stored by default)
%                |-> fname        (file name)
%                |-> fpath        (file path)
%                |-> Param        (parameter values used by AEDES_READFID to read the data)
%          |-> FTDATA             (real fourier transformed image data)
%          |-> KSPACE             (complex k-space, empty by default)
%          |-> PROCPAR            (parameters from procpar file)
%          |-> PHASETABLE         (phasetable)
%
%        The DATA structure is returned as empty (without HDR, FTDATA,
%        KSPACE, and PROCPAR fields) if an error has occurred during
%        reading. The error message is returned in the second output
%        argument msg. If AEDES_READFID is called with only one output argument
%        (i.e. without MSG), the error message is echoed to the workspace.
%        
%        The first input argument is either a string containing full path to
%        the FID-file or the header structure HDR. Only the data file header
%        can be read by giving a string 'header' as the second input
%        argument.
%
%        By default the k-space is not returned, i.e. the field KSPACE is
%        empty. The returned data can be adjusted by using the 'return'
%        property and values 1, 2, or 3 (see below for more information).
%
%        The supported property-value pairs in AEDES_READFID (property strings
%        are not case sensitive):
%
%        Property:        Value:                Description:
%        *********        ******                ************
%        'Return'      :  [ {1} | 2 | 3 | 4 ]   % 1=return only ftdata (default)
%                                               % 2=return only k-space
%                                               % 3=return both ftdata & kspace
%                                               % 4=return raw kspace
%
%        'FastRead'    :  [{'off'} | 'on' ]     % Enables an experimental 
%                                               % method for very fast reading 
%                                               % of fid-files. This is not as
%                                               % robust as the older
%                                               % method and can consume a lot
%                                               % of memory.
%
%        'FileChunkSize' : [ megabytes ]        % Chunk size in megabytes for 
%                                               % reading fid-files
%                                               % (default=500 MB).
%                                               % Lowering the chunk size 
%                                               % helps to conserve memory 
%                                               % when reading large files, 
%                                               % but is slower. This
%                                               % property has no effect if
%                                               % FastRead='off'.
%                                               
%
%        'OutputFile'  :  filename              % Writes the slices into
%                                               % NIfTI-format files
%                                               % using the given filename.
%
%        'DCcorrection':  [ {'off'} | 'on' ]    % 'on'=use DC correction
%                                               % 'off'=don't use DC correction
%                                               % (default)
%
%        'Procpar'     :  (procpar-structure)   % Input procpar
%                                               % structure. If this
%                                               % property is not set the
%                                               % procpar structure 
%                                               % will be read using
%                                               % AEDES_READPROCPAR
%
%        'ZeroPadding' :  ['off'|'on'|{'auto'}] % 'off' = force
%                                               % zeropadding off 
%                                               % 'on' = force
%                                               % zeropadding on (force
%                                               % images to be square)
%                                               % 'auto' = autoselect
%                                               % using relative FOV
%                                               % dimensions (default) 
%
%        'PhaseTable'  : (custom phase table)   % User-specified
%                                               % line order for
%                                               % k-space. The phase table
%                                               % is obtained from the file
%                                               % specified in
%                                               % procpar.petable if
%                                               % necessary.
%
%        'sorting'      : [ 'off' | {'on'} ]    % 'off' =disable k-space
%                                               % sorting, i.e. ignore the
%                                               % phase table and/or arrays.
%                                               % 'on' =sort k-space using
%                                               % phase table and/or array
%                                               % information if necessary
%                                               % (default)
%
%        'wbar'        : [ {'on'} | 'off' ]     % 'on'  = show waitbar
%                                               % 'off' = don't show waitbar
%
%        'FlipKspace'  : [ {'off'} | 'LR' | 
%                             'UD' | 'LRUD' ]   % 'off' = don't flip (default)
%                                               % 'LR' = left/right
%                                               % 'UD' = up/down
%                                               % 'LRUD' = left/right and
%                                               % up/down
%
%        'FlipInd'     : [ {'all'} | 'alt' | 
%                          (custom vector)  ]   % 'all' = flip all slices
%                                               % (default)
%                                               % 'alt' = flip every
%                                               % second slice
%                                               % (custom vector) = A
%                                               % custom vector containing
%                                               % indices to the flipped slices
%
%        'Precision'   : ['single'|{'double'}]  % Precision of the
%                                               % outputted data.
%                                               % 'single' = 32-bit float
%                                               % 'double' = 64-bit float
%
%        'OrientImages': [ {'on'} | 'off' ]     % Orient FTDATA after
%                                               % reading the data using
%                                               % PROCPAR.orient property.
%
%        'RemoveEPIphaseIm' : [{'on'}|'off']    % Remove the phase image
%                                               % from EPI data. This
%                                               % option only affects EPI
%                                               % images.
%
%        'EPIBlockSize' : [integer value]       % Block size (number of
%                                               % volumes) used for
%                                               % processing multireceiver
%                                               % EPI data. Default=100
%
%        'EPIPhasedArrayData' : ['on'|{'off'}]  % Return data from
%                                               % individual receivers from 
%                                               % phased array EPI files. 
%
% Examples:
%        [DATA,msg]=aedes_readfid(filename)             % Read image data from 'filename'
%        [DATA,msg]=aedes_readfid(filename,'header')    % Read only data file header
%
%        [DATA,msg]=aedes_readfid(filename,'return',1)  % Return only image data (default)
%        [DATA,msg]=aedes_readfid(filename,'return',2)  % Return only k-space
%        [DATA,msg]=aedes_readfid(filename,'return',3)  % Return both image data and k-space
%        
%        % Read first data file header and then image data and k-space
%        [DATA,msg]=aedes_readfid(filename,'header')
%        [DATA,msg]=aedes_readfid(DATA.HDR,'return',3)
%
%        % Read VNMR-data with default options and save slices in NIfTI
%        % format
%        DATA=aedes_readfid('example.fid','saved_slices.nii');
%
%        % If you want to use non-default options and also write
%        % NIfTI-files, use the property/value formalism, for example
%        DATA=aedes_readfid('example.fid','ZeroPadding','off',...
%                     'OutputFile','saved_slices.nii');
%
% See also:
%        AEDES_READFIDPREFS, AEDES_READPROCPAR, AEDES_READPHASETABLE, 
%        AEDES_DATA_READ, AEDES_WRITE_NIFTI, AEDES

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




if nargout==2
  Dat.ShowErrors = false;
  msg_out='';
else
  Dat.ShowErrors = true;
end

%% ---------------------
% Defaults
Dat.ReturnKSpace  = false;
Dat.ReturnFTData  = true;
Dat.DCcorrection  = false;
Dat.ZeroPadding = 2;
Dat.Sorting = true;
Dat.UsePhaseTable = true;
Dat.FastDataRead = true;
Dat.precision = 'single';
Dat.FileChunkSize = 500; % Chunk size (in MB) for FastRead

%% Other defaults
Dat.ShowWaitbar   = true;
procpar=[];
Dat.phasetable=[];
Dat.FlipKspace = 0;
Dat.FlipInd = 'all';
Dat.OutputFile = false;
Dat.ReturnRawKspace = false;
Dat.ReorientEPI = false;
Dat.RemoveEPIphaseIm = true;
Dat.EPIBlockSize = 100;
Dat.EPIPhasedArrayData = false;
Dat.OrientImages = true;

%% -------------------------------------------------------------------


%% Set data format label
DATA.DataFormat = 'vnmr';

%% Parse input arguments
if nargin==0 || isempty(filename)
  
  %% Get default directory
  try
    defaultdir = getpref('Aedes','GetDataFileDir');
  catch
    defaultdir = [pwd,filesep];
  end
  
  %% If no input arguments are given, ask for a file
  [f_name, f_path, f_index] = uigetfile( ...
       {'fid;FID','Varian FID-files (fid)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Select VNMR (Varian) FID file',defaultdir);
  if ( all(f_name==0) || all(f_path==0) ) % Cancel is pressed
    DATA=[];
    msg_out='Canceled';
    return
  end
  ReadHdr = true;
  ReadData = true;
  filename=[f_path,f_name];
  
  %% Set default directory to preferences
  setpref('Aedes','GetDataFileDir',f_path)
  
end

if nargin==1
  if isstruct(filename)
    hdr = filename;
    filename = [hdr.fpath,hdr.fname];
    ReadHdr = false;
    ReadData = true;
% $$$     ReturnKSpace = false;
% $$$     ReturnFTData = true;
    
  elseif ischar(filename)
    ReadHdr = true;
    ReadData = true;
% $$$     ReturnKSpace = false;
% $$$     ReturnFTData = true;
  end
elseif nargin==2
  if strcmpi(varargin{1},'header')
    ReadHdr = true;
    ReadData = false;
% $$$     ReturnKSpace = false;
  elseif ischar(varargin{1})
    ReadHdr = true;
    ReadData = true;
    Dat.OutputFile = varargin{1};
  else
    if ~Dat.ShowErrors
      DATA=[];
      msg_out=sprintf('Invalid second input argument "%s".',varargin{1});
      return
    else
      error('Invalid second input argument "%s".',varargin{1})
    end
  end
else
  if isstruct(filename)
    hdr = filename;
    filename = [hdr.fpath,hdr.fname];
    ReadHdr = false;
    ReadData = true;
  elseif isempty(filename)
    [f_name, f_path, f_index] = uigetfile( ...
        {'fid;FID','Varian FID-files (fid)'; ...
         '*.*', 'All Files (*.*)'}, ...
        'Select VNMR (Varian) FID file');
    if ( all(f_name==0) || all(f_path==0) ) % Cancel is pressed
      DATA=[];
      msg_out='Canceled';
      return
    end
    ReadHdr = true;
    ReadData = true;
    filename=[f_path,f_name];
  else
    ReadHdr = true;
    ReadData = true;
  end
  
  for ii=1:2:length(varargin)
    switch lower(varargin{ii})
     case 'return'
      if length(varargin)<2
        DATA=[];
        if ~Dat.ShowErrors
          msg_out='"Return" value not specified!';
        else
          error('"Return" value not specified!')
        end
        return
      else
        if varargin{ii+1}==1
          Dat.ReturnKSpace = false;
          Dat.ReturnFTData = true;
        elseif varargin{ii+1}==2
          Dat.ReturnKSpace = true;
          Dat.ReturnFTData = false;
        elseif varargin{ii+1}==3
          Dat.ReturnKSpace = true;
          Dat.ReturnFTData = true;
		elseif varargin{ii+1}==4
		  Dat.ReturnRawKspace = true;
		  Dat.ReturnKSpace = true;
          Dat.ReturnFTData = false;
        end
      end
      
     case {'dc','dccorrection','dccorr'}
      if strcmpi(varargin{ii+1},'on')
        Dat.DCcorrection = true;
      else
        Dat.DCcorrection = false;
      end
      
     case 'procpar'
      procpar=varargin{ii+1};
      
      case 'zeropadding'
      if ischar(varargin{ii+1})
        if strcmpi(varargin{ii+1},'on')
          Dat.ZeroPadding = 1; % on
        elseif strcmpi(varargin{ii+1},'off')
          Dat.ZeroPadding = 0; % off
        else
          Dat.ZeroPadding = 2; % auto
        end
      else
        % Undocumented
        Dat.ZeroPadding = 3; % Custom
        Dat.CustomZeroPadding = varargin{ii+1};
      end
      
      case 'phasetable'
      Dat.phasetable = varargin{ii+1};
      
      case 'sorting'
	  if strcmpi(varargin{ii+1},'on')
        Dat.UsePhaseTable = true;
        Dat.Sorting = true;
      else
        Dat.UsePhaseTable = false;
        Dat.Sorting = false;
	  end
      
	 case 'fastread'
     if strcmpi(varargin{ii+1},'on')
       Dat.FastDataRead = true;
     else
       Dat.FastDataRead = false;
     end
     
      case 'filechunksize'
        Dat.FileChunkSize = round(varargin{ii+1});
        
      case 'wbar'
        if strcmpi(varargin{ii+1},'on')
          Dat.ShowWaitbar = 1;
        else
          Dat.ShowWaitbar = 0;
        end
	  
	 case 'precision'
	   if strcmpi(varargin{ii+1},'single')
		 Dat.precision = 'single';
	   end
		
	  case 'flipkspace'
        
        flipkspace = varargin{ii+1};
        if strcmpi(flipkspace,'off')
          Dat.FlipKspace=0;
        elseif strcmpi(flipkspace,'LR')
          Dat.FlipKspace=1;
        elseif strcmpi(flipkspace,'UD')
          Dat.FlipKspace=2;
        elseif strcmpi(flipkspace,'LRUD')
          Dat.FlipKspace=3;
        else
          DATA=[];
          if ~Dat.ShowErrors
            msg_out = 'Bad value for property FlipKspace!';
          else
            error('Bad value for property FlipKspace!')
          end
          return
        end
        
      case 'flipind'
        Dat.FlipInd = varargin{ii+1};
        
      case 'outputfile'
        Dat.OutputFile = varargin{ii+1};
        
      case 'reorientepi'
        if strcmpi(varargin{ii+1},'on')
          Dat.ReorientEPI = true;
        end
        
      case 'removeepiphaseim'
        if strcmpi(varargin{ii+1},'on')
          Dat.RemoveEPIphaseIm = true;
        end
      case 'epiblocksize'
        Dat.EPIBlockSize = round(varargin{ii+1});
        
      case 'epiphasedarraydata'
        if strcmpi(varargin{ii+1},'on')
          Dat.EPIPhasedArrayData = true;
        else
          Dat.EPIPhasedArrayData = false;
        end
      case 'orientimages'
        if strcmpi(varargin{ii+1},'off')
          Dat.OrientImages = false;
        end
     otherwise
      DATA=[];
      if ~Dat.ShowErrors
        msg_out = ['Unknown property "' varargin{ii} '"'];
      else
        error(['Unknown property "' varargin{ii} '"'])
      end
      return
    end
  end
end

% Parse filename
[fpath,fname,fext]=fileparts(filename);
if ~strcmpi(fname,'fid')
  if isempty(fname)
    fpath = [fpath,filesep];
  else
	if isempty(fpath)
	  fpath = [pwd,filesep,fname,fext,filesep];
	else
	  fpath = [fpath,filesep,fname,fext,filesep];
	end
  end
  fname = 'fid';
else
  fpath = [fpath,filesep];
end

%% Read procpar if not given as input argument
if isempty(procpar) || nargin~=2
  [procpar,proc_msg]=aedes_readprocpar([fpath,'procpar']);
  if isempty(procpar)
    DATA=[];
    if ~Dat.ShowErrors
      msg_out=proc_msg;
    else
      error(proc_msg)
    end
    return
  end
end

%% Read phasetable if necessary
if isfield(procpar,'petable') && isempty(Dat.phasetable) && ...
    ~isempty(procpar.petable{1}) && ~strcmpi(procpar.petable{1},'n') && ...
    Dat.Sorting
  % Look in preferences for tablib-directory
  try
    tabpath = getpref('Aedes','TabLibDir');
  catch
    % If the table path was not found in preferences, try to use aedes
    % directory
    tmp_path = which('aedes');
    if ~isempty(tmp_path)
      aedes_path=fileparts(tmp_path);
      tabpath = [aedes_path,filesep,'tablib',filesep];
    else
      % If all else fails, look in the current directory
      tabpath = [pwd,filesep,'tablib',filesep];
    end
  end
  [phasetable,msg] = aedes_readphasetable([tabpath,procpar.petable{1}]);
  
  % If petable cannot be found, check if it is in procpar...
  if isempty(phasetable) && isfield(procpar,'pe_table')
    phasetable = {'t1',procpar.pe_table(:)};
  elseif isempty(phasetable) && isfield(procpar,'pelist') && ...
      ~isempty(procpar.pelist) && isnumeric(procpar.pelist)
    phasetable = {'t1',reshape(procpar.pelist,procpar.etl,[]).'};
  end
  
  % If the sequence is a fast spin echo, try to construct the phasetable
  if isempty(phasetable) && isfield(procpar,'petable') && ...
      strncmpi(procpar.petable{1},'fse',3)
    err_msg = 'Could not find phasetable!';
    if ~Dat.ShowErrors
      msg_out=err_msg;
    else
      error(err_msg)
    end
    return
    %phasetable = {'t1',l_CreateFsemsPhaseTable(procpar)};
  end
  
  %% Abort and throw error if phasetable cannot be read and the FID-file
  % has not been sorted
  if isempty(phasetable) && ~isfield(procpar,'flash_converted')
    DATA=[];
    if ~Dat.ShowErrors
      msg_out={['Could not find the required phase table "' ...
                procpar.petable{1} '" in the following folders'],...
               ['"' tabpath '".']};
    else
      error('Could not find the required phase table "%s" in %s',...
            procpar.petable{1},tabpath)
    end
    return
  elseif ( isempty(phasetable) && isfield(procpar,'flash_converted') ) || ...
      ~Dat.UsePhaseTable
    %% If the FID-file has been sorted, the reading can continue but
    % throw a warning.
    fprintf(1,'Warning: aedes_readfid: Could not find phasetable "%s" in "%s"!\n',procpar.petable{1},tabpath)
    Dat.phasetable=[];
  else
    Dat.phasetable = phasetable{1,2};
  end
end

% Convert phasetable indices to MATLAB indices
if ~isempty(Dat.phasetable)
  Dat.phasetable=Dat.phasetable+(-min(min(Dat.phasetable))+1);
else
  Dat.UsePhaseTable = false;
end

%% Open FID-file
[file_fid,msg]=fopen([fpath,fname],'r','ieee-be');
if file_fid < 0
  DATA=[];
  if ~Dat.ShowErrors
    msg_out=['Could not open file "' filename '" for reading.'];
  else
    error('Could not open file "%s" for reading.',filename)
  end
  return
end

% Read only header
if ReadHdr && ~ReadData
  [hdr,msg_out]=l_ReadDataFileHeader(file_fid);
  if ~isempty(msg_out)
    DATA=[];
    fclose(file_fid);
    if Dat.ShowErrors
      error(msg_out)
    end
    return
  else
    DATA.HDR.FileHeader = hdr.FileHeader;
    DATA.FTDATA=[];
    DATA.KSPACE=[];
    DATA.PROCPAR=[];
    DATA.PHASETABLE=[];
  end
elseif ~ReadHdr && ReadData % Header structure given as input argument
                            % Read only data. 
  [hdr,data,kspace,msg_out]=l_ReadBlockData(file_fid,hdr,Dat,procpar);
  if ~isempty(msg_out)
    DATA=[];
    fclose(file_fid);
    if Dat.ShowErrors
      error(msg_out)
    end
    return
  else
    DATA.HDR.FileHeader=hdr.FileHeader;
    DATA.HDR.BlockHeaders = hdr.BlockHeaders;
    DATA.FTDATA=data;
    DATA.KSPACE=kspace;
    DATA.PROCPAR=procpar;
    DATA.PHASETABLE=Dat.phasetable;
  end
elseif ReadHdr && ReadData  % Read headers and data
  [hdr,msg_out]=l_ReadDataFileHeader(file_fid);
  [hdr,data,kspace,msg_out]=l_ReadBlockData(file_fid,hdr,Dat,procpar);
  if ~isempty(msg_out)
    DATA=[];
    fclose(file_fid);
    if Dat.ShowErrors
      error(msg_out)
    end
    return
  else
    DATA.HDR.FileHeader=hdr.FileHeader;
    DATA.HDR.BlockHeaders = hdr.BlockHeaders;
    DATA.FTDATA=data;
    DATA.KSPACE=kspace;
    DATA.PROCPAR=procpar;
    DATA.PHASETABLE=Dat.phasetable;
  end
end

% Set file name and path to the HDR structure
DATA.HDR.fname=fname;
DATA.HDR.fpath=fpath;

% Set parameter values
DATA.HDR.Param.ReturnKSpace = Dat.ReturnKSpace;
DATA.HDR.Param.ReturnFTData = Dat.ReturnFTData;
if Dat.ZeroPadding==0
  DATA.HDR.Param.ZeroPadding = 'off';
elseif Dat.ZeroPadding==1
  DATA.HDR.Param.ZeroPadding = 'on';
else
  DATA.HDR.Param.ZeroPadding = 'auto';
end
if ~Dat.DCcorrection
  DATA.HDR.Param.DCcorrection = 'off';
else
  DATA.HDR.Param.DCcorrection = 'on';
end
if Dat.Sorting==0
  DATA.HDR.Param.Sorting = 'off';
else
  DATA.HDR.Param.Sorting = 'on';
end
if Dat.FlipKspace==0
  DATA.HDR.Param.FlipKspace = 'off';
elseif Dat.FlipKspace==1
  DATA.HDR.Param.FlipKspace = 'LR';
elseif Dat.FlipKspace==2
  DATA.HDR.Param.FlipKspace = 'UD';
elseif Dat.FlipKspace==3
  DATA.HDR.Param.FlipKspace = 'LRUD';
end
DATA.HDR.Param.FlipInd = Dat.FlipInd;


% Close file
fclose(file_fid);

%% Write slices to NIfTI files
if ischar(Dat.OutputFile)
  if isempty(Dat.OutputFile)
    [fp,fn,fe]=fileparts(DATA.HDR.fpath(1:end-1));
    fprefix=fn;
    niipath = [pwd,filesep];
  else
    [fp,fn,fe]=fileparts(Dat.OutputFile);
    if strcmpi(fe,'.nii')
      fprefix=fn;
    else
      fprefix=[fn,fe];
    end
    if isempty(fp)
      niipath = [pwd,filesep];
    else
      niipath = [fp,filesep];
    end
  end
  
  % Create file names
  filenames={};
  for ii=1:size(DATA.FTDATA,3)
    filenames{ii}=sprintf('%s%03d%s',[fprefix,'_'],ii,'.nii');
  end
  nFiles = length(filenames);
    
  h=aedes_wbar(0,sprintf('Saving slice 1/%d in NIfTI format...',nFiles));
  for ii=1:nFiles
    aedes_wbar(ii/nFiles,h,sprintf('Saving slice %d/%d in NIfTI format...',ii, ...
                             nFiles))
    [done,msg]=aedes_write_nifti(DATA.FTDATA(:,:,ii),...
                           [niipath,filenames{ii}],'DataType','single',...
                           'FileType',2);
  end
  delete(h)
  
  if ~done
    warning('Error occurred while writing NIfTI-files. Could not write file(s)!')
  end
  
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Data File (Main) Header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hdr,msg_out]=l_ReadDataFileHeader(file_fid)

msg_out='';
%% Read Data File Header
try
  FH.nblocks = fread(file_fid,1,'long');   % Number of blocks in file
  FH.ntraces = fread(file_fid,1,'long');   % Number of traces per block
  FH.np = fread(file_fid,1,'long');        % Number of elements per trace
  FH.ebytes = fread(file_fid,1,'long');    % Number of bytes per element
  FH.tbytes = fread(file_fid,1,'long');    % Number of bytes per trace
  FH.bbytes = fread(file_fid,1,'long');    % Number of bytes per block
  FH.vers_id = fread(file_fid,1,'short');  % Software version, file_id status bits
  FH.status = fread(file_fid,1,'short');   % Status of whole file
  FH.nbheaders = fread(file_fid,1,'long'); % Number of block headers per block
  
  hdr.FileHeader = FH;
  
  %% Parse status bits
  hdr.FileHeader.status=[];
  tmp_str = {'S_DATA',...          % 0 = no data, 1 = data
             'S_SPEC',...          % 0 = FID, 1 = spectrum
             'S_32',...            % 0 = 16 bit, 1 = 32 bit integer
             'S_FLOAT',...         % 0 = integer, 1 = floating point
             'S_COMPLEX',...       % 0 = real, 1 = complex
             'S_HYPERCOMPLEX',...  % 1 = hypercomplex
             '',...                % Unused bit
             'S_ACQPAR',...        % 0 = not Acqpar, 1 = Acqpar
             'S_SECND',...         % 0 = first FT, 1 = second FT
             'S_TRANSF',...        % 0 = regular, 1 = transposed
             '',...                % Unused bit
             'S_NP',...            % 1 = np dimension is active
             'S_NF',...            % 1 = nf dimension is active
             'S_NI',...            % 1 = ni dimension is active
             'S_NI2',...           % 1 = ni2 dimension is active
             ''...                 % Unused bit 
            };
  status_bits = fliplr(double(dec2bin(FH.status,16))==49);
  for ii=1:length(tmp_str)
    if ~isempty(tmp_str{ii})
      hdr.FileHeader.status.(tmp_str{ii}) = status_bits(ii);
    end
  end
catch
  hdr=[];
  msg_out=['Error occurred while reading file header from "' filename '".'];
  return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Block Headers and Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hdr,data,kspace,msg_out]=l_ReadBlockData(file_fid,hdr,Dat,procpar)

hdr.BlockHeaders = [];
%hdr.HypercmplxBHeaders = [];
data=[];
kspace=[];
count = [];
msg_out='';

if isfield(procpar,'seqcon')
  procpar.seqcon=procpar.seqcon{1};
else
  procpar.seqcon='nnnnn';
end

%% Determine acquisition type from seqcon parameter
if all(procpar.seqcon=='n')          % 1D spectroscopy
  Dat.AcqType=1;
elseif all(procpar.seqcon(4:5)=='n') % 2D imaging
  Dat.AcqType=2;
elseif procpar.seqcon(5)=='n'        % 3D imaging
  Dat.AcqType=3;
else                         % 4D imaging
  Dat.AcqType=4;
end

%% Double-check that AcqType is not 1D spectral
if isfield(procpar,'nv') && procpar.nv==0
  Dat.AcqType=1;
end

%% Use some heuristical approach to "triple-check" that the data is not
%  1D spectral
if ~isempty(strfind(procpar.seqfil{1},'STEAM')) || ...
      ~isempty(strfind(procpar.seqfil{1},'steam')) || ...
      ~isempty(strfind(procpar.seqfil{1},'LASER')) || ...
      ~isempty(strfind(procpar.seqfil{1},'laser')) || ...
      ~isempty(strfind(procpar.seqfil{1},'PRESS')) || ...
      ~isempty(strfind(procpar.seqfil{1},'press')) || ...
      ~isempty(strfind(procpar.seqfil{1},'1PULSE')) || ...
      ~isempty(strfind(procpar.seqfil{1},'1pulse')) || ...
      ~isempty(strfind(procpar.seqfil{1},'2PULSE')) || ...
      ~isempty(strfind(procpar.seqfil{1},'2pulse'))
  Dat.AcqType=1;
end

UsePhaseTable=Dat.UsePhaseTable;


%% If the data has been converted with flashc, the seqcon parameter needs
% to be changed
if isfield(procpar,'flash_converted')
  if isfield(procpar,'ni') && procpar.ni>1
    procpar.seqcon(3)='s';
  elseif ((procpar.seqcon(2)=='c') && (procpar.seqcon(3)=='c'))
    procpar.seqcon(2)='s';
  end
  UsePhaseTable=false; % Do not try to order data if flashc has been run
end

%% Set number of transients
if ~isfield(procpar,'nt')
  nt=1;
else
  nt=procpar.nt;
end

%% Determine if the arrayed acquisition was used
if isfield(procpar,'array') && ~isempty(procpar.array{1})
  
  if length(procpar.array)==1 && ~iscell(procpar.array{1}) && ...
      strcmp(procpar.array{1},'pad') && all(procpar.pad==0)
    % Skip the array parsing if the array is a dummy "pad"-array...
    Dat.isAcqArrayed = false;
    Dat.ArrayLength = 1;
  else
    Dat.isAcqArrayed = true;
    Dat.ArrayLength = [];
    
    % Determine array length
    for ii=1:length(procpar.array)
      if iscell(procpar.array{ii})
        Dat.ArrayLength(ii) = length(procpar.(procpar.array{ii}{1}));
      else
        Dat.ArrayLength(ii) = length(procpar.(procpar.array{ii}));
      end
    end
    Dat.ArrayLength = prod(Dat.ArrayLength);
  end
else
  Dat.isAcqArrayed = false;
  Dat.ArrayLength = 1;
end

%% Determine if the data is EPI data
if ( isfield(procpar,'readres') && isfield(procpar,'phaseres') ) || ...
    ( isfield(procpar,'apptype') && strcmp(procpar.apptype{1},'imEPI') )
  Dat.isEPIdata = true;
else
  Dat.isEPIdata = false;
end

BlockHeadStatusLabels = {'S_DATA',...       % 0 = no data, 1 = data
                    'S_SPEC',...          % 0 = FID, 1 = spectrum
                    'S_32',...            % 0 = 16 bit, 1 = 32 bit integer
                    'S_FLOAT',...         % 0 = integer, 1 = floating point
                    'S_COMPLEX',...       % 0 = real, 1 = complex
                    'S_HYPERCOMPLEX',...  % 1 = hypercomplex
                    '',...                % Unused bit
                    'MORE_BLOCKS',...     % 0 = absent, 1 = present
                    'NP_CMPLX',...        % 0 = real, 1 = complex
                    'NF_CMPLX',...        % 0 = real, 1 = complex
                    'NI_CMPLX',...        % 0 = real, 1 = complex
                    'NI2_CMPLX',...       % 0 = real, 1 = complex
                    '',...                % Unused bit
                    '',...                % Unused bit
                    '',...                % Unuesd bit
                    ''...                 % Unused bit 
                   };

BlockHeadModeLabels = {'NP_PHMODE',...   % 1 = ph mode
                    'NP_AVMODE',...    % 1 = av mode
                    'NP_PWRMODE',...   % 1 = pwr mode
                    '',...             % Unused bit
                    'NF_PHMODE',...    % 1 = ph mode
                    'NF_AVMODE',...    % 1 = av mode
                    'NF_PWRMODE',...   % 1 = pwr mode
                    '',...             % Unused bit
                    'NI_PHMODE',...    % 1 = ph mode
                    'NI_AVMODE',...    % 1 = av mode
                    'NI_PWRMODE',...   % 1 = pwr mode
                    '',...             % Unused bit
                    'NI2_PHMODE',...   % 1 = ph mode
                    'NI2_AVMODE',...   % 1 = av mode
                    'NI2_PWRMODE',...  % 1 = pwr mode
                    ''...              % Unused bit 
                   };


% The nbheaders -field is not correct in some cases. Thus, this field
% cannot be trusted and the real nbheaders has to be calculated from the
% byte counts.
tmp_bytes=hdr.FileHeader.bbytes-hdr.FileHeader.tbytes*hdr.FileHeader.ntraces;
header_bytes=28;
if rem(tmp_bytes,header_bytes)~=0
  msg_out = 'Block header byte count does not match with file header';
  return
else
  nbheaders = tmp_bytes/28;
end
%nbheaders = hdr.FileHeader.nbheaders;

%% Allocate space for k-space
% kspace=zeros(hdr.FileHeader.np/2,...
%              hdr.FileHeader.ntraces,hdr.FileHeader.nblocks);

if ~Dat.FastDataRead
  if any(Dat.AcqType==[1 2])
    switch procpar.seqcon(2:3)
      case {'cc','sc'}
        kspace = complex(zeros(hdr.FileHeader.np/2,...
          hdr.FileHeader.ntraces,...
          hdr.FileHeader.nblocks,Dat.precision));
      otherwise
        kspace = complex(zeros(hdr.FileHeader.np/2,...
          hdr.FileHeader.nblocks,...
          hdr.FileHeader.ntraces,Dat.precision));
    end
  else
    kspace = complex(zeros(hdr.FileHeader.np/2,...
      hdr.FileHeader.ntraces,...
      hdr.FileHeader.nblocks,Dat.precision));
  end
else
  %kspace = [];
   kspace = complex(zeros(hdr.FileHeader.np/2*hdr.FileHeader.ntraces,...
     hdr.FileHeader.nblocks,Dat.precision));
end

%% - The older robust (but also slower) way for reading the data. 
%% When the blocksize is relatively small, this is also quite fast.
if ~Dat.FastDataRead

  % Initialize waitbar
  if Dat.ShowWaitbar
	wb_h = aedes_wbar(0/hdr.FileHeader.nblocks,...
	  {['Reading ',num2str(Dat.AcqType),'D VNMR data (seqcon: "' procpar.seqcon '")'],...
	  ['(Processing data block ' ...
	  num2str(0) '/' num2str(hdr.FileHeader.nblocks) ')']});
  end

  %% Read data and block headers
  for ii=1:hdr.FileHeader.nblocks
	%% Update waitbar
	if Dat.ShowWaitbar
	  aedes_wbar(ii/hdr.FileHeader.nblocks,...
		wb_h,...
		{['Reading ',num2str(Dat.AcqType),'D VNMR data (seqcon: "' procpar.seqcon '")'],...
		['(Processing data block ' ...
		num2str(ii) '/' num2str(hdr.FileHeader.nblocks) ')']})
	end

	%% Read block header and hypercomplex header
	for kk=1:nbheaders
	  %% Read block header
	  if kk==1
		hdr.BlockHeaders.scale = fread(file_fid,1,'short');  % Scaling factor
		tmp_status = fread(file_fid,1,'short'); % Status of data in block
		hdr.BlockHeaders.status = [];
		hdr.BlockHeaders.index = fread(file_fid,1,'short');  % Block index
		tmp_mode = fread(file_fid,1,'short');   % Mode of data in block
		hdr.BlockHeaders.mode = [];
		hdr.BlockHeaders.ctcount = fread(file_fid,1,'long'); % ct value for FID
		hdr.BlockHeaders.lpval = fread(file_fid,1,'float');  % f2 (2D-f1) left phase in phasefile
		hdr.BlockHeaders.rpval = fread(file_fid,1,'float');  % f2 (2D-f1) right phase in phasefile
		hdr.BlockHeaders.lvl = fread(file_fid,1,'float');    % level drift correction
		hdr.BlockHeaders.tlt = fread(file_fid,1,'float');    % tilt drift correction

		%% Parse status and mode bits
		status_bits = fliplr(double(dec2bin(tmp_status,16))==49);
		mode_bits = fliplr(double(dec2bin(tmp_mode,16))==49);
		for tt=1:length(BlockHeadStatusLabels)
		  if ~isempty(BlockHeadStatusLabels{tt})
			hdr.BlockHeaders.status.(BlockHeadStatusLabels{tt}) = status_bits(tt);
		  end
		  if ~isempty(BlockHeadModeLabels{tt})
			hdr.BlockHeaders.mode.(BlockHeadModeLabels{tt}) = mode_bits(tt);
		  end
		end


	  else %% Read hypercomplex header
		fread(file_fid,1,'short'); % Spare
		hdr.BlockHeaders.HCHstatus = fread(file_fid,1,'short');
		fread(file_fid,1,'short'); % Spare
		fread(file_fid,1,'short'); % Spare
		fread(file_fid,1,'long'); % Spare
		hdr.BlockHeaders.HCHlpval1 = fread(file_fid,1,'float');
		hdr.BlockHeaders.HCHrpval1 = fread(file_fid,1,'float');
		fread(file_fid,1,'float'); % Spare
		fread(file_fid,1,'float'); % Spare
	  end
	end
	
	%% Check block index to be sure about the data type
	if hdr.BlockHeaders.index~=ii
	  fprintf(1,'Warning: Index mismatch in "%s" block %d\n',fopen(file_fid),ii);
	  
	  % Use information from the file header instead of the BlockHeader if
	  % there is a mismatch in blockheader index...
	  useFileHeader = true;
	else
	  useFileHeader = false;
	end
	
	%% Determine data precision
	if useFileHeader
	  if hdr.FileHeader.status.S_FLOAT==1
		prec_str = ['single=>',Dat.precision];
	  elseif hdr.FileHeader.status.S_32==1 ...
		  && hdr.FileHeader.status.S_FLOAT==0
		prec_str = ['int32=>',Dat.precision];
	  elseif hdr.FileHeader.status.S_32==0 ...
		  && hdr.FileHeader.status.S_FLOAT==0
		prec_str = ['int16=>',Dat.precision];
	  end
	  
	else
	  if hdr.BlockHeaders.status.S_FLOAT==1
		prec_str = ['single=>',Dat.precision];
	  elseif hdr.BlockHeaders.status.S_32==1 ...
		  && hdr.BlockHeaders.status.S_FLOAT==0
		prec_str = ['int32=>',Dat.precision];
	  elseif hdr.BlockHeaders.status.S_32==0 ...
		  && hdr.BlockHeaders.status.S_FLOAT==0
		prec_str = ['int16=>',Dat.precision];
	  end
	end

	% Read k-space
	tmp=fread(file_fid,...
	  [hdr.FileHeader.np,hdr.FileHeader.ntraces],...
	  prec_str);


	% Store complex block and perform DC correction
	if ~Dat.DCcorrection || ( nt(1)>1 )
	  complex_block = complex(tmp(1:2:end,:),tmp(2:2:end,:));
	else
	  complex_block = complex(tmp(1:2:end,:)-hdr.BlockHeaders.lvl,...
		tmp(2:2:end,:)-hdr.BlockHeaders.tlt);
	end


	%% Store and order k-space values
	if any(Dat.AcqType==[1 2])
	  switch procpar.seqcon(2:3)
		case {'cc','sc'}
		  kspace(:,:,ii) = complex_block;
		otherwise
		  kspace(:,ii,:) = complex_block;
	  end
	else
	  kspace(:,:,ii) = complex_block;
	end

	% Do not save blockheaders by default. When reading data files with a lot of
	% blocks (e.g. over 1000) the processing of the DATA structure can be
	% slowed down considerably. If you for some reason want to save also the
	% block headers in the DATA structure comment out the code line below.
	hdr.BlockHeaders = [];
  end % for ii=1:hdr.

  if Dat.ShowWaitbar
	delete(wb_h)
  end
  
else
  %% -------------------------------------------------------------------
  %% Fast Method for reading data. This may fail with some datas and can 
  %% also require a relatively large amount of memory. This
  %% method should be used for EPI datas that contain a large number
  %% of block headers...
  
  % Check the size of the FID-file
  d=dir(fopen(file_fid));
  file_sz = d.bytes/1024/1024; % File size in MB
  if file_sz<Dat.FileChunkSize
    nBlocks = 1;
  else
    nBlocks = ceil(file_sz/Dat.FileChunkSize); % Read data in 500 MB blocks
  end
  
  % Initialize waitbar
  if Dat.ShowWaitbar
    if nBlocks==1
      wb_h = aedes_calc_wait(['Reading ',num2str(Dat.AcqType),...
        'D VNMR data (seqcon: "' procpar.seqcon '")']);
    else
      wb_h = aedes_wbar(1/nBlocks,{['Reading ',num2str(Dat.AcqType),...
        'D VNMR data (seqcon: "' procpar.seqcon '")'],...
        sprintf('Reading block 0/%d',nBlocks)});
    end
  end
  
  % The first block header is read and it is assumed that the values in 
  % the other block headers don't change.  
  hdr.BlockHeaders.scale = fread(file_fid,1,'short');  % Scaling factor
  tmp_status = fread(file_fid,1,'short'); % Status of data in block
  hdr.BlockHeaders.status = [];
  hdr.BlockHeaders.index = fread(file_fid,1,'short');  % Block index
  tmp_mode = fread(file_fid,1,'short');   % Mode of data in block
  hdr.BlockHeaders.mode = [];
  hdr.BlockHeaders.ctcount = fread(file_fid,1,'long'); % ct value for FID
  hdr.BlockHeaders.lpval = fread(file_fid,1,'float');  % f2 (2D-f1) left phase in phasefile
  hdr.BlockHeaders.rpval = fread(file_fid,1,'float');  % f2 (2D-f1) right phase in phasefile
  hdr.BlockHeaders.lvl = fread(file_fid,1,'float');    % level drift correction
  hdr.BlockHeaders.tlt = fread(file_fid,1,'float');    % tilt drift correction
  
  %% Parse status and mode bits
  status_bits = fliplr(double(dec2bin(tmp_status,16))==49);
  mode_bits = fliplr(double(dec2bin(tmp_mode,16))==49);
  for tt=1:length(BlockHeadStatusLabels)
    if ~isempty(BlockHeadStatusLabels{tt})
      hdr.BlockHeaders.status.(BlockHeadStatusLabels{tt}) = status_bits(tt);
    end
    if ~isempty(BlockHeadModeLabels{tt})
      hdr.BlockHeaders.mode.(BlockHeadModeLabels{tt}) = mode_bits(tt);
    end
  end
  
  %% Determine data precision
  if hdr.BlockHeaders.status.S_FLOAT==1
    prec_str = ['single=>',Dat.precision];
    prec = 4; % Precision in bytes
  elseif hdr.BlockHeaders.status.S_32==1 ...
      && hdr.BlockHeaders.status.S_FLOAT==0
    prec_str = ['int32=>',Dat.precision];
    prec = 4;
  elseif hdr.BlockHeaders.status.S_32==0 ...
      && hdr.BlockHeaders.status.S_FLOAT==0
    prec_str = ['int16=>',Dat.precision];
    prec = 2;
  end
  
  % Seek the file back to the beginning of the first block header
  fseek(file_fid,-28,0);
  
  % Determine the number of values that will result from block header(s)
  nVals = (nbheaders*28)/prec;
  
  nbh = floor(hdr.FileHeader.nblocks/nBlocks);
  szh = nVals+hdr.FileHeader.np*hdr.FileHeader.ntraces;
  for ii=1:nBlocks
    if nBlocks~=1
      aedes_wbar(ii/nBlocks,wb_h,{['Reading ',num2str(Dat.AcqType),...
        'D VNMR data (seqcon: "' procpar.seqcon '")'],...
        sprintf('Reading block %d/%d',ii,nBlocks)});
    end
    
    % Read the whole data including block headers etc...
    if ii==nBlocks
      tmp = fread(file_fid,inf,prec_str);
    else
      tmp = fread(file_fid,nbh*szh,prec_str);
    end
    tmp=reshape(tmp,nVals+hdr.FileHeader.np*hdr.FileHeader.ntraces,[]);
    tmp(1:nVals,:)=[];
    
    if ii==nBlocks
      inds = ((ii-1)*nbh+1):size(kspace,2);
    else
      inds = ((ii-1)*nbh+1):ii*nbh;
    end
    
    % Do DC-correction if necessary
    if ~Dat.DCcorrection || ( nt(1)>1 )
      kspace(:,inds)=complex(tmp(1:2:end,:,:),tmp(2:2:end,:,:));
    else
      kspace(:,inds)=complex(tmp(1:2:end,:,:)-hdr.BlockHeaders.lvl,...
        tmp(2:2:end,:,:)-hdr.BlockHeaders.tlt);
    end
  end
  clear tmp
  
  % Transform to 3D matrix
  kspace = reshape(kspace,hdr.FileHeader.np/2,...
    hdr.FileHeader.ntraces,hdr.FileHeader.nblocks);
  
  %% Store and order k-space values
  if any(Dat.AcqType==[1 2])
    switch procpar.seqcon(2:3)
      case {'cc','sc'}
        %kspace(:,:,ii) = complex_block;
      otherwise
        %kspace(:,ii,:) = complex_block;
        kspace = permute(kspace,[1 3 2]);
    end
  else
    %kspace(:,:,ii) = complex_block;
  end
  
  if Dat.ShowWaitbar
    delete(wb_h)
    pause(0.1)
  end
  
end

% Remove singleton dimensions from kspace
kspace = squeeze(kspace);

% Check if raw kspace should be returned
if Dat.ReturnRawKspace
  %% Delete waitbar
  if Dat.ShowWaitbar && ishandle(wb_h)
	delete(wb_h)
  end
  return
end

%% Support for RASER sequence ---------------------------
if isfield(procpar,'teType')
  
  if strcmpi(procpar.sing_sh,'y') && strcmpi(procpar.teType,'c')
	
    % Separate reference scans from the data
    if isfield(procpar,'refscan') && strcmp(procpar.refscan,'y')
      kspace=permute(reshape(kspace,procpar.np/2,procpar.ne/2,2,[]),[1 2 4 3]);
      noise = kspace(:,:,:,1);
      kspace = kspace(:,:,:,2);
      % DC-correction
      dc_level = squeeze(mean(mean(noise)));
      for ii=1:size(kspace,3)
        kspace(:,:,ii)=kspace(:,:,ii)-dc_level(ii);
      end
    end
    
    if strcmpi(procpar.readType,'s')
      sw = procpar.sw;
      dwell = cumsum(ones(1,size(kspace,1))*(1/sw));
    else
      error('This feature has not been implemented yet!!!')
    end
    dwell = dwell*1000;
    %kspace = permute(kspace,[2 1 3]);
    
    
    
    % Phase correction
    dims = size(kspace);
    nI = size(kspace,3);
    nv = size(kspace,4);
    ne = size(kspace,2);
    np = size(kspace,1);
    g = 42.58;
    G = procpar.gvox2;
    Ds = procpar.vox2/100;
    pos2 = procpar.pos2/10;
    fudge = 0;
    fudge2 = 1;
   
    tt_dwell=repmat(dwell(:),1,ne);
    tt_ne = repmat((1:ne),length(dwell),1);
    tt_np = repmat((1:np)',1,ne);
    
    i=sqrt(-1);
    phi = (-2.*pi.*g.*G.*Ds.*tt_dwell.*(tt_ne./ne - 1/2 - pos2/Ds))+fudge.*tt_np;
    phi = repmat(phi,[1,1,size(kspace,3)]);
    kspace = abs(kspace).*exp(i*(angle(kspace)+phi));
    
    % Flip every other column
    sz(1)=size(kspace,1);
    sz(2)=size(kspace,2);
    sz(3)=size(kspace,3);
    sz(4)=size(kspace,4);
    kspace=reshape(kspace,sz(1),prod(sz(2:4)));
    kspace(:,1:2:end) = flipud(kspace(:,1:2:end));
    kspace = reshape(kspace,sz);
    
    %kspace = reshape(permute(kspace,[2 1 3]),sz(2),[],1);
    %kspace(:,1:2:end) = flipud(kspace(:,1:2:end));
    %kspace=permute(reshape(kspace,sz(2),[],sz(3)),[2 1 3]);
    
  else
    data = [];
    error('RASER sequence of this type has not been implemeted yet!')
  end
  
  % Reshape into 4D matrix
  kspace = reshape(kspace,size(kspace,1),size(kspace,2),[],size(kspace,3));
  
  % Reorient if requested
  if Dat.OrientImages
    kspace = flipdim(kspace,2);
  end
  
  data = abs(fftshift(fft(fftshift(fft(kspace,[],3),3),[],1),1));
  
  % Delete kspace if not returned
  if ~Dat.ReturnKSpace
    kspace=[];
  end
  
  return
  
end

%% Support for fat/water measurements ----------------------
if isfield(procpar,'seqfil') && strcmpi(procpar.seqfil,'ge3d_csi2')
  
  % Split kspace into fat and water (in 4th dimesion)
  kspace=reshape(permute(reshape(kspace,256,2,[]),[1 3 4 2]),256,128,[],2);
  
  % Fourier transform data
  if any(Dat.ZeroPadding==[1 2])
    data_sz = [procpar.np/2,procpar.nf,procpar.nv2,2];
    data = zeros(data_sz,class(kspace));
    data(:,:,:,1) = abs(fftshift(fftn(kspace(:,:,:,1),data_sz(1:3))));
    data(:,:,:,2) = abs(fftshift(fftn(kspace(:,:,:,2),data_sz(1:3))));
  else
    data = zeros(size(kspace),class(kspace));
    data(:,:,:,1) = abs(fftshift(fftn(kspace(:,:,:,1))));
    data(:,:,:,2) = abs(fftshift(fftn(kspace(:,:,:,2))));
  end
  
  % Delete kspace if not returned
  if ~Dat.ReturnKSpace
    kspace=[];
  end
  
  return
end

% Check the number of receivers (PI stuff)
if isfield(procpar,'rcvrs') && ~isempty(procpar.rcvrs) && ...
    length(procpar.rcvrs{1})>1
  % Multiple receivers used
  if isfield(procpar,'apptype') && strcmp(procpar.apptype{1},'imEPI')
    % Store multiple receiver data in EPI measurements in 5th dimension
    % and calculate sum-of-squares image
    nRcv = length(find(procpar.rcvrs{1}=='y'));
    if ~isfield(procpar,'epiref_type') || strcmpi(procpar.epiref_type,'single')
      nRef = 1;
    elseif strcmpi(procpar.epiref_type,'triple')
      nRef = 3;
    end
    nVols = size(kspace,3)/nRcv-nRef;
    if Dat.EPIPhasedArrayData
      data = zeros(procpar.nv,procpar.np/2,procpar.ns,nVols+nRef,nRcv,'single');
    else
      data = zeros(procpar.nv,procpar.np/2,procpar.ns,nVols+nRef,'single');
    end
    kssz=size(kspace);
    blksz = Dat.EPIBlockSize; % Process EPI data in 100 volume blocks (default)
    nBlocks = ceil((size(kspace,3)/nRcv-nRef)/blksz);
    lnum = length(num2str(nBlocks));
    lnumstr = num2str(lnum);
    bsl = lnum*2+1;
    fprintf(1,'Processing data in blocks of %d volumes\n',blksz)
    fprintf(1,['Processing block...%0',lnumstr,'d/%0',lnumstr,'d'],1,nBlocks);
    for ii=1:nBlocks
      fprintf(1,repmat('\b',1,bsl));
      fprintf(1,['%0',lnumstr,'d/%0',lnumstr,'d'],ii,nBlocks);
      tmp_data = [];
      for kk=1:nRcv
        inds_ref = kk:nRcv:nRcv*nRef;
        inds_im = (nRcv*nRef+kk):nRcv:kssz(3);
        inds = cat(2,inds_ref,inds_im(((ii-1)*blksz+1):min(ii*blksz,nVols)));
        tmp_kspace = l_ReconstructKspace(kspace(:,:,inds),procpar,Dat);
        tmp_data(:,:,:,:,kk) = fftshift(fftshift(fft(fft(tmp_kspace,[],1),[],2),1),2);
      end
      if Dat.EPIPhasedArrayData
        data_block = abs(tmp_data);
      else
        data_block = sqrt(sum(tmp_data.*conj(tmp_data),5));
      end
      if ii==1
        data(:,:,:,1:size(tmp_data,4),:) = data_block;
      elseif ii==nBlocks
        data_block(:,:,:,1:nRef,:)=[];
        data(:,:,:,((ii-1)*blksz+1+nRef):(nVols+nRef),:) = data_block;
      else
        data_block(:,:,:,1:nRef,:)=[];
        data(:,:,:,((ii-1)*blksz+1+nRef):(ii*blksz+nRef),:) = data_block;
      end
    end
    fprintf(1,'\n')
    
    % Remove reference image if requested
    if Dat.isEPIdata && Dat.RemoveEPIphaseIm
      data(:,:,:,1:nRef,:)=[];
    end
    
    if Dat.OrientImages && ~isempty(procpar) && ...
        isfield(procpar,'orient') && any(Dat.AcqType==[2 3 4])
      orient = procpar.orient{1};
      if any(strcmpi(orient,{'xyz','trans90','cor90','sag90'}))
        data = flipdim(aedes_rot3d(data,1,3),2);
      elseif strcmpi(orient,'oblique')
        data = flipdim(flipdim(data,1),2);
      else
        data = flipdim(flipdim(data,1),2);
      end
    end
    
  else
    nRcv = length(find(procpar.rcvrs{1}=='y'));
    data = [];
    kspace2 = [];
    for ii=1:nRcv
      tmp_kspace = l_ReconstructKspace(kspace(:,:,ii),procpar,Dat);
      kspace2(:,:,:,:,ii)=tmp_kspace;
      if Dat.ReturnFTData
        tmp_data = l_CalculateFFT(tmp_kspace,procpar,Dat);
        data(:,:,:,:,ii) = tmp_data;
      end
    end
    kspace = kspace2;
    kspace2=[];
  end
else
  % Only one receiver
  kspace = l_ReconstructKspace(kspace,procpar,Dat);
  data = l_CalculateFFT(kspace,procpar,Dat);
end

% Delete kspace if not returned
if ~Dat.ReturnKSpace
  kspace=[];
end

%===============================================
% Reconstruct kspace for different sequences
%===============================================
function kspace=l_ReconstructKspace(kspace,procpar,Dat)


%% Flip images for certain measurements
if Dat.FlipKspace~=0
  if ischar(Dat.FlipInd)
    if strcmpi(Dat.FlipInd,'all')
      flipind = 1:size(kspace,3);
    elseif strcmpi(Dat.FlipInd,'alt')
      flipind = 2:2:size(kspace,3);
    end
  else
    flipind = Dat.FlipInd;
  end
  
  for ii=flipind
    if Dat.FlipKspace==1
      kspace(:,:,ii)=fliplr(kspace(:,:,ii));
    elseif Dat.FlipKspace==2
      kspace(:,:,ii)=flipud(kspace(:,:,ii));
    else
      kspace(:,:,ii)=flipud(fliplr(kspace(:,:,ii)));
    end
  end
end

% Handle arrayed and RARE sequences
if Dat.isAcqArrayed && Dat.AcqType~=1 && Dat.Sorting && ~Dat.isEPIdata
  %if ~strcmpi(procpar.seqcon(2:3),'cc') && ( isempty(Dat.phasetable) | ...
  %        isfield(procpar,'flash_converted') )
  if (( isempty(Dat.phasetable) && ~strcmpi(procpar.seqcon(2:3),'sc') ) && ...
      ~(strcmpi(procpar.seqcon(2:3),'cc') && Dat.ArrayLength==size(kspace,3))) || ...
      isfield(procpar,'flash_converted')
                                  
    % Sort uncompressed arrayed data
    ks_order = 1:procpar.nv*Dat.ArrayLength;
    tmp = 0:procpar.nv:(procpar.ne-1)*procpar.nv;
    tmp=repmat(tmp,procpar.nv*Dat.ArrayLength,1);
    ks_order = reshape(ks_order,Dat.ArrayLength,procpar.nv).';
    ks_order = ks_order(:);
    ks_order = repmat(ks_order,1,procpar.ne);
    ks_order = ks_order+tmp;
    ks_order = ks_order.';
    
    if procpar.ns==1
      kspace=kspace(:,ks_order);
    else
      kspace = kspace(:,ks_order,:);
    end
    
    % Reshape into 3D matrix
    kspace=reshape(kspace,[size(kspace,1) procpar.nv Dat.ArrayLength* ...
                        procpar.ns]);
  else
    %% Sort RARE type data
    if Dat.UsePhaseTable && ~isempty(Dat.phasetable)
      %Dat.phasetable = Dat.phasetable';
      if Dat.isAcqArrayed && Dat.ArrayLength>1 && strcmpi(procpar.seqcon(2:3),'cs')
        kspace = permute(reshape(reshape(kspace,size(kspace,1),[]),size(kspace,1),Dat.ArrayLength,[]),[1 3 2]);
      else
        Dat.phasetable = Dat.phasetable.';
      end
      kspace(:,Dat.phasetable(:),:) = kspace;
    end
  end
elseif Dat.isEPIdata
  %% Support for EPI-measurements
  
  
  % EPI data is measured with old VNMR
  if isfield(procpar,'readres') && isfield(procpar,'phaseres')
    
    % Number of slices
    tmp_ns=length(procpar.pss);
    
    if isfield(procpar,'navecho') && strcmpi(procpar.navecho{1},'y')
      tmp_nv = procpar.nv-procpar.nseg;
    else
      tmp_nv = procpar.nv;
    end
    kspace = reshape(kspace,[size(kspace,1) ...
      size(kspace,2)/tmp_nv tmp_nv]);
    kspace = permute(kspace,[1 3 2]);
    
    % Reshape to 4D matrix
    kspace = reshape(kspace,[size(kspace,1) size(kspace,2) ...
      tmp_ns size(kspace,3)/tmp_ns]);
    
    
  elseif isfield(procpar,'apptype') && strcmp(procpar.apptype{1},'imEPI')
    % If EPI data is measured with new VNMRj system
    
    % Number of slices
    ns = length(procpar.pss);
    
    % Reshape kspace to 4D matrix
    if isfield(procpar,'fract_ky') && procpar.fract_ky~=procpar.nv/2
      kspace = reshape(kspace,procpar.np/2,procpar.nv/2+procpar.fract_ky,...
        ns,[]);
      if Dat.ZeroPadding~=0
        kspace(procpar.np/2,procpar.nv,1)=0;
      end
    else
      kspace = reshape(kspace,procpar.np/2,procpar.nv,...
        ns,[]);
    end
    
    % Flip even kspace lines
    if ~Dat.FastDataRead
      for tt=2:2:size(kspace,2)
        kspace(:,tt,:,:) = flipdim(kspace(:,tt,:,:),1);
      end
    else
      kspace(:,2:2:end,:,:) = flipdim(kspace(:,2:2:end,:,:),1);
    end
    
    
    % Sort centric measurements
    if isfield(procpar,'ky_order') && strcmpi(procpar.ky_order,'c')
      kspace(:,Dat.phasetable(:),:,:)=kspace;
    end
    %kspace = flipdim(kspace,1);
    
    % EPI phase correction -------------------------
    if ~isfield(procpar,'epiref_type') || strcmpi(procpar.epiref_type,'single')
      % Single reference pointwise phase correction
      
      % Get the reference image
      ref_im = kspace(:,:,:,1);
      
      % Do phase correction.
      phase_e = exp(-sqrt(-1)*angle(fft(ref_im,[],1)));
      for kk=2:size(kspace,4)
        kspace(:,:,:,kk) = ifft(fft(kspace(:,:,:,kk),[],1).*phase_e,[],1);
      end
    elseif strcmpi(procpar.epiref_type,'triple')
      % Triple reference pointwise phase correction
      
      % Get the reference images
      ref1 = kspace(:,:,:,1);
      phase1 = exp(-sqrt(-1)*angle(fft(ref1,[],1)));
      ref2 = flipdim(kspace(:,:,:,3),1);
      phase2 = exp(-sqrt(-1)*angle(fft(ref2,[],1)));
      im1 = flipdim(kspace(:,:,:,2),1);
      
      % Correct phase for reversed read gradient image
      rev_phase = fft(im1,[],1).*phase2;
      for kk=4:size(kspace,4)
        kspace(:,:,:,kk)=ifft(rev_phase+fft(kspace(:,:,:,kk),[],1).*phase1);
      end
    end
  end
  
else
  % This section contains various "chewing gum and ironwire" type
  % patches that hopefully work...
  
  if strcmp(procpar.seqcon,'nncnn')
    kspace = permute(kspace,[1 3 2]);
  elseif strcmp(procpar.seqcon,'nccnn') && length(size(kspace))==2 && ...
      procpar.ns>=1 && Dat.AcqType~=1
    if ~isempty(Dat.phasetable)
      kssz = size(kspace);
      phsz = size(Dat.phasetable);
      kspace=permute(reshape(kspace,procpar.np/2,...
        phsz(2),...
        kssz(2)/phsz(2)),[1 3 2]);
      kspace = permute(reshape(kspace,...
        procpar.np/2,procpar.ns,kssz(2)/procpar.ns),[1 3 2]);
      if Dat.Sorting
        kspace(:,Dat.phasetable(:),:)=kspace;
      end
    else
      kspace = reshape(kspace,...
        [procpar.np/2 procpar.ns ...
        size(kspace,2)/procpar.ns]);
      kspace=permute(kspace,[1 3 2]);
    end
    
  elseif strcmp(procpar.seqcon,'nscsn') && length(size(kspace))==3 && ...
	  Dat.AcqType~=1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Support for 3D fast spin-echo
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if ~isempty(Dat.phasetable)
	  for ii=1:size(kspace,3)
		tmp = kspace(:,:,ii);
		kssz = size(tmp);
		phsz = size(Dat.phasetable);
		tmp=permute(reshape(tmp,procpar.np/2,...
		  phsz(2),...
		  kssz(2)/phsz(2)),[1 3 2]);
		tmp = permute(reshape(tmp,...
		  procpar.np/2,procpar.ns,kssz(2)/procpar.ns),[1 3 2]);
		if Dat.Sorting
		  tmp(:,Dat.phasetable(:),:)=tmp;
		end
		kspace(:,:,ii)=tmp;
	  end
	end
	
	
  elseif Dat.AcqType~=1 && isfield(procpar,'nv')
    if length(size(kspace))==2 && ...
              ( size(kspace,1)~=(procpar.np/2) || size(kspace,2)~=procpar.nv )
      %% Reshape the kspace to the appropriate size
      kspace=reshape(kspace,[procpar.np/2 procpar.nv size(kspace,2)/ ...
                          procpar.nv]);
    elseif length(size(kspace))==3
      kspace = reshape(kspace,procpar.np/2,procpar.nv,[]);
      if Dat.Sorting && ~isempty(Dat.phasetable)
        kspace(:,Dat.phasetable(:),:) = kspace; 
      end
    end
  end
  
  
  %%% Support for Teemu's ASE3D-data
  %if strcmpi(procpar.seqcon,'ncccn')
  %  %% Reshape the kspace to the appropriate size
  %  kspace=reshape(kspace,[procpar.np/2 procpar.nv size(kspace,2)/procpar.nv]);
  %end
end



%=========================================
% Fourier Transform Data
%=========================================
function data=l_CalculateFFT(kspace,procpar,Dat)
data=[];

% Return image/spectral data --------------------------------
if Dat.ReturnFTData
  %% Fourier transform spectral data
  if Dat.AcqType==1
	wb_h = aedes_wbar(0,'Fourier transforming spectral data');
    %data=zeros(size(kspace),'single');
	data=kspace;
    sz=size(kspace,2)*size(kspace,3);
    count=1;
    for kk=1:size(kspace,3)
      for ii=1:size(kspace,2)
        %% Update waitbar
        if Dat.ShowWaitbar
          aedes_wbar(count/sz,...
               wb_h,...
               {'Fourier transforming spectral data',...
                ['(Processing spectrum ' num2str(count) '/' num2str(sz) ')']})
		end
        data(:,ii,kk) = abs(fftshift(fft(data(:,ii,kk))));
        count=count+1;
      end
    end
  
  %% Fourier transform image data
  else
    
    %% Zeropadding auto
    if Dat.ZeroPadding==2
      
      if Dat.isEPIdata && isfield(procpar,'readres') && isfield(procpar,'phaseres')
        
        %% Determine the image size for EPI data for procpar "readres"
        %% and "phaseres" fields
        data_sz = [procpar.readres procpar.phaseres];
        if isequal(data_sz,[size(kspace,1),size(kspace,2)])
          DoZeroPadding=false;
        else
          DoZeroPadding=true;
        end
        data_sz = [procpar.readres ...
          procpar.phaseres size(kspace,3) size(kspace,4)];
      else
        
        %% Check if zeropadding is necessary.
        ks_sz_sorted = sort([size(kspace,1),size(kspace,2)]);
        fov_sz_sorted = sort([procpar.lro,procpar.lpe]);
        sliceRelativeDim = ks_sz_sorted(2)/ks_sz_sorted(1);
        FOVrelativeDim = fov_sz_sorted(2)/fov_sz_sorted(1);
        if FOVrelativeDim~=sliceRelativeDim
          ind=find([size(kspace,1),size(kspace,2)]==ks_sz_sorted(1));
          data_sz = size(kspace);
          data_sz(ind) = round((fov_sz_sorted(1)/fov_sz_sorted(2))*ks_sz_sorted(2));
          DoZeroPadding=true;
        else
          data_sz = size(kspace);
          DoZeroPadding=false;
        end
      end
      
      %% Force zeropadding on. Force images to be square.
    elseif Dat.ZeroPadding==1
      ks_sz_max = max(size(kspace,1),size(kspace,2));
      data_sz = size(kspace);
      data_sz(1:2)=ks_sz_max;
      DoZeroPadding=true;
      
      %% Force zeropadding off
    else
      data_sz = size(kspace);
      DoZeroPadding=false;
    end
    
    %% Fourier transform 2D and EPI image data
    if Dat.AcqType==2 || Dat.isEPIdata
      sz=size(kspace,3);
      
      if DoZeroPadding
        data=zeros(data_sz,class(kspace));
        data(1:size(kspace,1),1:size(kspace,2),:,:)=kspace;
      else
        data=kspace;
      end
      if ~Dat.ReturnKSpace
        kspace = [];
      end
      if Dat.ShowWaitbar
        [wb_h,txh]=aedes_calc_wait('Fourier transforming image data...');
      end
     
      %data=permute(fftshift(fftshift(abs(fft(permute(fft(data),[2 1 3 4]))),1),2),[2 1 3 4]);
      if ~Dat.FastDataRead
        for tt=1:size(data,4)
          data(:,:,:,tt)=fftshift(fftshift(abs(fft(fft(data(:,:,:,tt),[],1),[],2)),1),2);
        end
      else
        data=fftshift(fftshift(abs(fft(fft(data,[],1),[],2)),1),2);
      end
    else
      %% Fourier transform 3D image data
      if Dat.ShowWaitbar
        [wb_h,txh]=aedes_calc_wait('Fourier transforming 3D image data...');
      end
      
      if DoZeroPadding
        % Check zeropadding in the 3rd dimension if zeropadding is set to
        % "auto"
        if Dat.ZeroPadding==2
          if isfield(procpar,'lpe2') && isfield(procpar,'lpe')
            data_sz(3) = max(1,round((procpar.lpe2/procpar.lpe)*data_sz(2)*Dat.ArrayLength));
          end
        end
        data=zeros(data_sz,class(kspace));
        data = abs(fftshift(fftn(kspace,data_sz)));
      else
        data=zeros(data_sz,class(kspace));
        data = abs(fftshift(fftn(kspace)));
      end
    end
    
    % Remove reference image if requested
    if Dat.isEPIdata && Dat.RemoveEPIphaseIm
      data = data(:,:,:,2:end);
    end
    
    % Reorient images and remove the reference image if requested
    if Dat.OrientImages && ~isempty(procpar) && ...
        isfield(procpar,'orient') && any(Dat.AcqType==[2 3 4])
      orient = procpar.orient{1};
      if any(strcmpi(orient,{'xyz','trans90','cor90','sag90'}))
        data = flipdim(aedes_rot3d(data,1,3),2);
      elseif strcmpi(orient,'oblique')
        data = flipdim(flipdim(data,1),2);
      else
        data = flipdim(flipdim(data,1),2);
      end
    end
    
   
    
% $$$   %% Fourier transform 4D image data
% $$$   elseif Dat.AcqType==4
% $$$     data = abs(fftshift(fftshift(fft2(kspace),2),1));
  end
  
  %% Delete waitbar
  if Dat.ShowWaitbar && ishandle(wb_h)
    delete(wb_h)
  end
end





%========================================================
% A function for creating phasetables for fse and fsems
%========================================================
function phasetable=l_CreateFsemsPhaseTable(procpar)

phasetable = [];
if ~isfield(procpar,'etl') || ~isfield(procpar,'kzero')
  return
end
etl = procpar.etl;
kzero = procpar.kzero;
nv = procpar.nv;

t = (-nv/2+1):(nv/2);
t = t(:);
tt = flipud(t);
tt = tt(:);

phasetable = [reshape(t(1:nv/2),[],etl);...
  flipud(reshape(tt(1:nv/2),[],etl))];
phasetable = circshift(fliplr(phasetable),[0 kzero-1]);


%% - EOF - %%
return