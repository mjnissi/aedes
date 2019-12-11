function [DATA,msg_out] = aedes_readvnmr(filename,varargin)
% AEDES_READVNMR - Read VNMR (Varian) FID-files 
%   
%
% Synopsis: 
%        [DATA,msg]=aedes_readvnmr(filename,'PropertyName1',value1,'PropertyName2',value2,...)
%        [DATA,msg]=aedes_readvnmr(filename,'header')
%        [DATA,msg]=aedes_readvnmr(filename,output_filename)
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
%        argument msg. If AEDES_READVNMR is called with only one output argument
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
%        The supported property-value pairs in AEDES_READVNMR (property strings
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
%        'SumOfSquares' :  [ {1} | 2 | 3 ]      % 1=Return only the 
%                                               % sum-of-squares image
%                                               % for multireceiver
%                                               % data (default).
%                                               % 2=Return both SoQ and
%                                               % individual receiver data
%                                               % 3=Return only individual
%                                               % receiver data
%                                               % NOTE: This property has
%                                               % no effect on single
%                                               % receiver data.
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
%        'SortPSS'     :  [ 'off' | {'on'} ]    % Sort slices in 2D stacks 
%                                               % using pss. Turn this 'off' 
%                                               % if interleaved slice
%                                               % order is to be preserved.
%                                               % The original pss is saved
%                                               % in procpar.pss_orig
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
%        'Precision'   : [{'single'}|'double']  % Precision of the
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
%                                               % EPI data. Default=300
%
%        'EPIPhasedArrayData' : ['on'|{'off'}]  % Return data from
%                                               % individual receivers from 
%                                               % phased array EPI files.
%
%        'EPI_PC' : [{'auto'}|'pointwise'|      % Phase correction for EPI. 
%                     'triple'|'off']           % 'auto' = Choose correction
%                                               % based on procpar.
%                                               % 'pointwise' = Pointwise 
%                                               % single reference.
%                                               % 'triple' = Use triple
%                                               % reference correction
%                                               % 'off' = Do not perform
%                                               % phase correction.
%                                                         
%
% Examples:
%        [DATA,msg]=aedes_readvnmr(filename)             % Read image data from 'filename'
%        [DATA,msg]=aedes_readvnmr(filename,'header')    % Read only data file header
%
%        [DATA,msg]=aedes_readvnmr(filename,'return',1)  % Return only image data (default)
%        [DATA,msg]=aedes_readvnmr(filename,'return',2)  % Return only k-space
%        [DATA,msg]=aedes_readvnmr(filename,'return',3)  % Return both image data and k-space
%        
%        % Read first data file header and then image data and k-space
%        [DATA,msg]=aedes_readvnmr(filename,'header')
%        [DATA,msg]=aedes_readvnmr(DATA.HDR,'return',3)
%
%        % Read VNMR-data with default options and save slices in NIfTI
%        % format
%        DATA=aedes_readvnmr('example.fid','saved_slices.nii');
%
%        % If you want to use non-default options and also write
%        % NIfTI-files, use the property/value formalism, for example
%        DATA=aedes_readvnmr('example.fid','ZeroPadding','off',...
%                     'OutputFile','saved_slices.nii');
%
% See also:
%        AEDES_READFIDPREFS, AEDES_READPROCPAR, AEDES_READPHASETABLE, 
%        AEDES_DATA_READ, AEDES_WRITE_NIFTI, AEDES

% This function is a part of Aedes - A graphical tool for analyzing 
% medical images
%
% Copyright (C) 2010 Juha-Pekka Niskanen <Juha-Pekka.Niskanen@uef.fi>
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
Dat.ReturnRawKspace = false;
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
Dat.OutputFile = false;
Dat.ReorientEPI = false;
Dat.RemoveEPIphaseIm = true;
Dat.EPIBlockSize = 300;
Dat.EPIPhasedArrayData = false;
Dat.EPI_PC = 'auto';
Dat.OrientImages = true;
Dat.SumOfSquares = 1;
Dat.UseCustomRecon = true;
Dat.SortPSS = true;

% Default custom reconstruction file directory
tmp = which('aedes');
if ~isempty(tmp)
  tmp_path=fileparts(tmp);
  recon_dir = [tmp_path,filesep,'vnmr_recon',filesep];
else
  % If all else fails, look in the current directory
  recon_dir = [pwd,filesep,'vnmr_recon',filesep];
end

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
  elseif ischar(filename)
    ReadHdr = true;
    ReadData = true;
  end
elseif nargin==2
  if strcmpi(varargin{1},'header')
    ReadHdr = true;
    ReadData = false;
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
        
      case {'sumofsquares','soq'}
        Dat.SumOfSquares = varargin{ii+1};
         
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
        elseif strcmpi(varargin{ii+1},'double')
          Dat.precision = 'double';
        else
          warning('Unsupported precision "%s". Using single...',varargin{ii+1});
          Dat.precision = 'single';
        end
        
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
			case 'epi_pc'
				
				if any(strcmpi(varargin{ii+1},...
						{'off','auto','triple','pointwise'}))
					Dat.EPI_PC = varargin{ii+1};
				else
					if ~Dat.ShowErrors
						msg_out=sprintf('Unknown EPI phase correction type "%s".',varargin{ii+1});
					else
						error('Unknown EPI phase correction type "%s".',varargin{ii+1})
					end
					return
				end
      case 'orientimages'
        if strcmpi(varargin{ii+1},'off')
          Dat.OrientImages = false;
				end
				
			case 'sortpss'
				 if strcmpi(varargin{ii+1},'off')
					 Dat.SortPSS = false;
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

% Don't use DC correction if number of transients is greater that 1
if procpar.nt>1
  Dat.DCcorrection = false;
end

%% Read phasetable if necessary
if Dat.Sorting
	
  % Look in preferences for tablib-directory
  if isfield(procpar,'petable') && ~strcmpi(procpar.petable{1},'n')
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
	else
		phasetable = [];
	end
		
  % If petable cannot be found, check if it is in procpar...
  if isempty(phasetable) && isfield(procpar,'pe_table')
    phasetable = {'t1',procpar.pe_table(:)};
  elseif isempty(phasetable) && isfield(procpar,'pelist') && ...
      ~isempty(procpar.pelist) && isnumeric(procpar.pelist) && ...
      length(procpar.pelist)>1
    phasetable = {'t1',reshape(procpar.pelist,procpar.etl,[]).'};
	end
  
	if ~isempty(phasetable)
		Dat.phasetable = phasetable{1,2};
	else
		Dat.phasetable = [];
	end
	
  %% Abort and throw error if phasetable cannot be read and the FID-file
  % has not been sorted
%   if isempty(phasetable) && ~isfield(procpar,'flash_converted')
%     DATA=[];
%     if ~Dat.ShowErrors
%       msg_out={['Could not find the required phase table "' ...
%                 procpar.petable{1} '" in the following folders'],...
%                ['"' tabpath '".']};
%     else
%       error('Could not find the required phase table "%s" in %s',...
%             procpar.petable{1},tabpath)
%     end
%     return
%   elseif ( isempty(phasetable) && isfield(procpar,'flash_converted') ) || ...
%       ~Dat.UsePhaseTable
%     %% If the FID-file has been sorted, the reading can continue but
%     % throw a warning.
%     fprintf(1,'Warning: aedes_readfid: Could not find phasetable "%s" in "%s"!\n',procpar.petable{1},tabpath)
%     Dat.phasetable=[];
%   else
%     Dat.phasetable = phasetable{1,2};
%   end
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
  [hdr,kspace,msg_out]=l_ReadBlockData(file_fid,hdr,Dat);
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
    DATA.FTDATA=[];
    DATA.KSPACE=[];
    DATA.PROCPAR=procpar;
    DATA.PHASETABLE=Dat.phasetable;
  end
end

% Close file
fclose(file_fid);


% Set file name and path to the HDR structure
DATA.HDR.fname=fname;
DATA.HDR.fpath=fpath;

% Set parameter values
if ~Dat.DCcorrection
  DATA.HDR.Param.DCcorrection = 'off';
else
  DATA.HDR.Param.DCcorrection = 'on';
end

% Reconstruct kspace =======================================
if Dat.ReturnRawKspace
  % If asked to return raw unsorted kspace, return immediately
  DATA.KSPACE=kspace;
  return
else
  
  % Check if the sequence used has a custom reconstruct code
  recon_func=l_GetReconFunc(recon_dir);
  recon_func_ind = 0;
  if ~isempty(recon_func)
    % Get sequence name
    seqfil = procpar.seqfil;
    
    % Check if any of the custom reconstruction files would
    % like to do the reconstruction...
    for ii=1:length(recon_func)
      try
        list=recon_func{ii}();
      catch
        warning('Error in custom reconsruction function "%s", skipping...',func2str(recon_func{ii}));
        continue
      end
      if any(strcmp(seqfil,list))
        recon_func_ind = ii;
        break
      end
    end
  end
  
  if recon_func_ind==0
    Dat.UseCustomRecon = false;
  end
  
  if Dat.UseCustomRecon
    if Dat.ShowWaitbar
      wbh=aedes_calc_wait(sprintf('%s\n%s',...
        ['Using custom function ',upper(func2str(recon_func{recon_func_ind}))],...
        ['to reconstruct sequence ',procpar.seqfil{1}]));
      drawnow
    end
    [kspace,data,msg_out]=recon_func{recon_func_ind}(kspace,Dat,procpar);
    if Dat.ShowWaitbar
      close(wbh)
    end
    if isempty(data) && Dat.ReturnFTData
      % Fourier transform data if not done in custom reconstruction code
      [data,msg_out]=l_CalculateFFT(kspace,Dat,procpar);
    end
  else
    % Use the default reconstruction code
    [kspace,msg_out,procpar]=l_ReconstructKspace(kspace,Dat,procpar);
    DATA.PROCPAR = procpar;
		
    % Fourier transform data
    if Dat.ReturnFTData
      [data,msg_out]=l_CalculateFFT(kspace,Dat,procpar);
    end
  end
end
if Dat.ReturnKSpace
  DATA.KSPACE = kspace;
else
  DATA.KSPACE = [];
end
if Dat.ReturnFTData
  DATA.FTDATA = data;
else
  DATA.FTDATA = [];
end
clear kspace data


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
function [hdr,kspace,msg_out]=l_ReadBlockData(file_fid,hdr,Dat)

hdr.BlockHeaders = [];
%hdr.HypercmplxBHeaders = [];
kspace=[];
count = [];
msg_out='';


BlockHeadStatusLabels = {...
  'S_DATA',...       % 0 = no data, 1 = data
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

BlockHeadModeLabels = {...
  'NP_PHMODE',...   % 1 = ph mode
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

% Allocate space for kspace
% kspace = complex(zeros(hdr.FileHeader.np/2,...
%           hdr.FileHeader.ntraces*hdr.FileHeader.nblocks,...
%           Dat.precision));

kspace = complex(zeros(hdr.FileHeader.np/2,...
	hdr.FileHeader.ntraces*hdr.FileHeader.nblocks,Dat.precision));


%% - The older robust (but also slower) way for reading the data.
%% When the blocksize is relatively small, this is also quite fast.
if ~Dat.FastDataRead
  
  % Initialize waitbar
  if Dat.ShowWaitbar
    wb_h = aedes_wbar(0/hdr.FileHeader.nblocks,...
      {['Reading VNMR data.'],...
      ['(Processing data block ' ...
      num2str(0) '/' num2str(hdr.FileHeader.nblocks) ')']});
  end

  %% Read data and block headers
  for ii=1:hdr.FileHeader.nblocks
    %% Update waitbar
    if Dat.ShowWaitbar
      aedes_wbar(ii/hdr.FileHeader.nblocks,...
        wb_h,...
        {['Reading VNMR data.'],...
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
    if ~Dat.DCcorrection
      complex_block = complex(tmp(1:2:end,:),tmp(2:2:end,:));
    else
      complex_block = complex(tmp(1:2:end,:)-hdr.BlockHeaders.lvl,...
        tmp(2:2:end,:)-hdr.BlockHeaders.tlt);
    end
    
    inds = ((ii-1)*hdr.FileHeader.ntraces+1):(ii*hdr.FileHeader.ntraces);
    kspace(:,inds) = complex_block;

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
      wb_h = aedes_calc_wait(['Reading VNMR data.']);
    else
      wb_h = aedes_wbar(1/nBlocks,{['Reading VNMR data.'],...
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
      aedes_wbar(ii/nBlocks,wb_h,{['Reading VNMR data.'],...
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
    tmp=reshape(tmp,hdr.FileHeader.np,[]);

    if ii==nBlocks
      inds = ((ii-1)*nbh*hdr.FileHeader.ntraces+1):size(kspace,2);
    else
      inds = ((ii-1)*nbh*hdr.FileHeader.ntraces+1):ii*nbh*hdr.FileHeader.ntraces;
    end

    % Do DC-correction if necessary
    if ~Dat.DCcorrection
      kspace(:,inds)=complex(tmp(1:2:end,:,:),tmp(2:2:end,:,:));
    else
      kspace(:,inds)=complex(tmp(1:2:end,:,:)-hdr.BlockHeaders.lvl,...
        tmp(2:2:end,:,:)-hdr.BlockHeaders.tlt);
    end
  end
  clear tmp

  if Dat.ShowWaitbar
    delete(wb_h)
    pause(0.1)
  end

end

% ==================================================================
% Reconstruct k-space
% ==================================================================
function [kspace,msg_out,procpar]=l_ReconstructKspace(kspace,Dat,procpar)

msg_out = '';

% Get the strategic parameters
seqcon = procpar.seqcon{1};
seqfil = procpar.seqfil{1};
np = procpar.np;
nv = procpar.nv;
nv2 = procpar.nv2;
nv3 = procpar.nv3;
ns = procpar.ns;
if isfield(procpar,'etl')
  etl = procpar.etl;
else
  etl = 1;
end
pss = procpar.pss;
nt = procpar.nt;
nf = procpar.nf;
ne = procpar.ne;
if isfield(procpar,'flash_converted')
  % Don't try to sort already sorted data... 
  Dat.Sorting = false;
  seqcon(2:3) = 'sc';
end
if isfield(procpar,'nseg')
	nseg = procpar.nseg;
else
	nseg = 1;
end


% Check dimensions
if nv==0
  % 1D-image (specrum)
  AcqType = 1;
elseif nv2==0
  % 2D-image
  AcqType = 2;
elseif nv3==0
  % 3D-image
  AcqType = 3;
else
  AcqType = 4;
end

% Check number of receivers
if isfield(procpar,'rcvrs') && length(procpar.rcvrs{1})>1
  nRcvrs = length(find(procpar.rcvrs{1}=='y'));
else
  nRcvrs = 1;
end


% Check for arrayed acquisition
if isfield(procpar,'array') && ~isempty(procpar.array{1})
  %if length(procpar.array)==1 && ~iscell(procpar.array{1}) && ...
  %    strcmp(procpar.array{1},'pad') && all(procpar.pad==0)
  %  % Skip the array parsing if the array is a dummy "pad"-array...
  %  isAcqArrayed = false;
  %  ArrayLength = 1;
  %else
    isAcqArrayed = true;
    ArrayLength = [];
    
    % Determine array length
    for ii=1:length(procpar.array)
      if iscell(procpar.array{ii})
        ArrayLength(ii) = length(procpar.(procpar.array{ii}{1}));
      else
        ArrayLength(ii) = length(procpar.(procpar.array{ii}));
      end
    end
    ArrayLength = prod(ArrayLength);
  %end
else
  isAcqArrayed = false;
  ArrayLength = 1;
end


% Reconstruct k-space ----------------------
if AcqType==1
  % Reconstruct 1D data ...
elseif AcqType==2
  
  % Reconstruct 2D data
  if strcmpi(seqcon,'ncsnn')
    kspace = reshape(kspace,[np/2,etl,ns,nRcvrs,ArrayLength,nv/etl]);
    kspace = permute(kspace,[1 2 6 3 5 4]);
    kspace = reshape(kspace,[np/2,nv,ns*ArrayLength,1,nRcvrs]);
  elseif strcmpi(seqcon,'nscnn')
    if isfield(procpar,'flash_converted')
      kspace = reshape(kspace,[np/2,nRcvrs,ArrayLength,ns,nv]);
      kspace = permute(kspace,[1 5 4 3 2]);
      kspace = reshape(kspace,[np/2,nv,ns*ArrayLength,1,nRcvrs]);
    else
      kspace = reshape(kspace,[np/2,nRcvrs,nv,ArrayLength,ns]);
      kspace = permute(kspace,[1 3 5 4 2]);
      kspace = reshape(kspace,[np/2,nv,ns*ArrayLength,1,nRcvrs]);
    end
  elseif strcmpi(seqcon,'nccnn')
    kspace = reshape(kspace,[np/2,etl,ns,nv/etl,ArrayLength,nRcvrs]);
    kspace = permute(kspace,[1 2 4 3 5 6]);
    kspace = reshape(kspace,[np/2,nv,ns*ArrayLength,1,nRcvrs]);
  end
elseif AcqType==3
  % Reconstruct 3D data
  if strcmpi(seqcon,'nscsn')
    kspace = reshape(kspace,[np/2,etl,nv/etl,nRcvrs,ArrayLength*nv2]);
    kspace = permute(kspace,[1 2 3 5 4]);
    kspace = reshape(kspace,[np/2,nv,nv2,ArrayLength,nRcvrs]);
	end
	if strcmpi(seqcon,'ccccn')
		kspace = reshape(kspace,[np/2,etl,nRcvrs,ne,nv/etl,ArrayLength*nv2]);
		kspace = permute(kspace,[1 2 5 6 4 3]);
		kspace = reshape(kspace,[np/2,nv,nv2,ne*ArrayLength,nRcvrs]);
	end
else
  % Reconstruct nD data ... to be written
end


% Sort data using phasetable ------------------
if Dat.Sorting && ~isempty(Dat.phasetable)
  Dat.phasetable = Dat.phasetable.';
  kspace(:,Dat.phasetable(:),:,:,:)=kspace;
end

% Sort interleaved 2D data using pss
[sorted_pss,I_pss]=sort(procpar.pss);
if ~ismember(procpar.pss,sorted_pss,'rows') && Dat.SortPSS
	if isAcqArrayed
		sz = size(kspace);
		kspace = reshape(kspace,sz(1),sz(2),length(I_pss),[]);
		kspace = kspace(:,:,I_pss,:,:);
		kspace = reshape(kspace,sz(1),sz(2),sz(3),[]);
	else
		kspace = kspace(:,:,I_pss,:,:);
	end
	procpar.pss_orig = procpar.pss;
	procpar.pss = sorted_pss;
end

% ==================================================================
% Fourier transform k-space
% ==================================================================
function [data,msg_out]=l_CalculateFFT(kspace,Dat,procpar)

msg_out = '';
data = [];


% Check dimensions
if procpar.nv==0
  % 1D-image
  AcqType = 1;
elseif procpar.nv2==0
  % 2D-image
  AcqType = 2;
elseif procpar.nv3==0
  % 3D-image
  AcqType = 3;
else
  AcqType = 4;
end

% Check number of receivers
if isfield(procpar,'rcvrs') && length(procpar.rcvrs{1})>1
  nRcvrs = length(find(procpar.rcvrs{1}=='y'));
else
  nRcvrs = 1;
end

% If zeropadding is requested, calculate the padded size
if Dat.ZeroPadding~=0
  if Dat.ZeroPadding==1
    % Zeropad to square
    if AcqType==1
      padSize = procpar.np/2;
    elseif AcqType==2
      padSize = ones(1,2)*procpar.np/2;
      padSize(3) = size(kspace,3);
    else
      padSize = ones(1,3)*procpar.np/2;
    end
  else
    % Zeropadding is on "auto", i.e. zeropad to FOV
    lpe = procpar.lpe;
    lpe2 = procpar.lpe2;
    lro = procpar.lro;
    if AcqType==2
      % 2D data
      padSize = [procpar.np/2 ...
        procpar.np/2*(lpe/lro) ...
        size(kspace,3)];
    elseif AcqType==3 && lpe2~=0
      % 3D data
      padSize = [procpar.np/2 ...
        procpar.np/2*(lpe/lro) ...
        procpar.np/2*(lpe2/lro)];
    end
  end
  ks_sz = [size(kspace,1) ...
    size(kspace,2) ...
    size(kspace,3)];
	padSize = round(padSize);
  if any(padSize>ks_sz)
    %kspace(padSize) = 0;
		kspace(padSize(1),padSize(2),padSize(3)) = 0;
  end
else
  padSize = [size(kspace,1) ...
    size(kspace,2) ...
    size(kspace,3)];
end

% Allocate space for Fourier transformed data
if nRcvrs>1 && any(Dat.SumOfSquares==[1 2])
  data_sz = [padSize,size(kspace,4),size(kspace,5)+1];
  data = zeros(data_sz,Dat.precision);
else
  data = zeros(size(kspace),Dat.precision);
end
%data = [];
%if strcmpi(Dat.precision,'single')
%  data = single(data);
%end

% Fourier transform data
if nRcvrs>1 && any(Dat.SumOfSquares==[1 2])
  ind = [2:size(data,5)];
else
  ind = [1:size(data,5)];
end
if AcqType==1
  data(:,:,:,:,ind) = abs(fftshift(ifft(kspace,[],1),1));
elseif AcqType==2
  data(:,:,:,:,ind) = abs(fftshift(fftshift(ifft(ifft(kspace,[],1),[],2),1),2));
elseif AcqType==3
  data(:,:,:,:,ind) = abs(fftshift(fftshift(fftshift(ifft(ifft(ifft(kspace,[],1),[],2),[],3),1),2),3));
end

% Calculate sum-of-squares image
if nRcvrs>1 && any(Dat.SumOfSquares==[1 2])
  % Calculate sum-of-squares
  data(:,:,:,:,1) = sqrt(sum(data(:,:,:,:,ind).^2,5));
  data=abs(data);
  if Dat.SumOfSquares==1
    % Remove individual receiver data
    data=data(:,:,:,:,1);
  end
end

% ==================================================================
% Find custom functions for VNMR k-space reconstruction
% ==================================================================
function recon_func=l_GetReconFunc(recon_dir)

recon_func = {};

if ~isdir(recon_dir)
  return
end

dir_struct=dir(recon_dir);
recon_files = {dir_struct(~[dir_struct(:).isdir]).name};
if isempty(recon_files)
  return
end

% Remove files that don't have .m extension
ind = regexpi(recon_files,'\.m$');
recon_files = {recon_files{~cellfun('isempty',ind)}};

currentDir = pwd;
try
  cd(recon_dir);
  for ii=1:length(recon_files)
    recon_func{ii}=str2func(recon_files{ii}(1:end-2));
  end
  cd(currentDir);
catch
  cd(currentDir);
end



