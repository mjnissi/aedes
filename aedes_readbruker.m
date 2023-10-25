function DATA = aedes_readbruker(filename,varargin)
% AEDES_READBRUKER - Read Bruker file formats (raw FID-files and 
%                    reconstructed 2DSEQ files) 
%   
%
% Synopsis: 
%        DATA=aedes_readbruker(filename,'PropertyName1',value1,'PropertyName2',value2,...)
%        DATA=aedes_readbruker(filename,'header')
%
% Description:
%        The function reads Bruker FID- and 2DSEQ files and returns a structure
%        DATA with fields DataFormat, HDR, FTDATA, KSPACE, PROCPAR, and
%        PHASETABLE. The fields of the DATA structure are constructed as
%        follows: 
%        
%        DATA
%          |-> DataFormat         (string identifier for data format 'bruker_raw'/'bruker_reco')
%          |-> HDR
%                |-> FileHeader   (data file header)
%                |-> BlockHeaders (Varian FID file data block headers)
%                |-> fname        (file name)
%                |-> fpath        (file path)
%                |-> Param        (parameter values used by AEDES_READBRUKER to read the data)
%          |-> FTDATA             (real fourier transformed image data)
%          |-> KSPACE             (complex k-space, empty by default)
%          |-> PROCPAR            (parameters from procpar file, NOTE: Varian FID files only)
%          |-> PHASETABLE         (phasetable)
%
%        The DATA structure is returned as empty (without the above described fields) 
%        if an error has occurred during reading.
%        
%        The first input argument is either a string containing full path to
%        the FID/2DSEQ-file or an empty string. Only the data file headers
%        can be read by giving a string 'header' as the second input
%        argument.
%
%        By default the k-space is not returned, i.e. the field KSPACE is
%        empty. The returned data can be adjusted by using the 'return'
%        property and values 1, 2, or 3 (see below for more information).
%        The return option is ignored for 2DSEQ files.
%
%        The supported property-value pairs in AEDES_READBRUKER (property strings
%        are not case sensitive):
%
%        Property:        Value:                Description:
%        *********        ******                ************
%        'Return'      :  [{1}| 2 | 3 | 4 | 5 ] % 1=return only ftdata (default)
%                                               % 2=return only k-space
%                                               % 3=return both ftdata & kspace
%                                               % 4=return raw kspace
%                                               % 5=return complex ftdata & kspace (SumOfSquares off)
%
%        'wbar'        : [ {'on'} | 'off' ]     % Show/hide waitbar
%
%        'Precision'   : [{'double'}|'single']  % The precision for the
%                                               % reconstructed and scaled
%                                               % data. (default='double')
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
%        'force_4d'    : [ {'on'} | 'off' ]     % Force data that has 5 or
%                                               % more dimensions to be 4D.
%                                               % Default 'on'. This only
%                                               % applies to 2dseq files.
%                                               % NOTE: Aedes
%                                               % cannot handle data with 5
%                                               % or more dimensions.
%
%        'SumOfSquares':  [ {1} | 2 ]           % 1=Return only the 
%                                               % sum-of-squares image
%                                               % for multireceiver
%                                               % data (default).
%                                               % 2=Return individual
%                                               % receiver data 
%                                               % NOTE: This property has 
%                                               % no effect on single 
%                                               % receiver data or when 
%                                               % reading 2DSEQ files.
%
%        NOTE: The support for reconstructing raw bruker data is
%        EXPERIMENTAL and should not normally be used. If you don't need
%        kspace information USE THE 2DSEQ FILES INSTEAD OF THE FID FILES!
%
% Examples:
%        DATA=aedes_readbruker(filename)   % Read image data from 'filename'
%
% See also:
%        AEDES, AEDES_READVNMR

% This function is a part of Aedes - A graphical tool for analyzing 
% medical images
%
% Copyright (C) 2011 Juha-Pekka Niskanen <Juha-Pekka.Niskanen@uef.fi>
% 
% Department of Applied Physics, Department of Neurobiology
% University of Eastern Finland, Kuopio, FINLAND
%
% This program may be used under the terms of the GNU General Public
% License version 2.0 as published by the Free Software Foundation
% and appearing in the file LICENSE.TXT included in the packaging of
% this program.
%
% This program is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
% WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
% 
% 
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EDITS, Tuomainen, TV
%Department of Applied Physics, University of Eastern Finland 2020

%Added check with number of receivers to number of repetitions on line ~1030
%To ensure that RAREVTR data is correctly read (not
%averaged over number of reps), tvt2020

%changed 3 to length(ks_sz) TVT2020, and 1:length(ks_sz) to padSize
%on line  ~~998 changed
%from     any(padSize>ks_sz(1:3))
%to       any(padSize(1:min(length(ks_sz),3))>ks_sz(1:min(length(ks_sz),3)))
%the 1: in 1:min(length(ks_sz),3)) was added such that MGE acq can be read
%without the 1:, the RAREVTR acq were read correctly tho.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Defaults
DATA = [];
isRawData = true;
Dat.showWbar = true;
Dat.precision = 'double';
Dat.return = 3;
Dat.zerofill = 'auto';
Dat.SumOfSquares = 1;
Dat.Force4D = true;
Dat.UseCustomRecon = true;

% Default custom reconstruction file directory
tmp = which('aedes');
if ~isempty(tmp)
  tmp_path=fileparts(tmp);
  recon_dir = [tmp_path,filesep,'bruker_recon',filesep];
else
  % If all else fails, look in the current directory
  recon_dir = [pwd,filesep,'bruker_recon',filesep];
end

if nargin == 0 || isempty(filename)
	[fname,fpath,findex] = uigetfile({'fid;fid.*;2dseq','Bruker file formats (fid, 2dseq)';...
		'*.*','All Files (*.*)'},'Select Bruker file formats');
	if isequal(fname,0)
		% Canceled
		return
	end
	filename = [fpath,fname];
end

% Parse filename
[fp,fn,fe] = fileparts(filename);
if strcmpi(fn,'fid')
	isRawData = true;
	data_format = aedes_getdataformat(filename);
	if strcmpi(data_format,'vnmr')
		error('aedes_readbruker cannot read Varian FID-files. Use aedes_readvnmr instead.');
	end
elseif strcmpi(fn,'2dseq')
	isRawData = false;
else
	error('Only Bruker FID and 2DSEQ files are supported.');
end

%keyboard;

% Check varargin length and check if header only requested
if strcmpi(varargin{1},'header')
    % read header only
    if isRawData
    	[hdr,msg] = l_ReadHeaderFid(filename);
    	if isempty(hdr)
    		error(msg);
    	end

    else
    	[hdr,msg] = l_ReadHeader2dseq(filename);
    	if isempty(hdr)
    		error(msg);
    	end
    end
    DATA=hdr;
    return;

elseif rem(length(varargin),2)==1
    error('Input arguments may only include filename and param/value pairs.')
end

% Parse varargin
for ii=1:2:length(varargin)
	param = varargin{ii};
	if ~ischar(param) || isempty(param)
		error('Parameters must be strings.')
	end
	value = varargin{ii+1};
	switch lower(param)
		case 'wbar'
			if strcmpi(value,'on')
				Dat.wbar = true;
			elseif strcmpi(value,'off')
				Dat.wbar = false;
			elseif islogical(value)
				Dat.wbar = value;
			else
				error('Value for WBAR parameter can be either ''on'' or ''off''.')
			end
		case 'precision'
			Dat.precision = value;
		case 'sumofsquares'
			if ~isscalar(value) || ~any(value==[1 2])
				error('SumOfSquares property can be 1 (return sum-of-squares) or 2 (return individual receiver data).');
			end
			Dat.SumOfSquares = value;
		case 'zeropadding'
			 if ischar(varargin{ii+1}) && any(strcmpi(value,{'auto','on','off'}))
          if strcmpi(varargin{ii+1},'on')
            Dat.zerofill = lower(varargin{ii+1}); % on
          elseif strcmpi(varargin{ii+1},'off')
            Dat.zerofill = lower(varargin{ii+1}); % off
					else
						Dat.zerofill = lower(varargin{ii+1}); % auto
					end
			 else
				 error('ZeroPadding property has to be a string ''auto'', ''on'' or ''off''.')
			 end
		case 'return'
			if ~isnumeric(value) || ~isscalar(value) || isempty(value) || ...
					isnan(value) || ~isreal(value) || (value-floor(value)) ~= 0 || ...
					value < 1 || value > 5
				error('Invalid return value. Return value can be an integer in the range 1..5')
			end
			Dat.return = value;
            if value==5
                Dat.SumOfSquares = 2;
                warning('Complex image data requested, SumOfSquares turned off (value 2)')
            end
		case 'force_4d'
			if strcmpi(value,'on')
				Dat.Force4D = true;
			elseif strcmpi(value,'off')
				Dat.Force4D = false;
			elseif islogical(value)
				Dat.Force4D = value;
			else
				error('Value for Force_4D parameter can be either ''on'' or ''off''.')
			end
		otherwise
			error('Unknown parameter "%s".',param)
	end
end

if isRawData
	[hdr,msg] = l_ReadHeaderFid(filename);
	if isempty(hdr)
		error(msg);
	end
	
	% Check if the sequence used has a custom reconstruct code
  recon_func=l_GetReconFunc(recon_dir);
  recon_func_ind = 0;
  if ~isempty(recon_func)
    % Get sequence name
    pulprog = hdr.acqp.PULPROG;
    
    % Check if any of the custom reconstruction files would
    % like to do the reconstruction...
    for ii=1:length(recon_func)
      try
        list=recon_func{ii}();
      catch
        warning('Error in custom reconsruction function "%s", skipping...',func2str(recon_func{ii}));
        continue
      end
      if any(strcmp(pulprog,list))
        recon_func_ind = ii;
        break
      end
    end
  end
  
  if recon_func_ind==0
    Dat.UseCustomRecon = false;
	end
	
	% Read data from the fid file
	[data,kspace,pe1_table,pe2_table,msg] = l_ReadDataFid(filename,hdr,Dat);
	if ~isempty(msg)
		error(msg)
	end
	DATA.DataFormat = 'bruker_raw';
	
  if Dat.UseCustomRecon
    if Dat.showWbar
      wbh=aedes_calc_wait(sprintf('%s\n%s',...
        ['Using custom function ',upper(func2str(recon_func{recon_func_ind}))],...
        ['to reconstruct data measured with pulse program ',pulprog]));
      drawnow
		end
		
    [kspace,data,msg]=recon_func{recon_func_ind}(kspace,Dat,hdr);
    if Dat.showWbar
      close(wbh)
		end
	end
else
	[hdr,msg] = l_ReadHeader2dseq(filename);
	if isempty(hdr)
		error(msg);
	end
	[data,kspace,msg] = l_ReadData2dseq(filename,hdr,Dat);
	DATA.DataFormat = 'bruker_reco';
	pe1_table = [];
	pe2_table = [];
end
if isempty(data) && isempty(kspace)
	DATA = [];
	error(msg);
end

DATA.HDR.FileHeader = hdr;
DATA.HDR.BlockHeader = [];
DATA.HDR.fname = [fn,fe];
DATA.HDR.fpath = [fp,filesep];
if Dat.return == 5
    DATA.FTDATA = abs(data);
    DATA.FTDATACOMPLEX = data;
else
    DATA.FTDATA = data;
end
if Dat.return == 1
	DATA.KSPACE = [];
	clear kspace;
else
	DATA.KSPACE = kspace;
end
DATA.PROCPAR = [];
DATA.PHASETABLE = {pe1_table,pe2_table};


% - Subfunctions -

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read header files for raw data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hdr,msg] = l_ReadHeaderFid(filename)

hdr = [];
msg = '';

% Parse filename
[fp,fn,fe] = fileparts(filename);

% Read Header files -------------------------

% Read ACQP file
acqp_file = [fp,filesep,'acqp'];
if exist(acqp_file,'file')~=2
	msg = sprintf('Cannot find ACQP file in "%s".',fp);
	return
end
try
	hdr.acqp = aedes_readjcamp(acqp_file);
catch
	hdr = [];
	msg = sprintf('Error while reading "%s". The error was "%s".',...
		acqp_file,lasterr);
	return
end

% Read Method file
method_file = [fp,filesep,'method'];
if exist(method_file,'file')~=2
	msg = sprintf('Cannot find METHOD file in "%s".',fp);
	return
end
try
	hdr.method = aedes_readjcamp(method_file);
catch
	hdr = [];
	msg = sprintf('Error while reading "%s". The error was "%s".',...
		method_file,lasterr);
	return
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read header files for reconstructed 2dseq data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hdr,msg] = l_ReadHeader2dseq(filename)

hdr = [];
msg = '';

% Parse filename
[fp,fn,fe] = fileparts(filename);

% Read Header files -------------------------

% Read visu_pars file
visu_pars_file = [fp,filesep,'visu_pars'];
if exist(visu_pars_file,'file')==2
	try
		hdr.visu_pars = aedes_readjcamp(visu_pars_file);
	catch
		hdr = [];
		msg = sprintf('Error while reading "%s". The error was "%s".',...
			visu_pars_file,lasterr);
		return
	end
else
	% Warn if visu_pars file was not found
	warning('AEDES_READBRUKER:VisuParsNotFound',...
		'VISU_PARS file was not not found in "%s". The 2dseq file might not be read correctly.',fp);
end

% Read d3proc file. This is however deprecated and the Visu parameters
% should be used instead...
d3proc_file = [fp,filesep,'d3proc'];
if exist(d3proc_file,'file')==2
	try
		hdr.d3proc = aedes_readjcamp(d3proc_file);
	catch
		warning('Could not read D3PROC file "%s".',d3proc_file);
	end
end

% Read reco file
reco_file = [fp,filesep,'reco'];
if exist(reco_file,'file')==2
	try
		hdr.reco = aedes_readjcamp(reco_file);
	catch
		hdr = [];
		msg = sprintf('Error while reading "%s". The error was "%s".',...
			reco_file,lasterr);
		return
	end
end

% Read procs file
procs_file = [fp,filesep,'procs'];
if exist(procs_file,'file')==2
	try
		hdr.procs = aedes_readjcamp(procs_file);
	catch
		warning('Could not read PROCS file "%s".',procs_file);
	end
end

% Read also acqp and method files if they exist
ind = find(fp==filesep);
if length(ind) >= 2
	fp_fid = fp(1:ind(end-1));
	
	% Try to read ACQP file
	acqp_file = [fp_fid,filesep,'acqp'];
	if exist(acqp_file,'file')==2
		hdr.acqp = aedes_readjcamp(acqp_file);
	end
	
	% Try to read Method file
	method_file = [fp_fid,filesep,'method'];
	if exist(method_file,'file')==2
		hdr.method = aedes_readjcamp(method_file);
	end
end

% Check that either VISU_PARS, RECO or D3PROC file is found
if ~isfield(hdr,'visu_pars') && ~isfield(hdr,'reco') && ...
		~isfield(hdr,'d3proc')
	hdr = [];
	msg = 'Could not read VISU_PARS, RECO or D3PROC files.';
	return
end

% Check for Functional Imaging Tool files
if exist([fp,filesep,'fun'],'dir') == 7
	
	% To be implemented...
	
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read raw Bruker data (fid files) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,kspace,pe1_table,pe2_table,msg] = l_ReadDataFid(filename,hdr,Dat)

data = [];
kspace = [];
msg = '';

% Get word type
if strcmpi(hdr.acqp.GO_raw_data_format,'GO_32BIT_SGN_INT')
	precision = ['int32=>',Dat.precision];
elseif strcmpi(hdr.acqp.GO_raw_data_format,'GO_16BIT_SGN_INT')
	precision = ['int16=>',Dat.precision];
elseif strcmpi(hdr.acqp.GO_raw_data_format,'GO_16BIT_FLOAT')
	precision = ['single=>',Dat.precision];
else
	msg = sprintf('Unknown word type "%s".',hdr.acqp.GO_raw_data_format);
	return
end

% Get byteorder
if strcmpi(hdr.acqp.BYTORDA,'little')
	byteorder = 'ieee-le';
else
	byteorder = 'ieee-be';
end

% Get block size
if strcmpi(hdr.acqp.GO_block_size,'continuous')
	fixed_blocksize = false;
elseif strcmpi(hdr.acqp.GO_block_size,'Standard_KBlock_format')
	fixed_blocksize = true;
else
	msg = sprintf('Unknown block size identifier "%s".',hdr.acqp.GO_block_size);
	return
end

% Open fid file for reading
fid = fopen(filename,'r',byteorder);
if fid < 0
	hdr.acqp = [];
	hdr.method = [];
	msg = sprintf('Could not open file "%s" for reading.',filename);
	return
end

% Read fid file
tmp = fread(fid,inf,precision);
fclose(fid);
kspace = complex(tmp(1:2:end),tmp(2:2:end));

% Handle block size
if strcmp(hdr.acqp.GO_block_size,'Standard_KBlock_Format')
	blk_size = round(length(kspace)/1024)*1024;
else
	blk_size = hdr.acqp.ACQ_size(1)/2;
end
kspace = reshape(kspace,blk_size,[]);

% Reconstruct k-space
if Dat.return ~= 4 && ~Dat.UseCustomRecon
	% Issue a warning that reconstruction of raw data is experimental
	warning('AEDES_READBRUKER:RawDataReconstructionExperimental',...
		'The reconstruction of raw bruker data is EXPERIMENTAL and may not work correctly.');
	
	[kspace,pe1_table,pe2_table,msg] = l_ReconstructKspace(kspace,hdr,Dat);
else
	[pe1_table,pe2_table] = l_GetPeTable(hdr);
end
if isempty(kspace)
	kspace = [];
	return
end

% Fourier transform k-space
if any(Dat.return==[1 3 5]) && ~Dat.UseCustomRecon
	[data,msg,kspace] = l_doFFT(kspace,hdr,Dat);
	if isempty(data)
		kspace = [];
		data = [];
		return
	end
else
	data = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read reconstructed Bruker data (2dseq files) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,kspace,msg] = l_ReadData2dseq(filename,hdr,Dat)

data = [];
kspace = [];
msg = '';

% Check which files to use
isVisuPars = false;
isReco = false;
isD3proc = false;
if isfield(hdr,'visu_pars')
	isVisuPars = true;
end
if isfield(hdr,'reco')
	isReco = true;
end
if isfield(hdr,'d3proc')
	isD3proc = true;
end

% Get word type
if isVisuPars
	% Use visu_pars information first
	if strcmpi(hdr.visu_pars.VisuCoreWordType,'_8BIT_UNSGN_INT')
		precision = ['uint8=>',Dat.precision];
	elseif strcmpi(hdr.visu_pars.VisuCoreWordType,'_16BIT_SGN_INT')
		precision = ['int16=>',Dat.precision];
	elseif strcmpi(hdr.visu_pars.VisuCoreWordType,'_32BIT_SGN_INT')
		precision = ['int32=>',Dat.precision];
	elseif strcmpi(hdr.visu_pars.VisuCoreWordType,'_32BIT_FLOAT')
		precision = ['single=>',Dat.precision];
	else
		msg = sprintf('Unknown word type "%s".',hdr.visu_pars.VisuCoreWordType);
		return
	end
elseif isReco
	% If visu_pars was not found, use reco information
	if strcmpi(hdr.reco.RECO_wordtype,'_8BIT_UNSGN_INT')
		precision = ['uint8=>',Dat.precision];
	elseif strcmpi(hdr.reco.RECO_wordtype,'_16BIT_SGN_INT')
		precision = ['int16=>',Dat.precision];
	elseif strcmpi(hdr.reco.RECO_wordtype,'_32BIT_SGN_INT')
		precision = ['int32=>',Dat.precision];
	elseif strcmpi(hdr.reco.RECO_wordtype,'_32BIT_FLOAT')
		precision = ['single=>',Dat.precision];
	else
		msg = sprintf('Unknown word type "%s".',hdr.reco.RECO_wordtype);
		return
	end
else
	% No VISU_PARS or RECO information available. Fall back to the depricated
	% D3PROC information and try to guess the word type from file size...
	sz = [hdr.d3proc.IM_SIX,...
		hdr.d3proc.IM_SIY,...
		hdr.d3proc.IM_SIZ,...
		hdr.d3proc.IM_SIT];
	d=dir(filename);
	nBytes = d.bytes;
	bytesPerPixel = nBytes/prod(sz);
	if bytesPerPixel == 1
		precision = ['uint8=>',Dat.precision];
	elseif bytesPerPixel == 2
		precision = ['int16=>',Dat.precision];
	else
		precision = ['int32=>',Dat.precision];
	end
end

% Get byteorder
if isVisuPars
	if strcmpi(hdr.visu_pars.VisuCoreByteOrder,'littleEndian')
		byteorder = 'ieee-le';
	else
		byteorder = 'ieee-be';
	end
elseif isReco
	if strcmpi(hdr.reco.RECO_byte_order,'littleEndian')
		byteorder = 'ieee-le';
	else
		byteorder = 'ieee-be';
	end
else
	% Use little endian if no information is available
	byteorder = 'ieee-le';
end

% Get image type
if isVisuPars
	if isfield(hdr.visu_pars,'VisuCoreFrameType')
		if strcmpi(hdr.visu_pars.VisuCoreFrameType,'COMPLEX_IMAGE')
			isDataComplex = true;
		elseif strcmpi(hdr.visu_pars.VisuCoreFrameType,'MAGNITUDE_IMAGE')
			isDataComplex = false;
        elseif strcmpi(hdr.visu_pars.VisuCoreFrameType,'REAL_IMAGE')
            isDataComplex = false;
		else
			msg = sprintf('Unknown image type "%s".',hdr.visu_pars.VisuCoreFrameType);
			return
		end
	else
		% If VisuCoreFrameType parameter was not found in visu_pars, assume
		% MAGNITUDE...
		isDataComplex = false;
	end
elseif isReco
	if strcmpi(hdr.reco.RECO_image_type,'COMPLEX_IMAGE')
		isDataComplex = true;
	elseif strcmpi(hdr.reco.RECO_image_type,'MAGNITUDE_IMAGE')
		isDataComplex = false;
	else
		msg = sprintf('Unknown image type "%s".',hdr.reco.RECO_image_type);
		return
	end
else
	% Assume magnitude image if information is not available
	isDataComplex = false;
end

% Open 2dseq file for reading
fid = fopen(filename,'r',byteorder);
if fid < 0
	msg = sprintf('Could not open 2dseq file "%s" for reading.');
	return
end

% Show waitbar
if Dat.showWbar
	wbh = aedes_calc_wait('Reading data from Bruker 2dseq file.');
end

% Check if word type is float
isWordFloat = false;
if isVisuPars && strcmpi(hdr.visu_pars.VisuCoreWordType,'_32BIT_FLOAT')
	isWordFloat = true;
elseif isReco && strcmpi(hdr.reco.RECO_wordtype,'_32BIT_FLOAT')
	isWordFloat = true;
end

% - Read data ---------------------------------
[data,count] = fread(fid,inf,precision);
fclose(fid);

% Try to get image dimensions
im_size = ones(1,4);
if isVisuPars
	for ii=1:length(hdr.visu_pars.VisuCoreSize)
		im_size(ii) = hdr.visu_pars.VisuCoreSize(ii);
	end
	im_size(length(hdr.visu_pars.VisuCoreSize)+1) = hdr.visu_pars.VisuCoreFrameCount;
else
	% Fall back to the depricated d3proc information, if visu_pars
	% information is not available...
	im_size = [hdr.d3proc.IM_SIX,...
		hdr.d3proc.IM_SIY,...
		hdr.d3proc.IM_SIZ,...
		hdr.d3proc.IM_SIT];
end

% Reshape to correct size
data = reshape(data,im_size);

% Map integer data to single or double values
slope = [];
offset = [];
if isfield(hdr,'visu_pars') && ~isWordFloat
	% The slope values in visu_pars are multiplicative inverses of the slope
	% values in RECO file...
	slope = hdr.visu_pars.VisuCoreDataSlope;
	offset = hdr.visu_pars.VisuCoreDataOffs;
	nDim = hdr.visu_pars.VisuCoreDim;
elseif isfield(hdr,'reco') && ~isWordFloat
	% Use RECO information if VISU_PARS is not available
	slope = 1/hdr.reco.RECO_map_slope;
	offset = hdr.reco.RECO_map_offset;
	nDim = length(hdr.reco.RECO_size);
elseif ~isWordFloat
	warning('AEDES_READBRUKER:ScalingInformationNotFound',...
		'Scaling information for 2dseq data was not found. Integer data is not scaled.');
end

% Calibrate data
if ~isempty(slope) && ~isempty(offset)
	for ii=1:length(slope)
		if any(nDim == [1 2])
			data(:,:,ii) = data(:,:,ii)*slope(ii)+offset(ii);
		else
			data(:,:,:,ii) = data(:,:,:,ii)*slope(ii)+offset(ii);
		end
	end
end

% Reshape frames using VisuFGOrderDesc information
if isVisuPars && isfield(hdr.visu_pars,'VisuFGOrderDesc')
	if hdr.visu_pars.VisuCoreFrameCount > 1 && ...
			size(hdr.visu_pars.VisuFGOrderDesc,1)>1
		nTotalDims = hdr.visu_pars.VisuCoreDim+hdr.visu_pars.VisuFGOrderDescDim;
		if size(hdr.visu_pars.VisuFGOrderDesc,1)==1
			frame_size = [hdr.visu_pars.VisuCoreSize,...
				hdr.visu_pars.VisuFGOrderDesc{1,1}];
		else
			FG_size = [hdr.visu_pars.VisuFGOrderDesc{:,1}];
			frame_size = [hdr.visu_pars.VisuCoreSize,FG_size];
			
			% Force data to be 4D...
			if Dat.Force4D
				frame_size = [frame_size(1:3) prod(frame_size(4:end))];
			end
		end
		data = reshape(data,frame_size);
	end
end

% Permute to correct orientation
%data = permute(data,[2 1 3 4]);
nDim = ndims(data);
data = permute(data,[2 1 3:nDim]);

if isDataComplex	
	kspace = complex(data(:,:,:,1),data(:,:,:,2));
	data = [];
end

% Close waitbar
if Dat.showWbar
	close(wbh);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruct k-space 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [kspace,pe1_table,pe2_table,msg] = l_ReconstructKspace(kspace,hdr,Dat)
%
% NOTE: The reconstruction of raw bruker data is EXPERIMENTAL and may work
% correctly.
%


msg = '';

% Don't try to reconstruct 1D data
nDims = hdr.acqp.ACQ_dim;
if nDims == 1
	pe1_table = l_GetPeTable(hdr);
	pe2_table = [];
	return
end

% Handle EPI data
if strfind(lower(hdr.acqp.PULPROG),'epi')    
	kspace = [];
	msg = 'EPI reconstruction is not supported.';
	return
	%epi_matrix_size = hdr.method.PVM_Matrix;
	%kspace = reshape(kspace,epi_matrix_size(1),epi_matrix_size(2),NSLICES,NR);
end

NI = hdr.acqp.NI; % Number of objects
NSLICES = hdr.acqp.NSLICES; % Number of slices
NR = hdr.acqp.NR; % Number of repetitions
NECHOES = hdr.acqp.NECHOES; % Number of echoes
phase_factor = hdr.acqp.ACQ_phase_factor; % scans belonging to a single image
im_size = hdr.acqp.ACQ_size;im_size(1)=im_size(1)/2;
order = hdr.acqp.ACQ_obj_order;

% Get number of receivers
if isfield(hdr,'method') && isfield(hdr.method,'PVM_EncNReceivers')
	nRcvrs = hdr.method.PVM_EncNReceivers;
else
	nRcvrs = 1; 
end

% Get phase table
usePeTable = true;
[pe1_table,pe2_table] = l_GetPeTable(hdr);
if isempty(pe1_table) && isempty(pe2_table)
	usePeTable = false;
end


% Reshape data so that all echoes are in correct planes
if nDims == 2
	% Handle 2D data
    rosize = numel(kspace)/nRcvrs/NI/im_size(2)/NR;
	kspace = reshape(kspace,rosize,...
		nRcvrs,phase_factor,NI,im_size(2)/phase_factor,NR);
	kspace = permute(kspace,[1 3 5 4 6 2]);
	kspace = reshape(kspace,rosize,im_size(2),NI,NR,nRcvrs);
    kspace = kspace(1:im_size(1),:,:,:,:,:);
elseif nDims == 3
	% Handle 3D data
    rosize = numel(kspace)/nRcvrs/NI/im_size(2)/im_size(3)/NR;
	kspace = reshape(kspace,rosize,...
		nRcvrs,phase_factor,NI,im_size(2)/phase_factor,im_size(3),NR);
	kspace = permute(kspace,[1 3 5 6 4 7 2]);
	kspace = reshape(kspace,rosize,im_size(2),im_size(3),NI,NR,nRcvrs);
    kspace = kspace(1:im_size(1),:,:,:,:,:);
else
	kspace = [];
	msg = sprintf('The number of dimensions "%d" is not supported.',nDims);
	return
end

% Sort echoes
if usePeTable
	if ~isempty(pe1_table)
		kspace = kspace(:,pe1_table,:,:,:);
	end
	if nDims == 3 && ~isempty(pe2_table)
		kspace = kspace(:,:,pe2_table,:,:);
	end
end

% Sort object order
if nDims == 2
	kspace(:,:,order+1,:,:) = kspace;
end

% Permute to correct orientation
if nDims == 2
    kspace = flipdim(flipdim(kspace,1),2);
elseif nDims == 3
    kspace = flipdim(flipdim(flipdim(kspace,1),2),3);
end

% Separate slices in case of 2D
if nDims == 2
    kspace = reshape(kspace,im_size(1),im_size(2),NI/NSLICES,NSLICES,NR,nRcvrs);
    kspace = permute(kspace,[1 2 4 3 5 6]);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check that the required reconstruction parameters are found in acqp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ok,msg] = l_CheckAcqpParams(hdr)

ok = false;
msg = '';

% Minimal list of ACQP parameters for data reconstruction
minParamList = {...
	'ACQ_dim',...
	'ACQ_size',...
	'NI',...
	'NR',...
	'BYTORDA',...
	'ACQ_word_size',...
	'ACQ_scan_size',...
	'ACQ_scan_shift',...
	'ACQ_phase_factor',...
	'ACQ_rare_factor',...
	'ACQ_phase_encoding_mode',...
	'ACQ_phase_enc_start',...
	'ACQ_spatial_size_1',...
	'ACQ_spatial_size_2',...
	'ACQ_spatial_phase_1',...
	'ACQ_spatial_phase_2',...
	'ACQ_obj_order',...
	'ACQ_mod',...
	'DSPFVS',...
	'DSPFIRM',...
	'DECIM'};

param = fieldnames(hdr.acqp);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourier transform data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,msg,kspace] = l_doFFT(kspace,hdr,Dat)

data = [];
msg = '';

% Check dimensions
if hdr.acqp.ACQ_dim == 1
  % 1D-image
  AcqType = 1;
elseif hdr.acqp.ACQ_dim == 2
  % 2D-image
  AcqType = 2;
elseif hdr.acqp.ACQ_dim == 3
  % 3D-image
  AcqType = 3;
else
  AcqType = 4;
end

% Get number of receivers
nRcvrs = size(kspace,5);
%THIS MAY BE ACTUALLY PROBLEMATIC AS IN RAREVTR THE NR IS CONNECTED TO THIS
%VALUE! Tuomainen TV 2020 (fixed in line ~1030, nRcvrs ~= NR check)

%keyboard;


% Do Zero filling if requested
if any(AcqType == [2 3])
	switch Dat.zerofill
		case 'auto'
			% Zeropadding is on "auto", i.e. zeropad to FOV
			fov1 = hdr.acqp.ACQ_fov(1);
			fov2 = hdr.acqp.ACQ_fov(2);
			if AcqType == 3
				fov3 = hdr.acqp.ACQ_fov(3);
			end
			if AcqType==2
				% 2D data
				padSize = round([hdr.acqp.ACQ_size(1)/2 ...
					hdr.acqp.ACQ_size(1)/2*(fov2/fov1) ...
					size(kspace,3)]);
			elseif AcqType==3
				% 3D data
				padSize = round([hdr.acqp.ACQ_size(1)/2 ...
					hdr.acqp.ACQ_size(1)/2*(fov2/fov1) ...
					hdr.acqp.ACQ_size(1)/2*(fov3/fov1)]);
            end

            
          %  keyboard;
            % MJN 01/2021: check if padSize vs. kspace size are "almost the same" and
            % don't pad if they are. Let's say, if difference is under 1%
            ksz=size(kspace);
            ksz=ksz(1:numel(padSize));
            if (1-mean(ksz./padSize))<0.01
                disp('Note -- padsize vs. voxel size almost equal -- NOT PADDING');
                padSize=ksz;
            end
                
		case 'on'
			% Zeropad to square
			if AcqType==1
				padSize = round(hdr.acqp.ACQ_size(1)/2);
			elseif AcqType==2
				padSize = round(ones(1,2)*hdr.acqp.ACQ_size(1)/2);
				padSize(3) = size(kspace,3);
			else
				padSize = ones(1,3)*hdr.acqp.ACQ_size(1)/2;
			end
		case 'off'
			padSize = [size(kspace,1) ...
				size(kspace,2) ...
				size(kspace,3)];
		otherwise
			warning('Unknown zerofill option "%s". Skipping zerofilling...',...
				Dat.zerofill);
			padSize = [size(kspace,1) ...
				size(kspace,2) ...
				size(kspace,3)];
	end
	
	% Pad with zeros
	% ks_sz = [size(kspace,1) size(kspace,2) size(kspace,3)];
    ks_sz = size(kspace);
	if any(padSize(1:min(length(ks_sz),3))>ks_sz(1:min(length(ks_sz),3))) %TVT2020
		% kspace(padSize(1),padSize(2),padSize(3))=0;
        
        ks_sz1 = ks_sz;
        ks_sz1(1:3) = padSize;
        halfpad = (ks_sz1(1:3)-ks_sz(1:3))/2;
        kspace1 = complex(zeros(ks_sz1));
        kspace1(1+floor(halfpad(1)):padSize(1)-ceil(halfpad(1)), ...
            1+floor(halfpad(2)):padSize(2)-ceil(halfpad(2)), ...
            1+floor(halfpad(3)):padSize(3)-ceil(halfpad(3)),:,:,:,:,:) = kspace;
        kspace = kspace1;
        clear kspace1
	end
end

%keyboard;

% Shift FOV to correct position using Fourier shift theorem
FOV = hdr.method.PVM_Fov;
matrix = hdr.method.PVM_Matrix;
FOVcenter = [hdr.method.PVM_SPackArrReadOffset hdr.method.PVM_SPackArrPhase1Offset hdr.method.PVM_SPackArrSliceOffset];
FOVcenter = FOVcenter(1:length(FOV));
shift = matrix.*FOVcenter./FOV; % FOV shift in pixels
for ii = 2:length(FOV) % Shift not necessary along read out direction, already accounted for in recorded data
    if ii == 3 % Sign change appears necessary for 3rd dimension
        shift(ii) = -shift(ii);
    end
    mod1 = ((1:matrix(ii))'-1-floor(matrix(ii)/2))/matrix(ii)*2*pi*round(shift(ii));
    mod1 = exp(1i*mod1);
    sizeks = size(kspace);
    sizeks(ii) = 1;
    kpermute = 1:length(sizeks);
    kpermute(1) = ii;
    kpermute(ii) = 1;
    mod1 = permute(mod1,kpermute);
    mod1 = repmat(mod1,sizeks);
    kspace = kspace.*mod1;
end


% Allocate space for Fourier transformed data
if nRcvrs>1 && Dat.SumOfSquares==1
    if hdr.acqp.NR ~= nRcvrs %TVT2020
  data_sz = [padSize,size(kspace,4),size(kspace,5)];
  data = zeros(data_sz,Dat.precision);
    else %TVT2020,  if the nRcvrs is mistaken in RAREVTR to be multiple receivers,
        %this has to be compared to the NR (repetitions, which now are with VTR and not averages of same signal...)
        %(so if nRcvrs == NR)
        data_sz = [padSize,size(kspace,5),size(kspace,4)];
        data = zeros(data_sz,Dat.precision);
    end
else
  data = zeros(size(kspace),Dat.precision);
end


if AcqType==1
  data(:,:,:,:,1:size(data,5)) = fftshift(fft(ifftshift(kspace,1),[],1),1);
  % Output magnitude image unless complex image requested
  if Dat.return~=5
      data = abs(data);
  end
elseif AcqType==2
  data(:,:,:,:,1:size(data,5)) = fftshift(fftshift(fft(fft(ifftshift(ifftshift(kspace,1),2),[],1),[],2),1),2);
  % Output magnitude image unless complex image requested
  if Dat.return~=5
      data = abs(data);
  end
elseif AcqType==3
  data(:,:,:,:,1:size(data,5)) = fftshift(fftshift(fftshift(fft(fft(fft(ifftshift(ifftshift(ifftshift(kspace,1),2),3),[],1),[],2),[],3),1),2),3);
  % Output magnitude image unless complex image requested
  if Dat.return~=5
      data = abs(data);
  end
end

%if A L L  E C H O E S  are collected with bruker MGE, then the whole thing
%is flipped with every-other volume. Thus this needs to be fixed;

if isfield(hdr.method,'EchoAcqMode')
if strcmpi(hdr.method.EchoAcqMode,'allEchoes')
    warning('aedes_readbruker.m: allEchoes in use, check if volumes are correct');
    data(:,:,:,2:2:end) = data(end:-1:1,:,:,2:2:end);
else
    warning('aedes_readbruker.m: EchoAcqMode exists, check if volumes are correct');
end
end

% Calculate sum-of-squares image
if nRcvrs>1 && Dat.SumOfSquares==1
  % Calculate sum-of-squares
  data = sqrt(sum(data(:,:,:,:,1:size(data,5)).^2,5));
  data=abs(data);
elseif nRcvrs>1 && Dat.SumOfSquares==2
	% Move individual receiver data to 4th dimension if possible...
	if size(data,4) == 1
		data = reshape(data,size(data,1),size(data,2),size(data,3),size(data,5));
	end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate phase table indexes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pt1,pt2] = l_GetPeTable(hdr)

if ~isfield(hdr,'acqp') || ~isfield(hdr,'method')
	error('Cannot determine phase table, because ACQP and/or METHODS fields are missing.')
end

pt1 = [];
pt2 = [];

if isfield(hdr.acqp,'ACQ_spatial_phase_1')
	[s,pt1] = sort(hdr.acqp.ACQ_spatial_phase_1);
elseif isfield(hdr.method,'PVM_EncSteps1')
	pt1 = hdr.method.PVM_EncSteps1+(-min(hdr.method.PVM_EncSteps1(:))+1);
end

if isfield(hdr.method,'PVM_EncSteps2')
	pt2 = hdr.method.PVM_EncSteps2+(-min(hdr.method.PVM_EncSteps2(:))+1);
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

















