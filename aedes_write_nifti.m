function [done,msg] = aedes_write_nifti(DATA,filename,varargin)
% AEDES_WRITE_NIFTI - Write data in NIfTI or Analyze 7.5 formats
%   
%
% Synopsis: 
%       [done,msg] = aedes_write_nifti(DATA,filename,prop1,value1,prop2,value2,...);
%
% Description:
%       The function writes "DATA" to a file "filename" in NIfTI
%       one file (default), NIfTI two file, or Analyze 7.5
%       format. Function returns a done=1, if the writing was
%       successful. Otherwise, done=0 and the second output argument msg
%       contains the error message.
%
%       DATA can be a valid Aedes DATA-structure or a 2-D, 3-D, or 4-D
%       matrix. Filename is the full path to the output file. If the
%       filename is given without path, the data file is written into the
%       working directory (pwd). If only one input argument is given or
%       the filename is an empty string, the output file name is prompted.
%
%       The possible property/value pairs are the following:
%
%       Property        Value ({ }=default)       Desc.
%       --------        --------                  --------
%
%       'FileType'      [ 0 | 1 | {2} ]           % 0 = Analyze 7.5
%                                                 % 1 = NIfTI (two file)
%                                                 % 2 = NIfTI (one file)
%                                                 % (default)
%
%       'DataType'      [ {0}    | 'uint8' |      % 0 = Determine
%                        'int16' | 'uint16'|      % datatype from data
%                        'int32' | 'uint32'|      % (default)
%                        'int64' | 'uint64'|      % 'str' = specify
%                        'single'| 'double' ]     % datatype
%
%       'VoxelSize'     vox_size                  % A vector
%                                                 % corresponding to the
%                                                 % voxel width in
%                                                 % dimension i:
%                                                 % vox_size(1)= x width
%                                                 % vox_size(2)= y width
%                                                 % vox_size(3)= z width
%                                                 % vox_size(4)= time width
%                                                 
%       'XYZUnits'      [ 'meter' | 'mm' |        % Units of the spatial
%                        'micron' | {'Unknown'} ] % x, y, and z dimensions
%                                               
%       'TimeUnits'     [ 'sec' | 'msec' |        % Units of the temporal
%                        'usec' | 'hz' |          % dimensions
%                         'ppm' | 'rad/sec' |
%                         {'Unknown'} ]
%
%       'RotMtx'        (3x3 or 3x4 matrix)       % User specified rotation
%                                                 % matrix. Stored using
%                                                 % the sform method. NOTE: 
%                                                 % For NIfTI only.  
%
%       'machine'       [ {'ieee-le'}|'ieee-be' ] % little-endian
%                                                 % (default) or
%                                                 % big-endian Byte
%                                                 % ordering
%                       
%       'HeaderOnly'    [ {0} | 1 ]               % 0 = write header and
%                                                 % data (default)
%                                                 % 1 = write only header
%
%       'Clim'          display min/max           % The displayed min and
%                       (1x2 vector, [min,max])   % max intensity values.
%                                                 % The NIfTI viewer may
%                                                 % or may not use these
%                                                 % values. Clim=[0 0] is
%                                                 % the default and means
%                                                 % that the whole range of
%                                                 % the data is shown.
%
%       'description'   (string)                  % A string to be written
%                                                 % in the description
%                                                 % field in NIfTI header.
%
%
% Examples:
%       DATA = aedes_readfid('testi.fid');
%
%       % Write one file NIfTI
%       aedes_write_nifti(DATA,'testi.nii');
%
%       % Write two file NIfTI
%       aedes_write_nifti(DATA,'testi.img','FileType',1)
%
%       % Write Analyze 7.5 *.hdr and *.img files
%       aedes_write_nifti(DATA,'testi.img','FileType',0)
%       
% See also:
%       AEDES, AEDES_READ_NIFTI

%       Acknowledgment: This function is modified under GNU license from
%       MRI_TOOLBOX developed by CNSP in Flinders University,
%       Australia. Some parts are also originally written by Jimmy Shen.
%       See the following links:
%       http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=8797&objectType=file
%       http://eeg.sourceforge.net/
%  
%       NIfTI data format specifications can be found here: 
%       http://nifti.nimh.nih.gov

% This function is a part of Aedes - A graphical tool for analyzing 
% medical images
%
% Copyright (C) 2006 Juha-Pekka Niskanen <Juha-Pekka.Niskanen@uku.fi>
% Copyright (C) 2006 Jimmy Shen
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


%% Default values
Dat.filetype=2;
Dat.machine='ieee-le';
Dat.datatype=0;
Dat.xyzunits = '';
Dat.timeunits = '';
Dat.voxelsize = [];
Dat.Clim=[0 0];
Dat.Descrip = '';
Dat.FlipOrient = 0;
Dat.RotMtx = [];
done=false;
writeOnlyHeader=0;


%% Parse input arguments
switch nargin
 case 0 % ---------------------------------
   done=0;
   msg=['Too few input arguments.'];
   return
 
 case {1, 2} % --------------------------------
   if isnumeric(DATA) || islogical(DATA)
     data = DATA;
     if islogical(data)
       data=uint8(data);
     end
     DATA=[];
     DATA.FTDATA = data;
     DATA.DataFormat = '';
   elseif isstruct(DATA)
     if ~isfield(DATA,'FTDATA')
       error('The field FTDATA not found from the inputted structure!')
     end
     if ~isfield(DATA,'DataFormat')
       DATA.DataFormat = '';
     end
   else
     done=0;
     msg=['First input argument must be structure or numeric type!'];
     return
   end
   
   %% Prompt for file name
   if ~exist('filename','var') || isempty(filename)
     [fname,fpath,findex]=uiputfile({'*.nii;*.NII',...
                         'NIfTI Files - One File Format (*.nii)';...
                   '*.hdr;*.HDR','NIfTI Files - Two File Format (*.hdr)';...
                   '*.hdr;*.HDR','Analyze 7.5 Files (*.hdr)';...
                   '*.*','All Files (*.*)'},'Save Data As');
     if isequal(fname,0) || isequal(fpath,0)
       % Cancelled
       done = false;
       msg = 'Cancelled';
       return
     else
       if findex==1
         Dat.filetype=2;
       elseif findex==2
         Dat.filetype=1;
       elseif findex==3
         Dat.filetype=0;
       else
         Dat.filetype=2;
       end
       filename = [fpath,fname];
     end
   end
  otherwise
    
    if isnumeric(DATA)
      data = DATA;
      DATA=[];
      DATA.FTDATA = data;
      DATA.DataFormat = '';
    elseif isstruct(DATA)
      if ~isfield(DATA,'FTDATA')
        error('The field FTDATA not found from the inputted structure!')
      end
      if ~isfield(DATA,'DataFormat')
        DATA.DataFormat = '';
      end
    else
      done=0;
      msg=['First input argument must be structure or numeric type!'];
      return
    end
    
  %% Parse parameter/value pairs
  for ii=1:2:length(varargin)
    switch lower(varargin{ii})
      case 'machine'
        Dat.machine=varargin{ii+1};
        
      case 'datatype'
        Dat.datatype=varargin{ii+1};
        
      case 'filetype'
        Dat.filetype=varargin{ii+1};
        
      case 'headeronly'
        writeOnlyHeader = varargin{ii+1};
        
      case 'xyzunits'
        Dat.xyzunits = varargin{ii+1};
        
      case 'timeunits'
        Dat.timeunits = varargin{ii+1};
        
      case 'voxelsize'
        Dat.voxelsize = varargin{ii+1};
        
      case 'clim'
        Dat.Clim = varargin{ii+1};
				
			case 'rotmtx'
				tmp_val = varargin{ii+1};
				sz = size(tmp_val);
				if isnumeric(tmp_val) && length(sz)==2 && sz(1)==3 && any(sz(2)==[3 4])
					Dat.RotMtx = varargin{ii+1};
				else
					error('Rotation matrix has to be given as a 3x3 or 3x4 matrix.')
				end
	  case 'description'
		Dat.Descrip = varargin{ii+1};
		
	  case 'FlipInplaneOrient'
		Dat.FlipOrient = varargin{ii+1};
        
     otherwise
      msg=['Unknown parameter "' varargin{ii} '".'];
      return
    end
  end
end

%% Get default datatype from DATA structure if not specified as an input
if ischar(Dat.datatype)
  % Datatype given as a string 'uint16', 'int32', 'single',...
  Dat.datatype=l_GetDataType(Dat.datatype);
elseif Dat.datatype==0
  Dat.datatype=l_GetDataType(DATA);
end

%% Parse filename
[fp,fn,fe]=fileparts(filename);
if ~isempty(fp)
  fp = [fp,filesep];
else
  fp = [pwd,filesep];
end
  
if ~any(strcmpi(fe,{'.nii','.hdr','.img'}))
  fn=[fn,fe];
end


%% Construct header structure
[done,msg,hdr] = l_ConstructValidHeader(DATA,Dat);
if ~done
  return
end

%% Open file(s) for writing
if any(Dat.filetype==[0 1]) % Analyze75 and two file NIfTI
  fid_hdr = fopen([fp,fn,'.hdr'],'w',Dat.machine);
  if fid_hdr<0
	done = false;
    msg={'Could not open file',['"',fp,fn,'.hdr"'],...
         'for writing'};
    return
  end
  fid_img = fopen([fp,fn,'.img'],'w',Dat.machine);
  if fid_hdr<0
	done = false;
    msg={'Could not open file',['"',fp,fn,'.img"'],...
         'for writing'};
    return
  end
  
  %% Write header
  [done,msg] = l_WriteNiftiHdr(hdr,fid_hdr,Dat);
  if ~done
    fclose(fid_hdr);
    return
  end
  fclose(fid_hdr);
  
  %% Write data
  [done,msg] = l_WriteNiftiData(DATA,hdr,fid_img,Dat);
  if ~done
    fclose(fid_img);
    return
  end
  fclose(fid_img);
elseif Dat.filetype==2 % On file NIfTI (default)
  fid = fopen([fp,fn,'.nii'],'w',Dat.machine);
  if fid<0
	done = false;
    msg={'Could not open file',['"',fp,fn,'.nii"'],...
         'for writing'};
    return
  end
  
  %% Write header
  [done,msg] = l_WriteNiftiHdr(hdr,fid,Dat);
  if ~done
    fclose(fid);
    return
  end
  
  %% Write data
  [done,msg] = l_WriteNiftiData(DATA,hdr,fid,Dat);
  if ~done
    fclose(fid);
    return
  end
  fclose(fid);
else
  msg={'Unsupported filetype.',...
       '0=Analyze75 format (*.hdr,*.img)',...
       '1=NIfTI (*.hdr, *.img)',...
       '2=NIfTI (*.nii)'};
end


%% All seems to be well and we can exit normally without errors
done=true;
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct header structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [done,msg,hdr] = l_ConstructValidHeader(DATA,Dat)

hdr = [];
msg='';
done=false;

%% If Analyze75 <-> NIfTI
try
  if isfield(DATA,'DataFormat') && any(strcmpi(DATA.DataFormat,{'Analyze75','NIfTI(1)','NIfTI(2)'}))
    hdr = DATA.HDR.FileHeader;
    
    %% NIfTI -> Analyze 7.5
    if Dat.filetype==0 && any(strcmpi(DATA.DataFormat,{'NIfTI(1)','NIfTI(2)'}))
	  
      hdr.hist.orient      = char(0);
      hdr.hist.originator=zeros(1,5);
      hdr.hist.originator(1)=round(hdr.dime.dim(2)/2);
      hdr.hist.originator(2)=round(hdr.dime.dim(3)/2);
      hdr.hist.originator(3)=round(hdr.dime.dim(4)/2);
      hdr.hist.generated   = ['Aedes',char([0 0 0 0 0])];
      hdr.hist.scannum     = char(zeros(1,10));
      hdr.hist.patient_id  = char(zeros(1,10));
      hdr.hist.exp_date    = char(zeros(1,10));
      hdr.hist.exp_time    = char(zeros(1,10));
      hdr.hist.hist_un0    = char(zeros(1,3));
      hdr.hist.views       = 0;
      hdr.hist.vols_added  = 0;
      hdr.hist.start_field = 0;
      hdr.hist.field_skip  = 0;
      hdr.hist.omax        = 0;
      hdr.hist.omin        = 0;
      hdr.hist.smax        = 0;
      hdr.hist.smin        = 0;
      hdr.hist.magic = '';
      hdr.dime.vox_offset=0;
	  hdr.dime.pixdim(1)=0;
      if any(Dat.Clim~=0)
        hdr.dime.cal_max = Dat.Clim(2);
        hdr.dime.cal_min = Dat.Clim(1);
      end
      
    %% Analyze 7.5 -> NIfTI
    elseif any(Dat.filetype==[1 2]) && strcmpi(DATA.DataFormat,'Analyze75')
      hdr.hist.qform_code=0;
      hdr.hist.sform_code=0;
      hdr.hist.quatern_b=0;
      hdr.hist.quatern_c=0;
      hdr.hist.quatern_d=0;
      hdr.hist.qoffset_x=0;
      hdr.hist.qoffset_y=0;
      hdr.hist.qoffset_z=0;
      hdr.hist.srow_x=zeros(1,4);
      hdr.hist.srow_y=zeros(1,4);
      hdr.hist.srow_z=zeros(1,4);
      hdr.hist.intent_name='';
      hdr.hist.originator=zeros(1,5);
      hdr.hist.originator(1)=round(hdr.dime.dim(2)/2);
      hdr.hist.originator(2)=round(hdr.dime.dim(3)/2);
      hdr.hist.originator(3)=round(hdr.dime.dim(4)/2);
      if Dat.filetype==1
        hdr.hist.magic='ni1';
        hdr.dime.vox_offset=0;
      elseif Dat.filetype==2
        hdr.hist.magic='n+1';
        hdr.dime.vox_offset=352;
        hdr.hist.sform_code = 1;
        hdr.hist.srow_x(1) = hdr.dime.pixdim(2);
        hdr.hist.srow_y(2) = hdr.dime.pixdim(3);
        hdr.hist.srow_z(3) = hdr.dime.pixdim(4);
        hdr.hist.srow_x(4) = (1-hdr.hist.originator(1))*hdr.dime.pixdim(2);
        hdr.hist.srow_y(4) = (1-hdr.hist.originator(2))*hdr.dime.pixdim(3);
        hdr.hist.srow_z(4) = (1-hdr.hist.originator(3))*hdr.dime.pixdim(4);
      end
      if any(Dat.Clim~=0)
        hdr.dime.cal_max = Dat.Clim(2);
        hdr.dime.cal_min = Dat.Clim(1);
      end
      
    %% NIfTI -> NIfTI or Analyze75 -> Analyze75
    else
      if any(Dat.Clim~=0)
        hdr.dime.cal_max = Dat.Clim(2);
        hdr.dime.cal_min = Dat.Clim(1);
	  end
	  if ~isempty(Dat.voxelsize)
		hdr.dime.pixdim(2:length(Dat.voxelsize)+1) = Dat.voxelsize;
	  end
	  if ~isempty(Dat.Descrip)
		hdr.hist.descrip=Dat.Descrip;
	  end
	  data_dim(1)=size(DATA.FTDATA,1);
	  data_dim(2)=size(DATA.FTDATA,2);
	  data_dim(3)=size(DATA.FTDATA,3);
	  data_dim(4)=size(DATA.FTDATA,4);
	  
	  % Swap 1st and 2nd dimesions
	  data_dim = [data_dim(2) data_dim(1) data_dim(3:end)];
	  
	  hdr.dime.dim=[length(size(DATA.FTDATA)) data_dim 1 1 1]; % Data dimensions:
	                                   % dim[0] = nbr of dimensions
	                                   % dim[i] = length of dimension i
	                                   % (i=1..7)
	  
	  % Make sure that datatype is valid
	  % Possible value pairs for DataType and BitPix 
	  d_type = [0 1 2 4 8 16 32 64 128 256 511 512 768 1024 1280 1536 1792 2048];
	  b_pix = [0 1 8 16 32 32 64 64 24 8 96 16 32 64 64 128 128 256];
	  hdr.dime.datatype = Dat.datatype;
	  hdr.dime.bitpix = b_pix(d_type==Dat.datatype);
	  
	end
% $$$       if hdr.hist.qform_code == 0 & hdr.hist.sform_code == 0
% $$$         hdr.hist.sform_code = 1;
% $$$         hdr.hist.srow_x(1) = hdr.dime.pixdim(2);
% $$$         hdr.hist.srow_y(2) = hdr.dime.pixdim(3);
% $$$         hdr.hist.srow_z(3) = hdr.dime.pixdim(4);
% $$$         hdr.hist.srow_x(4) = (1-hdr.hist.originator(1))*hdr.dime.pixdim(2);
% $$$         hdr.hist.srow_y(4) = (1-hdr.hist.originator(2))*hdr.dime.pixdim(3);
% $$$         hdr.hist.srow_z(4) = (1-hdr.hist.originator(3))*hdr.dime.pixdim(4);
% $$$       end
% $$$       
% $$$       %% Set the magic field
% $$$       switch filetype
% $$$        case 1 % NIfTI one file
% $$$         hdr.hist.magic = 'ni1';
% $$$         hdr.dime.vox_offset=0;
% $$$        case 2 % NIfTI two file
% $$$         hdr.hist.magic = 'n+1';
% $$$         hdr.dime.vox_offset=352;
% $$$        otherwise
% $$$         msg=['Unexpected filetype encountered while constructing header.'];
% $$$         return
% $$$       end
% $$$     end
  else
    %% For arbitrary DATA structure formats the header structure has to be
    %  constructed from scratch...
    
    % Possible value pairs for DataType and BitPix 
    d_type = [0 1 2 4 8 16 32 64 128 256 511 512 768 1024 1280 1536 1792 2048]; 
    b_pix = [0 1 8 16 32 32 64 64 24 8 96 16 32 64 64 128 128 256];
    
    % Header key
    hdr.hk.sizeof_hdr=348;
    hdr.hk.data_type='';
    hdr.hk.db_name='';
    hdr.hk.extents=0;
    hdr.hk.session_error=0;
    hdr.hk.regular='r';
    hdr.hk.dim_info=0;
    
    % Image dimensions
    data_dim(1)=size(DATA.FTDATA,1);
	data_dim(2)=size(DATA.FTDATA,2);
	data_dim(3)=size(DATA.FTDATA,3);
	data_dim(4)=size(DATA.FTDATA,4);
	
	% Swap 1st and 2nd dimesions
	data_dim = [data_dim(2) data_dim(1) data_dim(3:end)];
	
    hdr.dime.dim=[length(size(DATA.FTDATA)) data_dim 1 1 1]; % Data dimensions:
                                     % dim[0] = nbr of dimensions
                                     % dim[i] = length of dimension i (i=1..7)
    hdr.dime.intent_p1=0;
    hdr.dime.intent_p2=0;
    hdr.dime.intent_p3=0;
    hdr.dime.intent_code=0;
    hdr.dime.datatype = Dat.datatype;
    hdr.dime.bitpix = b_pix(d_type==Dat.datatype);
    hdr.dime.slice_start = 0;
    hdr.dime.pixdim = zeros(1,8); % pixdim[i] = voxel width along
                                  % dimension i
    if ~isempty(Dat.voxelsize)
      %hdr.dime.pixdim(1)=length(Dat.voxelsize);
      hdr.dime.pixdim(2:length(Dat.voxelsize)+1) = Dat.voxelsize;
    else
      hdr.dime.pixdim(2:4)=1;
    end
    hdr.dime.vox_offset=0;
    hdr.dime.scl_slope=0;
    hdr.dime.scl_inter=0;
    hdr.dime.slice_end=0;
    hdr.dime.slice_code=0;
    

	xyz_bits = '000'; % Unknown
	time_bits = '000'; % Unknown
    hdr.dime.xyzt_units=bin2dec(['00',time_bits,xyz_bits]);
    
    hdr.dime.cal_max=Dat.Clim(2);
    hdr.dime.cal_min=Dat.Clim(1);
    hdr.dime.slice_duration=0;
    hdr.dime.toffset=0;
    hdr.dime.glmax=round(max(DATA.FTDATA(:)));
    hdr.dime.glmin=round(min(DATA.FTDATA(:)));
    
    % Data history
	if ~isempty(Dat.Descrip)
	  hdr.hist.descrip=Dat.Descrip;
	else
	  hdr.hist.descrip=['Converted from ' DATA.DataFormat];
	end
    hdr.hist.aux_file='none';
    
    
    %% Analyze 7.5 style data history
    if Dat.filetype==0
      hdr.hist.orient      = char(0);
      hdr.hist.originator=zeros(1,5);
      hdr.hist.originator(1)=round(hdr.dime.dim(2)/2);
      hdr.hist.originator(2)=round(hdr.dime.dim(3)/2);
      hdr.hist.originator(3)=round(hdr.dime.dim(4)/2);
      hdr.hist.generated   = ['Aedes',char([0 0])];
      hdr.hist.scannum     = char(zeros(1,10));
      hdr.hist.patient_id  = char(zeros(1,10));
      hdr.hist.exp_date    = char(zeros(1,10));
      hdr.hist.exp_time    = char(zeros(1,10));
      hdr.hist.hist_un0    = char(zeros(1,3));
      hdr.hist.views       = 0;
      hdr.hist.vols_added  = 0;
      hdr.hist.start_field = 0;
      hdr.hist.field_skip  = 0;
      hdr.hist.omax        = 0;
      hdr.hist.omin        = 0;
      hdr.hist.smax        = 0;
      hdr.hist.smin        = 0;
      
      % Set the magic field (NIfTI identifier)
      hdr.hist.magic       = '';
    elseif any(Dat.filetype==[1 2])
      hdr.hist.qform_code=0;
      hdr.hist.sform_code=0;
      hdr.hist.quatern_b=0;
      hdr.hist.quatern_c=0;
      hdr.hist.quatern_d=0;
      hdr.hist.qoffset_x=0;
      hdr.hist.qoffset_y=0;
      hdr.hist.qoffset_z=0;
      hdr.hist.srow_x=zeros(1,4);
      hdr.hist.srow_y=zeros(1,4);
      hdr.hist.srow_z=zeros(1,4);
      hdr.hist.intent_name='';
      hdr.hist.originator=zeros(1,5);
      hdr.hist.originator(1)=round(hdr.dime.dim(2)/2);
      hdr.hist.originator(2)=round(hdr.dime.dim(3)/2);
      hdr.hist.originator(3)=round(hdr.dime.dim(4)/2);
      if Dat.filetype==1
        hdr.hist.magic='ni1';
      elseif Dat.filetype==2
        hdr.hist.magic='n+1';
        hdr.dime.vox_offset=352;
        hdr.hist.sform_code = 1;
        hdr.hist.srow_x(1) = hdr.dime.pixdim(2);
        hdr.hist.srow_y(2) = hdr.dime.pixdim(3);
        hdr.hist.srow_z(3) = hdr.dime.pixdim(4);
        hdr.hist.srow_x(4) = (1-hdr.hist.originator(1))*hdr.dime.pixdim(2);
        hdr.hist.srow_y(4) = (1-hdr.hist.originator(2))*hdr.dime.pixdim(3);
        hdr.hist.srow_z(4) = (1-hdr.hist.originator(3))*hdr.dime.pixdim(4);
      end
    end
	end
	
	% Input a custom rotation matrix to NIfTI file, if given
	if any(Dat.filetype==[1 2]) && ~isempty(Dat.RotMtx)
		hdr.hist.sform_code=1;
		hdr.hist.srow_x = Dat.RotMtx(1,1:3);
		hdr.hist.srow_y = Dat.RotMtx(2,1:3);
		hdr.hist.srow_z = Dat.RotMtx(3,1:3);
		if size(Dat.RotMtx,2)==3
			hdr.hist.srow_x(4) = (1-hdr.hist.originator(1))*hdr.dime.pixdim(2);
			hdr.hist.srow_y(4) = (1-hdr.hist.originator(2))*hdr.dime.pixdim(3);
			hdr.hist.srow_z(4) = (1-hdr.hist.originator(3))*hdr.dime.pixdim(4);
		else
			hdr.hist.srow_x(4) = Dat.RotMtx(1,4);
			hdr.hist.srow_y(4) = Dat.RotMtx(2,4);
			hdr.hist.srow_z(4) = Dat.RotMtx(3,4);
		end
	end
  
  % Units of spatial and temporal dimensions:
  % 0 = Unknown, 1 = meter, 2 = millimeter, 3 = micrometer, 8 = seconds
  % 16 = milliseconds, 24 = microseconds, 32 = Hertz, 40 = ppm,
  % 48 = rad/sec.
  % Bits 0..2 of xyzt_units specify the units of pixdim[1..3]. Bits
  % 3..5 of xyzt_units specify the units of pixdim[4] and are multiples
  % of 8.
  if not(isempty(Dat.xyzunits))
	if not(isfield(DATA,'DataFormat') && strcmpi(DATA.DataFormat,'Analyze75'))
	  if isempty(Dat.xyzunits) || strcmpi(Dat.xyzunits,'unknown')
		xyz_bits = '000';
	  elseif strcmpi(Dat.xyzunits,'meter')
		xyz_bits = dec2bin(1,3);
	  elseif strcmpi(Dat.xyzunits,'mm')
		xyz_bits = dec2bin(2,3);
	  elseif strcmpi(Dat.xyzunits,'micron')
		xyz_bits = dec2bin(3,3);
	  end
	  if isempty(Dat.timeunits) || strcmpi(Dat.timeunits,'unknown')
		time_bits = '000';
	  elseif strcmpi(Dat.timeunits,'sec')
		time_bits = dec2bin(1,3);
	  elseif strcmpi(Dat.timeunits,'msec')
		time_bits = dec2bin(2,3);
	  elseif strcmpi(Dat.timeunits,'usec')
		time_bits = dec2bin(3,3);
	  elseif strcmpi(Dat.timeunits,'hz')
		time_bits = dec2bin(4,3);
	  elseif strcmpi(Dat.timeunits,'ppm')
		time_bits = dec2bin(5,3);
	  elseif strcmpi(Dat.timeunits,'rad/sec')
		time_bits = dec2bin(6,3);
	  end
	  hdr.dime.xyzt_units=bin2dec(['00',time_bits,xyz_bits]);
	end
  end
  done=true;
catch
  msg=lasterr;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define datatype
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [datatype,bitpix]=l_GetDataType(DATA)

%  Set bitpix according to datatype
%
%  /*Acceptable values for datatype are*/ 
%
%     0 None                   (Unknown bit per voxel)   % DT_NONE, DT_UNKNOWN 
%     1 Binary                       (ubit1, bitpix=1)   % DT_BINARY 
%     2 Unsigned char       (uchar or uint8, bitpix=8)   % DT_UINT8, NIFTI_TYPE_UINT8 
%     4 Signed short                 (int16, bitpix=16)  % DT_INT16, NIFTI_TYPE_INT16 
%     8 Signed integer               (int32, bitpix=32)  % DT_INT32, NIFTI_TYPE_INT32 
%    16 Floating point   (single or float32, bitpix=32)  % DT_FLOAT32, NIFTI_TYPE_FLOAT32 
%    32 Complex, 2 float32     (Unsupported, bitpix=64)  % DT_COMPLEX64, NIFTI_TYPE_COMPLEX64
%    64 Double precision (double or float64, bitpix=64)  % DT_FLOAT64, NIFTI_TYPE_FLOAT64 
%   128 uint8 RGB                (Use uint8, bitpix=24)  % DT_RGB24, NIFTI_TYPE_RGB24 
%   256 Signed char          (schar or int8, bitpix=8)   % DT_INT8, NIFTI_TYPE_INT8 
%   511 Single RGB             (Use float32, bitpix=96)  % DT_RGB96, NIFTI_TYPE_RGB96
%   512 Unsigned short              (uint16, bitpix=16)  % DT_UNINT16, NIFTI_TYPE_UNINT16 
%   768 Unsigned integer            (uint32, bitpix=32)  % DT_UNINT32, NIFTI_TYPE_UNINT32 
%  1024 Signed long long             (int64, bitpix=64)  % DT_INT64, NIFTI_TYPE_INT64
%  1280 Unsigned long long          (uint64, bitpix=64)  % DT_UINT64, NIFTI_TYPE_UINT64 
%  1536 Long double, float128  (Unsupported, bitpix=128) % DT_FLOAT128, NIFTI_TYPE_FLOAT128 
%  1792 Complex128, 2 float64  (Unsupported, bitpix=128) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX12 
%  2048 Complex256, 2 float128 (Unsupported, bitpix=256) % DT_COMPLEX128,
%  NIFTI_TYPE_COMPLEX128
if ischar(DATA)
  action = DATA;
else
  action = class(DATA.FTDATA);
end

switch action
 case {'uint8','logical'}
  datatype = 2;
  bitpix = 8;
 case 'int16'
  datatype = 4;
  bitpix = 16;
 case 'uint16'
  datatype = 512;
  bitpix = 16;
 case 'int32'
  datatype = 8;
  bitpix = 32;
 case 'uint32'
  datatype = 768;
  bitpix = 32;
 case 'int64'
  datatype = 1024;
  bitpix = 64;
 case 'uint64'
  datatype = 1280;
  bitpix = 64;
 case {'single','float'}
  datatype = 16;
  bitpix = 32;
 case 'double'
  datatype = 64;
  bitpix = 64;
 otherwise
  % Default to float
  datatype = 16;
  bitpix = 32;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write NIfTI/Analyze75 header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [done,msg] = l_WriteNiftiHdr(hdr,fid,Dat)

done=false;
msg='';

%  Original header structures
%  struct dsr				/* dsr = hdr */
%       { 
%       struct header_key hk;            /*   0 +  40       */
%       struct image_dimension dime;     /*  40 + 108       */
%       struct data_history hist;        /* 148 + 200       */
%       };                               /* total= 348 bytes*/
try
  header_key(fid, hdr.hk);
  image_dimension(fid, hdr.dime);
  data_history(fid, hdr.hist,Dat.filetype);
catch
  msg=lasterr;
  return
end

%  check the file size is 348 bytes
%
fbytes = ftell(fid);
  
if ~isequal(fbytes,348),
  msg = sprintf('Header size is not 348 bytes.');
  return
end

done=true;
return;					% write_header


%---------------------------------------------------------------------
function header_key(fid, hk)

fseek(fid,0,'bof');

%  Original header structures    
%  struct header_key                      /* header key      */ 
%       {                                /* off + size      */
%       int sizeof_hdr                   /*  0 +  4         */
%       char data_type[10];              /*  4 + 10         */
%       char db_name[18];                /* 14 + 18         */
%       int extents;                     /* 32 +  4         */
%       short int session_error;         /* 36 +  2         */
%       char regular;                    /* 38 +  1         */
%       char dim_info;   % char hkey_un0;        /* 39 +  1 */
%       };                               /* total=40 bytes  */

fwrite(fid, hk.sizeof_hdr(1),    'int32');	% must be 348.

% data_type = sprintf('%-10s',hk.data_type);	% ensure it is 10 chars from left
% fwrite(fid, data_type(1:10), 'uchar');
pad = zeros(1, 10-length(hk.data_type));
hk.data_type = [hk.data_type  char(pad)];
fwrite(fid, hk.data_type(1:10), 'uchar');

% db_name   = sprintf('%-18s', hk.db_name);	% ensure it is 18 chars from left
% fwrite(fid, db_name(1:18), 'uchar');
pad = zeros(1, 18-length(hk.db_name));
hk.db_name = [hk.db_name  char(pad)];
fwrite(fid, hk.db_name(1:18), 'uchar');

fwrite(fid, hk.extents(1),       'int32');
fwrite(fid, hk.session_error(1), 'int16');
fwrite(fid, hk.regular(1),       'uchar');	% might be uint8

% fwrite(fid, hk.hkey_un0(1),    'uchar');
% fwrite(fid, hk.hkey_un0(1),    'uint8');
fwrite(fid, hk.dim_info(1),      'uchar');
return;					% header_key


%---------------------------------------------------------------------
function image_dimension(fid, dime)

%  Original header structures        
%  struct image_dimension
%       {                                /* off + size      */
%       short int dim[8];                /* 0 + 16          */
%       float intent_p1;   % char vox_units[4];   /* 16 + 4       */
%       float intent_p2;   % char cal_units[8];   /* 20 + 4       */
%       float intent_p3;   % char cal_units[8];   /* 24 + 4       */
%       short int intent_code;   % short int unused1;   /* 28 + 2 */
%       short int datatype;              /* 30 + 2          */
%       short int bitpix;                /* 32 + 2          */
%       short int slice_start;   % short int dim_un0;   /* 34 + 2 */
%       float pixdim[8];                 /* 36 + 32         */
%			/*
%				pixdim[] specifies the voxel dimensions:
%				pixdim[1] - voxel width
%				pixdim[2] - voxel height
%				pixdim[3] - interslice distance
%				pixdim[4] - volume timing, in msec
%					..etc
%			*/
%       float vox_offset;                /* 68 + 4          */
%       float scl_slope;   % float roi_scale;     /* 72 + 4 */
%       float scl_inter;   % float funused1;      /* 76 + 4 */
%       short slice_end;   % float funused2;      /* 80 + 2 */
%       char slice_code;   % float funused2;      /* 82 + 1 */
%       char xyzt_units;   % float funused2;      /* 83 + 1 */
%       float cal_max;                   /* 84 + 4          */
%       float cal_min;                   /* 88 + 4          */
%       float slice_duration;   % int compressed; /* 92 + 4 */
%       float toffset;   % int verified;          /* 96 + 4 */
%       int glmax;                       /* 100 + 4         */
%       int glmin;                       /* 104 + 4         */
%       };                               /* total=108 bytes */

fwrite(fid, dime.dim(1:8),        'int16');
fwrite(fid, dime.intent_p1(1),  'float32');
fwrite(fid, dime.intent_p2(1),  'float32');
fwrite(fid, dime.intent_p3(1),  'float32');
fwrite(fid, dime.intent_code(1),  'int16');
fwrite(fid, dime.datatype(1),     'int16');
fwrite(fid, dime.bitpix(1),       'int16');
fwrite(fid, dime.slice_start(1),  'int16');
fwrite(fid, dime.pixdim(1:8),   'float32');
fwrite(fid, dime.vox_offset(1), 'float32');
fwrite(fid, dime.scl_slope(1),  'float32');
fwrite(fid, dime.scl_inter(1),  'float32');
fwrite(fid, dime.slice_end(1),    'int16');
fwrite(fid, dime.slice_code(1),   'uchar');
fwrite(fid, dime.xyzt_units(1),   'uchar');
fwrite(fid, dime.cal_max(1),    'float32');
fwrite(fid, dime.cal_min(1),    'float32');
fwrite(fid, dime.slice_duration(1), 'float32');
fwrite(fid, dime.toffset(1),    'float32');
fwrite(fid, dime.glmax(1),        'int32');
fwrite(fid, dime.glmin(1),        'int32');
return;					% image_dimension


%---------------------------------------------------------------------
function data_history(fid, hist, filetype)

% Original header structures
%struct data_history       
%       {                                /* off + size      */
%       char descrip[80];                /* 0 + 80          */
%       char aux_file[24];               /* 80 + 24         */
%       short int qform_code;            /* 104 + 2         */
%       short int sform_code;            /* 106 + 2         */
%       float quatern_b;                 /* 108 + 4         */
%       float quatern_c;                 /* 112 + 4         */
%       float quatern_d;                 /* 116 + 4         */
%       float qoffset_x;                 /* 120 + 4         */
%       float qoffset_y;                 /* 124 + 4         */
%       float qoffset_z;                 /* 128 + 4         */
%       float srow_x[4];                 /* 132 + 16        */
%       float srow_y[4];                 /* 148 + 16        */
%       float srow_z[4];                 /* 164 + 16        */
%       char intent_name[16];            /* 180 + 16        */
%       char magic[4];   % int smin;     /* 196 + 4         */
%       };                               /* total=200 bytes */

% descrip     = sprintf('%-80s', hist.descrip);     % 80 chars from left
% fwrite(fid, descrip(1:80),    'uchar');
pad = zeros(1, 80-length(hist.descrip));
hist.descrip = [hist.descrip  char(pad)];
fwrite(fid, hist.descrip(1:80), 'uchar');

% aux_file    = sprintf('%-24s', hist.aux_file);    % 24 chars from left
% fwrite(fid, aux_file(1:24),   'uchar');
pad = zeros(1, 24-length(hist.aux_file));
hist.aux_file = [hist.aux_file  char(pad)];
fwrite(fid, hist.aux_file(1:24), 'uchar');

% Write NIfTI style data history
if any(filetype==[1 2])
  fwrite(fid, hist.qform_code,    'int16');
  fwrite(fid, hist.sform_code,    'int16');
  fwrite(fid, hist.quatern_b,   'float32');
  fwrite(fid, hist.quatern_c,   'float32');
  fwrite(fid, hist.quatern_d,   'float32');
  fwrite(fid, hist.qoffset_x,   'float32');
  fwrite(fid, hist.qoffset_y,   'float32');
  fwrite(fid, hist.qoffset_z,   'float32');
  fwrite(fid, hist.srow_x(1:4), 'float32');
  fwrite(fid, hist.srow_y(1:4), 'float32');
  fwrite(fid, hist.srow_z(1:4), 'float32');

  % intent_name = sprintf('%-16s', hist.intent_name);	% 16 chars from left
  % fwrite(fid, intent_name(1:16),    'uchar');
  pad = zeros(1, 16-length(hist.intent_name));
  hist.intent_name = [hist.intent_name  char(pad)];
  fwrite(fid, hist.intent_name(1:16), 'uchar');
  
  % magic	= sprintf('%-4s', hist.magic);		% 4 chars from left
  % fwrite(fid, magic(1:4),           'uchar');
  pad = zeros(1, 4-length(hist.magic));
  hist.magic = [hist.magic  char(pad)];
  fwrite(fid, hist.magic(1:4),        'uchar');
elseif filetype==0 % Write Analyze 7.5 style (old) data history
  fwrite(fid,hist.orient,'uchar');
  fwrite(fid,hist.originator,'int16');
  fwrite(fid,hist.generated,'uchar');
  fwrite(fid,hist.scannum,'uchar');
  fwrite(fid,hist.patient_id,'uchar');
  fwrite(fid,hist.exp_date,'uchar');
  fwrite(fid,hist.exp_time,'uchar');
  fwrite(fid,hist.hist_un0,'uchar');
  fwrite(fid,hist.views,'int32');
  fwrite(fid,hist.vols_added,'int32');
  fwrite(fid,hist.start_field,'int32');
  fwrite(fid,hist.field_skip,'int32');
  fwrite(fid,hist.omax,'int32');
  fwrite(fid,hist.omin,'int32');
  fwrite(fid,hist.smax,'int32');
  fwrite(fid,hist.smin,'int32');
  
  % Write also the magic field NIfTI/Analyze75 identifier
  pad=zeros(1,4-length(hist.magic));
  hist.magic=[hist.magic char(pad)];
  fseek(fid,344,'bof');
  fwrite(fid,hist.magic(1:4),'uchar');
end
return;					% data_history


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write NIfTI/Analyze75 file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [done,msg] = l_WriteNiftiData(DATA,hdr,fid,Dat)

done=false;
msg='';

switch double(hdr.dime.datatype),
 case   1,
  precision = 'ubit1';
 case   2,
  precision = 'uint8';
 case   4,
  precision = 'int16';
 case   8,
  precision = 'int32';
 case  16,
  precision = 'float32';
 case  64,
  precision = 'float64';
 case 128,
  precision = 'uint8';
 case 256 
  precision = 'int8';
 case 512 
  precision = 'uint16';
 case 768 
  precision = 'uint32';
 case 1024
  precision = 'int64';
 case 1280
  precision = 'uint64';
 otherwise
  msg='Unsupported datatype.';
  return
end

%% if original data is readed from an Analyze75 or NIfTI file, it is very
%% likely that the data is already calibrated using scl_slope and
%% scl_interp. Because of this, the data is "decalibrated".
% ---- On second thought, don't "decalibrate"
% if isfield(DATA,'DataFormat') && ...
%     any(strcmpi(DATA.DataFormat,{'analyze75','nifti(1)','nifti(2)'}))
%   if DATA.HDR.FileHeader.dime.scl_slope~=0
%     DATA.FTDATA = (DATA.FTDATA-DATA.HDR.FileHeader.dime.scl_inter)/...
%         DATA.HDR.FileHeader.dime.scl_slope;
%   end
% end


try
  ScanDim = double(hdr.dime.dim(5));		% t
  SliceDim = double(hdr.dime.dim(4));		% z
  RowDim   = double(hdr.dime.dim(3));		% y
  PixelDim = double(hdr.dime.dim(2));		% x
  SliceSz  = double(hdr.dime.pixdim(4));
  RowSz    = double(hdr.dime.pixdim(3));
  PixelSz  = double(hdr.dime.pixdim(2));
  
  x = 1:PixelDim;
  
  if Dat.filetype == 2
    skip_bytes = double(hdr.dime.vox_offset) - 348;
  else
    skip_bytes = 0;
  end
  
  if double(hdr.dime.datatype) == 128
    DATA.FTDATA = permute(DATA.FTDATA, [4 1 2 3 5]);
  end
  
  if skip_bytes
    fwrite(fid, ones(1,skip_bytes), 'uint8');
  end
  %flipdim(permute(img,[2 1 3 4]),1)
  fwrite(fid,permute(flipdim(DATA.FTDATA,1),[2 1 3 4]), precision);
  %   fwrite(fid, nii.img, precision, skip_bytes);        % error using skip
  %fclose(fid);
catch
  msg=lasterr;
end

done=true;
return;					% write_nii
