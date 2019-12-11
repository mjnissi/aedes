function [DATA,msg] = aedes_read_nifti(filename,varargin)
% AEDES_READ_NIFTI - Read NIfTI (*.nii) and Analyze 7.5 (*.hdr,*.img) files
%   
%
% Synopsis: 
%       [DATA,msg]=aedes_read_nifti(filename,opt)
%       [DATA,msg]=aedes_read_nifti(DATA.HDR)
%
% Description:
%       The function reads the NIfTI and Analyze 7.5 file formats and
%       returns data in DATA.FTDATA and header in DATA.HDR (for more
%       information about the structure of DATA, see the help of AEDES_READFID
%       function, i.e. type "help aedes_readfid"). If error has occurred during
%       reading, the DATA structure is empty and the second output
%       argument msg contains the error message. The function is capable
%       of reading NIfTI and Analyze 7.5 files in the two file
%       format (*.img and *.hdr) and the NIfTI on file format (*.nii).
%
%       First input argument can be either the full path to the NIfTI or
%       Analyze 7.5 file or the NIfTI/Analyze7.5 header structure. If the
%       second input argument is a string 'header', only the header of
%       the file is read. If the function is called without input
%       arguments, a file dialog will appear to browse for a file.
%
%
%       Acknowledgment: This function is modified under GNU license from
%       MRI_TOOLBOX developed by CNSP in Flinders University,
%       Australia. Some parts are also originally written by Jimmy Shen.
%       See the following links:
%       http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=8797&objectType=file
%       http://eeg.sourceforge.net/
%  
%       NIfTI data format specifications can be found here: 
%       http://nifti.nimh.nih.gov
%
% Examples:
%       [DATA,msg]=aedes_read_nifti(filename,'header'); % Read only header
%       [DATA,msg]=aedes_read_nifti(DATA.HDR);          % Read data
%       or
%       [DATA,msg]=aedes_read_nifti;                    % Browse for a file and
%                                                 % read both header and
%                                                 % data
%       or
%       [DATA,msg]=aedes_read_nifti(filename)           % Read both header and
%                                                 % data
%
% See also:
%       AEDES_READFID, AEDES, AEDES_WRITE_NIFTI

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


DATA=[];
data=[];
msg='';

ReadHdr = true;
ReadData = true;
machine=[];

%% Parse input arguments
switch nargin
 
  %% Ask for a file name, if not provided in the input argument
 case 0,
  [f_name, f_path, f_index] = uigetfile( ...
       {'*.nii;*.NII;*.nii.gz;*.hdr;*.HDR',...
        'NIfTI and Analyze 7.5 files (*.nii,*.nii.gz,*.hdr)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Select NIfTI or Analyze 7.5 files');
  if ( all(f_name==0) | all(f_path==0) ) % Cancel is pressed
    DATA=[];
    msg_out='Canceled';
    return
  end
  filename=[f_path,f_name];
  
 case 1,
  if isstruct(filename)
    hdr=filename;
    filename=[hdr.fpath,hdr.fname];
	machine = hdr.byteorder;
    ReadHdr = false;
    ReadData = true;
  else
    ReadHdr = true;
    ReadData = true;
  end
 case 2,
  % Read only header
  if strcmpi(varargin{1},'header')
    ReadHdr = true;
    ReadData = false;
  else
    msg=['Unknown option "' varargin{1} '".'];
    return
  end
 otherwise
  msg='Too many input arguments';
  return
end

%% Parse filename
[f_path,f_name,f_ext] = fileparts(filename);
if isempty(f_path)
  % If there is no path information in "filename", do some extra stuff
  % here...
  filename = which(filename);
  [f_path,f_name,f_ext] = fileparts(filename);
end

%% Check if the file is gzipped NIfTI (*.nii.gz)
gz_filename = '';
if strcmpi(f_ext,'.gz') && length(f_name)>3 && ...
	strcmpi(f_name(end-3:end),'.nii')
  
  % Try to gunzip the file to a temporary folder and read it from there
  try
	gunzip(filename,tempdir);
  catch
	error('Could not gunzip "%s" to "%s"',filename,tempdir)
  end
  gz_filename = fullfile(tempdir,f_name);
  [f_path,f_name,f_ext]=fileparts(gz_filename);
end

%% Read header
if ReadHdr
  [hdr,tmp_msg,machine]=l_ReadHdr([f_path,filesep,f_name],f_ext,machine);
  if isempty(hdr)
    DATA=[];
    if iscell(tmp_msg)
      msg = {'Error while reading data header.',...
             tmp_msg{:}};
    else
      msg = {'Error while reading data header.',...
             tmp_msg};
    end
    return
  end
end

%% Set dataformat string
if strcmp(hdr.FileHeader.hist.magic,'n+1')
  DATA.DataFormat = 'NIfTI(1)'; % NIfTI format (*.nii), one file
  filetype=2;
elseif strcmp(hdr.FileHeader.hist.magic,'ni1')
  DATA.DataFormat = 'NIfTI(2)'; % NIfTI format (*.hdr,*.img), two files
  filetype=1;
else
  DATA.DataFormat = 'Analyze75'; % Analyze 7.5 format (*.hdr,*.img), two
                                % files
  filetype=0;
end

%% Read data
if ReadData
  [data,hdr,tmp_msg]=l_ReadData(hdr,machine,[f_path,filesep,f_name,f_ext],filetype);
  if isempty(data)
    DATA=[];
    if iscell(tmp_msg)
      msg = {'Error while reading data.',...
             tmp_msg{:}};
    else
      msg = {'Error while reading data.',...
             tmp_msg};
    end
    return
  end
end

%% Set fields to DATA
if isempty(gz_filename)
  DATA.HDR.FileHeader = hdr.FileHeader;
  DATA.HDR.fname = [f_name,f_ext];
  DATA.HDR.fpath = [f_path,filesep];
  DATA.HDR.byteorder = machine;
else
  % In .nii.gz files, use the original filename
  [f_path,f_name,f_ext] = fileparts(filename);
  
  % Clean up -> remove the temporary gunzipped file
  try
	delete(gz_filename)
  catch
	warning(['Could not remove temporary file "%s". ',...
	  'Maybe you should see what''s wrong...'],gz_filename)
  end
  
  DATA.HDR.FileHeader = hdr.FileHeader;
  DATA.HDR.fname = [f_name,f_ext];
  DATA.HDR.fpath = [f_path,filesep];
  DATA.HDR.byteorder = machine;
end


% Set xyz and time units
tmp=dec2bin(DATA.HDR.FileHeader.dime.xyzt_units,8);

xyzunits = tmp(6:8);
timeunits = tmp(3:5);

% Units of spatial and temporal dimensions:
% 0 = Unknown, 1 = meter, 2 = millimeter, 3 = micrometer, 8 = seconds
% 16 = milliseconds, 24 = microseconds, 32 = Hertz, 40 = ppm,
% 48 = rad/sec.
% Bits 0..2 of xyzt_units specify the units of pixdim[1..3]. Bits
% 3..5 of xyzt_units specify the units of pixdim[4] and are multiples
% of 8.	  
if bin2dec(xyzunits)==0
  DATA.HDR.xyzunits = 'Unknown';
elseif bin2dec(xyzunits)==1
  DATA.HDR.xyzunits = 'meter';
elseif bin2dec(xyzunits)==2
  DATA.HDR.xyzunits = 'mm';
elseif bin2dec(xyzunits)==3
  DATA.HDR.xyzunits = 'micron';
end
if bin2dec(timeunits)==0
  DATA.HDR.timeunits = 'Unknown';
elseif bin2dec(timeunits)==1
  DATA.HDR.timeunits = 'sec';
elseif bin2dec(timeunits)==2
  DATA.HDR.timeunits = 'msec';
elseif bin2dec(timeunits)==3
  DATA.HDR.timeunits = 'usec';
elseif bin2dec(timeunits)==4
  DATA.HDR.timeunits = 'Hz';
elseif bin2dec(timeunits)==5
  DATA.HDR.timeunits = 'ppm';
elseif bin2dec(timeunits)==6
  DATA.HDR.timeunits = 'rad/sec';
end
if ReadData
  DATA.FTDATA = data;
  DATA.KSPACE=[];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Read NIfTI/Analyze75 header
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hdr,msg,machine]=l_ReadHdr(filename,ext,machine)

hdr=[];
msg='';

%% Try first little-endian byte ordering
if isempty(machine)
  machine = 'ieee-le';
end

new_ext = 0;

if strcmpi(ext,'.nii')
  new_ext = 1;
else
  new_ext = 0;
end

if any(strcmpi(ext,{'.hdr','.nii'}))
  fname = [filename,ext];
elseif strcmpi(ext,'.img')
  fname = [filename,'.hdr'];
end


fid = fopen(fname,'r',machine);
if fid < 0,
  hdr=[];
  msg = sprintf('Cannot open file %s.',fname);
  return
end

%% Try to determine if the header file is Analyze75 or NIfTI
fseek(fid,344,'bof');
tmp=fread(fid,4,'*char');
tmp=deblank(tmp');
if isempty(tmp)
  filetype=0; % Analyze 7.5
elseif strcmpi(tmp,'ni1')
  filetype=1; % Two file NIfTI
elseif strcmpi(tmp,'n+1')
  filetype=2; % One file NIfTI
else
  % Somethings wrong -> throw an error
  fclose(fid);
  hdr=[];
  msg = sprintf(['The file "%s" is not a valid NIfTI/Analyze 7.5 header ' ...
                 'file.'],fname);
  return
end

% Rewind file
fseek(fid,0,'bof');

% Read header
hdr.FileHeader = read_header(fid,filetype);
fclose(fid);





% Check if machine format was correct
if hdr.FileHeader.hk.sizeof_hdr ~= 348
  % first try reading the opposite endian to 'machine'
  switch machine,
   case 'ieee-le', machine = 'ieee-be';
   case 'ieee-be', machine = 'ieee-le';
  end
  
  fid = fopen(fname,'r',machine);
  if fid < 0,
    hdr=[];
    msg = sprintf('Cannot open file %s.',fname);
  else
    hdr.FileHeader = read_header(fid,filetype);
    fclose(fid);
  end
end


if hdr.FileHeader.hk.sizeof_hdr ~= 348
  % Now throw an error
  hdr=[];
  msg = sprintf('File "%s" is probably corrupted.',fname);
end

% $$$ if strcmp(hdr.hist.magic, 'n+1')
% $$$   filetype = 2;
% $$$ elseif strcmp(hdr.hist.magic, 'ni1')
% $$$   filetype = 1;
% $$$ else
% $$$   filetype = 0;
% $$$ end

return					% load_nii_hdr


%---------------------------------------------------------------------
function [ dsr ] = read_header(fid,filetype)

%  Original header structures
%  struct dsr
%       { 
%       struct header_key hk;            /*   0 +  40       */
%       struct image_dimension dime;     /*  40 + 108       */
%       struct data_history hist;        /* 148 + 200       */
%       };                               /* total= 348 bytes*/

dsr.hk   = header_key(fid);
dsr.dime = image_dimension(fid);
dsr.hist = data_history(fid,filetype);

%  For Analyze data format
%
% $$$ if ~strcmp(dsr.hist.magic, 'n+1') & ~strcmp(dsr.hist.magic, 'ni1')
% $$$   dsr.dime.scl_slope = 0;
% $$$   dsr.hist.qform_code = 0;
% $$$   dsr.hist.sform_code = 0;
% $$$ end

return					% read_header


%---------------------------------------------------------------------
function [ hk ] = header_key(fid)

fseek(fid,0,'bof');

%  Original header structures	
%  struct header_key                     /* header key      */ 
%       {                                /* off + size      */
%       int sizeof_hdr                   /*  0 +  4         */
%       char data_type[10];              /*  4 + 10         */
%       char db_name[18];                /* 14 + 18         */
%       int extents;                     /* 32 +  4         */
%       short int session_error;         /* 36 +  2         */
%       char regular;                    /* 38 +  1         */
%       char dim_info;   % char hkey_un0;        /* 39 +  1 */
%       };                               /* total=40 bytes  */
%
% int sizeof_header   Should be 348.
% char regular        Must be 'r' to indicate that all images and 
%                     volumes are the same size. 

hk.sizeof_hdr    = fread(fid, 1,'int32')';	% should be 348!
hk.data_type     = deblank(fread(fid,10,'*char')');
hk.db_name       = deblank(fread(fid,18,'*char')');
hk.extents       = fread(fid, 1,'int32')';
hk.session_error = fread(fid, 1,'int16')';
hk.regular       = fread(fid, 1,'*char')';
hk.dim_info      = fread(fid, 1,'char')';

return					% header_key


%---------------------------------------------------------------------
function [ dime ] = image_dimension(fid)

%  Original header structures    
%  struct image_dimension
%       {                                /* off + size      */
%       short int dim[8];                /* 0 + 16          */
%       /*
%           dim[0]      Number of dimensions in database; usually 4. 
%           dim[1]      Image X dimension;  number of *pixels* in an image row. 
%           dim[2]      Image Y dimension;  number of *pixel rows* in slice. 
%           dim[3]      Volume Z dimension; number of *slices* in a volume. 
%           dim[4]      Time points; number of volumes in database
%       */
%       float intent_p1;   % char vox_units[4];   /* 16 + 4       */
%       float intent_p2;   % char cal_units[8];   /* 20 + 4       */
%       float intent_p3;   % char cal_units[8];   /* 24 + 4       */
%       short int intent_code;   % short int unused1;   /* 28 + 2 */
%       short int datatype;              /* 30 + 2          */
%       short int bitpix;                /* 32 + 2          */
%       short int slice_start;   % short int dim_un0;   /* 34 + 2 */
%       float pixdim[8];                 /* 36 + 32         */
%	/*
%		pixdim[] specifies the voxel dimensions:
%		pixdim[1] - voxel width, mm
%		pixdim[2] - voxel height, mm
%		pixdim[3] - slice thickness, mm
%		pixdim[4] - volume timing, in msec
%					..etc
%	*/
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

dime.dim        = fread(fid,8,'int16')';
dime.intent_p1  = fread(fid,1,'float32')';
dime.intent_p2  = fread(fid,1,'float32')';
dime.intent_p3  = fread(fid,1,'float32')';
dime.intent_code = fread(fid,1,'int16')';
dime.datatype   = fread(fid,1,'int16')';
dime.bitpix     = fread(fid,1,'int16')';
dime.slice_start = fread(fid,1,'int16')';
dime.pixdim     = fread(fid,8,'float32')';
dime.vox_offset = fread(fid,1,'float32')';
dime.scl_slope  = fread(fid,1,'float32')';
dime.scl_inter  = fread(fid,1,'float32')';
dime.slice_end  = fread(fid,1,'int16')';
dime.slice_code = fread(fid,1,'char')';
dime.xyzt_units = fread(fid,1,'char')';
dime.cal_max    = fread(fid,1,'float32')';
dime.cal_min    = fread(fid,1,'float32')';
dime.slice_duration = fread(fid,1,'float32')';
dime.toffset    = fread(fid,1,'float32')';
dime.glmax      = fread(fid,1,'int32')';
dime.glmin      = fread(fid,1,'int32')';

return					% image_dimension


%---------------------------------------------------------------------
function [ hist ] = data_history(fid,filetype)

%  Original header structures /* (NIfTI format) */  
%  struct data_history  
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

hist.descrip     = deblank(fread(fid,80,'*char')');
hist.aux_file    = deblank(fread(fid,24,'*char')');

%% Read NIfTI fields
if any(filetype==[1 2])
  hist.qform_code  = fread(fid,1,'int16')';
  hist.sform_code  = fread(fid,1,'int16')';
  hist.quatern_b   = fread(fid,1,'float32')';
  hist.quatern_c   = fread(fid,1,'float32')';
  hist.quatern_d   = fread(fid,1,'float32')';
  hist.qoffset_x   = fread(fid,1,'float32')';
  hist.qoffset_y   = fread(fid,1,'float32')';
  hist.qoffset_z   = fread(fid,1,'float32')';
  hist.srow_x      = fread(fid,4,'float32')';
  hist.srow_y      = fread(fid,4,'float32')';
  hist.srow_z      = fread(fid,4,'float32')';
  hist.intent_name = deblank(fread(fid,16,'*char')');
  hist.magic       = deblank(fread(fid,4,'*char')');
else % Read Analyze 7.5 fields
  hist.orient      = fread(fid, 1,'*char');
  hist.originator  = fread(fid,5,'*int16')';
  hist.generated   = fread(fid,10,'*char')';
  hist.scannum     = fread(fid,10,'*char')';
  hist.patient_id  = fread(fid,10,'*char')';
  hist.exp_date    = fread(fid,10,'*char')';
  hist.exp_time    = fread(fid,10,'*char')';
  hist.hist_un0    = fread(fid, 3,'*char')';
  hist.views       = fread(fid, 1,'*int32');
  hist.vols_added  = fread(fid, 1,'*int32');
  hist.start_field = fread(fid, 1,'*int32');
  hist.field_skip  = fread(fid, 1,'*int32');
  hist.omax        = fread(fid, 1,'*int32');
  hist.omin        = fread(fid, 1,'*int32');
  hist.smax        = fread(fid, 1,'*int32');
  hist.smin        = fread(fid, 1,'*int32');
  
  % Read also the magic field
  fseek(fid,344,'bof');
  hist.magic       = deblank(fread(fid,4,'*char')');
end

fseek(fid,253,'bof');
hist.originator  = fread(fid, 5,'int16')';

return					% data_history



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Read NIfTI data
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,hdr,msg]=l_ReadData(hdr,machine,filename,filetype)

img_idx = [];
old_RGB = 0;
msg='';

%% Use little endian byte ordering by default
if isempty(machine)
  machine='ieee-le';
end

[fpath,fname,fext]=fileparts(filename);
if any(filetype==[0 1])
  filename = [fpath,filesep,fname,'.img'];
  
  %% Make sure that the *.img file is found
  if exist(filename,'file')~=2
    data=[];
    hdr=[];
    msg = {'Could not find data file',...
           ['"' filename '"']};
    return
  end
else
  filename = [fpath,filesep,fname,'.nii'];
end
  
try
  [data,hdr.FileHeader] = read_image(hdr.FileHeader,filetype,filename,machine,img_idx, ...
                          old_RGB);
catch
  data=[];
  hdr=[];
  msg = lasterr;
end

return


%---------------------------------------------------------------------
function [img,hdr] = read_image(hdr, filetype,filename,machine,img_idx,old_RGB)

%switch filetype
% case {0, 1}
%  fn = [fileprefix '.img'];
% case 2
%  fn = [fileprefix '.nii'];
%end

fid = fopen(filename,'r',machine);

if fid < 0,
  msg = sprintf('Cannot open file %s.',filename);
  error(msg);
end

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
%  1792 Complex128, 2 float64  (Unsupported, bitpix=128) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128 
%  2048 Complex256, 2 float128 (Unsupported, bitpix=256) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128 
%
switch hdr.dime.datatype
 case   1,
  hdr.dime.bitpix = 1;  precision = 'ubit1';
 case   2,
  hdr.dime.bitpix = 8;  precision = 'uint8';
 case   4,
  hdr.dime.bitpix = 16; precision = 'int16';
 case   8,
  hdr.dime.bitpix = 32; precision = 'int32';
 case  16,
  hdr.dime.bitpix = 32; precision = 'float32';
 case  64,
  hdr.dime.bitpix = 64; precision = 'float64';
 case 128,
  hdr.dime.bitpix = 24; precision = 'uint8';
 case 256 
  hdr.dime.bitpix = 8;  precision = 'int8';
 case 511 
  hdr.dime.bitpix = 96; precision = 'float32';
 case 512 
  hdr.dime.bitpix = 16; precision = 'uint16';
 case 768 
  hdr.dime.bitpix = 32; precision = 'uint32';
 case 1024
  hdr.dime.bitpix = 64; precision = 'int64';
 case 1280
  hdr.dime.bitpix = 64; precision = 'uint64';
 otherwise
  error('This datatype is not supported'); 
end

%  move pointer to the start of image block
%
switch filetype
 case {0, 1}
  fseek(fid, 0, 'bof');
 case 2
  fseek(fid, hdr.dime.vox_offset, 'bof');
end

%  Load whole image block for old Analyze format, or binary image,
%  or img_idx is empty; otherwise, load images that are specified
%  in img_idx
%
%  For binary image, we have to read all because pos can not be
%  seeked in bit and can not be calculated the way below.
%
if filetype == 0 | hdr.dime.datatype == 1 | isempty(img_idx)
  img = fread(fid,inf,sprintf('*%s',precision));
else
  img = [];
  
  for i=1:length(img_idx)
    
    %  Precision of value will be read in img_siz times, where
    %  img_siz is only the dimension size of an image, not the
    %  byte storage size of an image.
    %
    img_siz = prod(hdr.dime.dim(2:4));

    %  Position is seeked in bytes. To convert dimension size
    %  to byte storage size, hdr.dime.bitpix/8 will be
    %  applied.
    %
    pos = (img_idx(i) - 1) * img_siz * hdr.dime.bitpix/8;
    
    if filetype == 2
      fseek(fid, pos + hdr.dime.vox_offset, 'bof');
    else
      fseek(fid, pos, 'bof');
    end

    %  fread will read precision of value in img_siz times
    %
    img = [img fread(fid,img_siz,sprintf('*%s',precision))];
  end
end

fclose(fid);

%  Update the global min and max values
%hdr.dime.glmax = max(double(img(:)));
%hdr.dime.glmin = min(double(img(:)));


if isempty(img_idx)
  img_idx = 1:hdr.dime.dim(5);
end

if old_RGB & hdr.dime.datatype == 128 & hdr.dime.bitpix == 24
  img = squeeze(reshape(img, [hdr.dime.dim(2:3) 3 hdr.dime.dim(4) length(img_idx)]));
  img = permute(img, [1 2 4 3 5]);
elseif hdr.dime.datatype == 128 & hdr.dime.bitpix == 24
  img = squeeze(reshape(img, [3 hdr.dime.dim(2:4) length(img_idx)]));
  img = permute(img, [2 3 4 1 5]);
elseif hdr.dime.datatype == 511 & hdr.dime.bitpix == 96
  img = double(img);
  img = (img - min(img))/(max(img) - min(img));
  img = squeeze(reshape(img, [3 hdr.dime.dim(2:4) length(img_idx)]));
  img = permute(img, [2 3 4 1 5]);
else
  %img = squeeze(reshape(img, [hdr.dime.dim(2:4) length(img_idx)]));
  img = reshape(img, [hdr.dime.dim(2:4) length(img_idx)]);
  img = flipdim(permute(img,[2 1 3 4]),1);
end

if ~isempty(img_idx)
  hdr.dime.dim(5) = length(img_idx);
end

%% Scale data if requested
if isfield(hdr.dime,'scl_slope') && hdr.dime.scl_slope~=0
  if isfield(hdr.dime,'scl_inter')
	if hdr.dime.scl_slope~=1 && hdr.dime.scl_inter~=0
	  img = double(img);
	  img = hdr.dime.scl_slope*img+hdr.dime.scl_inter;
	end
  else
	if hdr.dime.scl_slope~=1
	  img = double(img);
	  img = hdr.dime.scl_slope*img;
	end
  end
end

return						% read_image
