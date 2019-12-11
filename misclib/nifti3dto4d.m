function nifti3dto4d(infiles,outfile,varargin)
% NIFTI3DTO4D - Convert a series of 3D NIfTI volumes to a single 4D NIfTI
%               file
%
% Synopsis: 
%	data = nifti3dto4d(infiles,outfile,param1,value1,param2,value2,...)
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

% Defaults ---------------------------------------
showWbar = true;
TR = [];

% Parse input arguments --------------------------
if nargin==0 || isempty(infiles)
  % Prompt files to read
  [fname,fpath,findex] = aedes_juigetfiles({'*.nii','NIfTI files (*.nii)';...
	'*.*','All Files (*.*)'},'Select volumes',pwd);
  
  if isequal(fname,0)
	% Canceled
	return
  end
  tmp={fpath{:};fname{:}};
  infiles = cell(1,size(tmp,2));
  for ii=1:size(tmp,2)
	infiles{ii} = [tmp{:,ii}];
  end
end

if nargin<2 || isempty(outfile)
  % Select output file
  [fn,fp] = uiputfile({'*.nii','NIfTI files (*.nii)';...
	'*.*','All Files (*.*)'},'Select output file',...
	'4D_nifti.nii');
  if isequal(fn,0)
    % Canceled
    return
  end
  outfile = fullfile(fp,fn);
end

% Read the first image and get header parameters
tmp_data = aedes_read_nifti(infiles{1});
TimeUnits = tmp_data.HDR.timeunits;
xyzunits = tmp_data.HDR.xyzunits;
byteorder = tmp_data.HDR.byteorder;
vox_size = tmp_data.HDR.FileHeader.dime.pixdim(2:4);

% Try to determine the TR value from the description field (SPM writes this
% information here when converting DICOMs, NOTE: THIS MAY FAIL!!!)
if isempty(TR)
  description = tmp_data.HDR.FileHeader.hist.descrip;
  ind=strfind(tmp_data.HDR.FileHeader.hist.descrip,'TR=');
  TR = str2num(tmp_data.HDR.FileHeader.hist.descrip(ind+3:ind+6))/1000;
end
  


data_sz(1)=size(tmp_data.FTDATA,1);
data_sz(2)=size(tmp_data.FTDATA,2);
data_sz(3)=size(tmp_data.FTDATA,3);

% Parse varargin
for ii=1:2:length(varargin)
  switch lower(varargin{ii})
	case 'tr'
	  TR = varargin{ii+1};
	case 'timeunits'
	  TimeUnits = varargin{ii+1};
	case 'wbar'
	  if strcmpi(varargin{ii+1},'off')
		showWbar = false;
	  end
	otherwise
	  error('Unknown property "%s"',lower(varargin{ii}))
  end
end
vox_size(end+1) = TR;

% Allocate space for the 4D nifti
DATA.DataFormat = tmp_data.DataFormat;
DATA.HDR.FileHeader = tmp_data.HDR.FileHeader;
DATA.FTDATA = zeros(data_sz(1),data_sz(2),data_sz(3),length(infiles),...
  class(tmp_data.FTDATA));
DATA.FTDATA(:,:,:,1) = tmp_data.FTDATA;

if showWbar
  wbh=aedes_wbar(1/length(infiles),...
	sprintf('Reading 3D NIfTI files...1/%d',length(infiles)));
end
for ii=2:length(infiles)
  % Read the nifti file
  tmp_data = aedes_read_nifti(infiles{ii});
  sz(1)=size(tmp_data.FTDATA,1);
  sz(2)=size(tmp_data.FTDATA,2);
  sz(3)=size(tmp_data.FTDATA,3);
  if ~isequal(sz,data_sz)
	error('Data size doesn''t match in file %s%s',infiles{ii})
  end
  DATA.FTDATA(:,:,:,ii) = tmp_data.FTDATA;
  if showWbar
	aedes_wbar(ii/length(infiles),wbh,...
	  sprintf('Reading 3D NIfTI files...%d/%d',ii,length(infiles)))
  end
end
if showWbar
  close(wbh)
  drawnow
end

% % High pass filtering
% [B2,A2]=butter(6,0.02,'high');
% counter=1;
% sz=size(DATA);
% nVox = prod(sz(1:3));
% fprintf(1,'\nHigh-pass filtering:\n');
% for ii=1:size(DATA,1)
%   for kk=1:size(DATA,2)
% 	for zz=1:size(DATA,3)
% 	  tmp = squeeze(double(DATA(ii,kk,zz,:)));
% 	  DATA(ii,kk,zz,:) = filtfilt(B2,A2,tmp)+mean(tmp);
% 	  if counter~=1
% 		fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
% 	  end
% 	  fprintf(1,'Filtering: %03d %%',round(100*counter/nVox))
% 	  counter=counter+1;
% 	end
%   end
% end
% fprintf(1,'\n');

% Write the resulting 4D NIfTI file
if showWbar
  cwh = aedes_calc_wait(sprintf('Writing 4D NIfTI file:\n%s',outfile));
end
aedes_write_nifti(DATA,outfile,...
  'VoxelSize',vox_size,...
  'XYZUnits',xyzunits,...
  'TimeUnits',TimeUnits,...
  'Description',description,...
  'machine',byteorder);
if showWbar
  delete(cwh)
end



