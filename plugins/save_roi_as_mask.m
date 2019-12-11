function save_roi_as_mask(DATA,ROI,AddInfo)
% Aedes plugin for saving ROI masks

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

if isempty(ROI)
  errordlg('ROI(s) not defined. Bailing out!','ROI(s) not defined','modal')
  return
end

if AddInfo.isDataMixed
  % Slice mask
  CurrentSlice = AddInfo.CurrentSlice;
  data = ROI(AddInfo.CurrentROI).voxels{CurrentSlice};
  data = uint8(data);
else
  % Volume mask
  data = ROI(AddInfo.CurrentROI).voxels{1}(:,:,:,AddInfo.CurrentVol);
  data = uint8(data);
end

% Get default directory
try,
  default_dir = getpref('Aedes','PutDataFileDir');
catch
  if isunix
	default_dir = getenv('HOME');
  else
	default_dir = getenv('USERPROFILE');
  end
  if not(strcmpi(default_dir(end),filesep))
	default_dir(end+1)=filesep;
  end
end

% Write the NIfTI file
[fname,fpath,findex]=uiputfile({'*.nii','NIfTI-files (*.nii)';...
  '*.*','All Files (*.*)'},...
  'Save As...',[default_dir,'mask.nii']);
if isequal(fname,0) || isequal(fpath,0)
  % Canceled
  return
end
setpref('Aedes','PutDataFileDir',fpath);
[fp,fn,fe]=fileparts([fpath,fname]);
filename = [fp,filesep,fn,'.nii'];

% If the data is VNMR data, set voxel size to 1mmx1mm
if strcmpi(DATA{1}.DataFormat,'vnmr')
  voxelsize = [1 1 DATA{1}.PROCPAR.thk];
  xyzunits = 'mm';
  timeunits = 'sec';
  datatype = 'uint8';
elseif strncmpi(DATA{1}.DataFormat,'nifti',5)
%   voxelsize = DATA{1}.HDR.FileHeader.dime.pixdim(2:4);
%   xyzunits = DATA{1}.HDR.xyzunits;
%   timeunits = DATA{1}.HDR.timeunits;
%  datatype = 'uint8';
else
  xyzunits = '';
  timeunits = '';
  voxelsize = [];
  datatype = 'uint8';
end
  
% Write data
if strncmpi(DATA{1}.DataFormat,'nifti',5)
  DATA{1}.FTDATA=data;
  datatype = 'uint8';
  [done,msg]=aedes_write_nifti(DATA{1},filename,...
	'DataType',datatype);
else
  [done,msg]=aedes_write_nifti(data,filename,...
	'VoxelSize',voxelsize,'XYZUnits',xyzunits,...
	'TimeUnits',timeunits);
end

if not(done)
  if not(iscell(msg))
	msg = {msg};
  end
  errordlg({'Could not save file "',filename,'"!',...
	'',msg{:}},'Error while saving mask!','modal')
else
  fprintf(1,'Successfully saved mask to "%s"...\n',filename);
end

