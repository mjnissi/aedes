function nifti4dto3d(infile,outfile)
% NIFTI4DTO3D - Convert a single 4D NIfTI file into a series of 3D NIfTI
% files
%
% Synopsis: 
%	data = nifti4dto3d(infile,outfile,param1,value1,param2,value2,...)
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
TR = 0;

if nargin==0
  % Prompt for a file
  [fname,fpath,findex] = uigetfile({'*.nii','4D NIfTI Files (*.nii)';...
	'*.*','All Files (*.*)'},'Select 4D NIfTI file');
  if isequal(fname,0)
	% Canceled
	return
  end
  infile = fullfile(fpath,fname);
end

if nargin<2
  % Prompt for the outfile
  [fname,fpath,findex] = uiputfile({'*.nii','4D NIfTI Files (*.nii)';...
	'*.*','All Files (*.*)'},...
	'Select folder and file name for 3D NIfTI files',...
	'3dniftifiles');
  if isequal(fname,0)
	% Canceled
	return
  end
  outfile = fullfile(fpath,fname);
end
[fp,fn,fe]=fileparts(outfile);
outfpath = [fp,filesep];
outfname = fn;

% Read the 4D NIfTI file
DATA = aedes_read_nifti(infile);
nVol = size(DATA.FTDATA,4);

% Construct the header struct
FileHeader = DATA.HDR.FileHeader;
FileHeader.dime.pixdim(5)=0;


% Write new 3D NIfTI volumes
if showWbar
  wbh=aedes_wbar(1/nVol,...
	sprintf('Writing 3D NIfTI files...1/%d',nVol));
end
for ii=1:size(DATA.FTDATA,4)
  
  % Generate file name
  filename = sprintf('%s%s_%04d.nii',outfpath,outfname,ii);
  
  % Write NIfTI files...
  tmp_data.FTDATA = DATA.FTDATA(:,:,:,ii);
  tmp_data.HDR.FileHeader = FileHeader;
  done=aedes_write_nifti(tmp_data,filename);
  
  if ~done
	fprintf(1,['Could not write ',filename,'\n'])
  end
  
  if showWbar
	aedes_wbar(ii/nVol,wbh,...
	  sprintf('Writing 3D NIfTI files...%d/%d',ii,nVol))
  end
end
if showWbar
  close(wbh)
end




