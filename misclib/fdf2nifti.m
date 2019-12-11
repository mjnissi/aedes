function fdf2nifti(indir,outname)
% FDF2NIFTI - Read Varian FDF files from a folder and output a single
% NIfTI image
%   
%
% Synopsis: 
%
% Description:
%
% Examples:
%
% See also:
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

% Defaults
promptInDir = true;
promptOutName = true;

if nargin == 1
	promptInDir = false;
elseif nargin == 2
	promptInDir = false;
	promptOutName = false;
elseif nargin > 2
	error('Too many input arguments.')
end

% Get default directory
if ispref('Aedes','GetDataFileDir')
	default_dir = getpref('Aedes','GetDataFileDir');
else
	default_dir = [pwd,filesep];
end

% Prompt for directory and file name
if promptInDir
	indir = uigetdir(default_dir,'Select .img directory');
	if isequal(indir,0)
		return
	end
end
if ~strcmp(indir(end),filesep)
	indir = [indir,filesep];
end
if promptOutName
	[fn,fp,fe] = uiputfile({'*.nii','NIfTI files (*.nii)'},...
		'Save NIfTI file as',[default_dir,'untitled.nii']);
	if isequal(fn,0)
		return
	end
	outname = [fp,fn];
end

% Get fdf files in the folder
d = dir(indir);
files = {d(~[d(:).isdir]).name};
fdf_files = regexp(files,'^.*\.fdf$');
if all(cellfun('isempty',regexp(files,'^(.*\.fdf)$')))
	error('Could not find any FDF-files in "%s".',indir)
end
fdf_files = files(~cellfun('isempty',regexp(files,'^(.*\.fdf)$')));

% Get indexes (slice, image echo)
tmp = regexprep(fdf_files,'^slice(\d+)image(\d+)echo(\d+)\.fdf$','$1 $2 $3;');
ind = str2num([tmp{:}]);

% Sort using image number
[tmp2,ind2] = sort(ind(:,2));
ind=ind(ind2,:);

nSlices = max(ind(:,1));
nImages = max(ind(:,2));
nEchos = max(ind(:,3));

% Read first fdf file for determining dimensions
data = aedes_readfdf([indir,fdf_files{1}]);
sz = size(data.FTDATA);

% Allocate space for NIfTI
nifti_data = zeros([sz(1),sz(2),nSlices,nImages],'single');


wbh = aedes_wbar(0/nImages,sprintf('Processing image %d/%d...',0,nImages));
for tt=1:nImages
	count = 1;
	for kk=1:nSlices
		for ii=1:nEchos
			fname = sprintf('slice%03dimage%03decho%03d.fdf',kk,tt,ii);
			data = aedes_readfdf([indir,fname]);
			nifti_data(:,:,count,tt) = data.FTDATA;
			count = count+1;
		end
	end
	aedes_wbar(tt/nImages,wbh,sprintf('Processing image %d/%d...',tt,nImages));
end
aedes_wbar(1,wbh,sprintf('Writing data in NIfTI format to\n%s',outname));


% Write the nifti data to file
aedes_write_nifti(nifti_data,outname);

% Close waitbar
close(wbh);




