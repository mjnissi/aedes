function write_difference_images(DATA,ROI,AddInfo)
% Aedes plugin for writing difference images for ASL

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

if AddInfo.isDataMixed
  if rem(length(DATA),2)~=0
	return
  end
  
  % Ask a directory to save the files
  dirname = uigetdir;
  if isequal(dirname,0)
	return
  else
	dirname = [dirname,filesep];
  end
  
  data_new = cell(1,length(DATA)/2);
  fnames = cell(1,length(DATA)/2);
  count = 1;
  for ii=1:2:length(DATA)
	data_new{count} = DATA{ii}.FTDATA-DATA{ii+1}.FTDATA;
	[fp1,fn1,fe1]=fileparts(DATA{ii}.HDR.fname);
	[fp2,fn2,fe2]=fileparts(DATA{ii+1}.HDR.fname);
	fnames{count}=[fn1,' - ',fn2,'.nii'];
	count=count+1;
  end
  
  
  % write *.nii files
  h=aedes_wbar(0,sprintf('Writing File...\n%s',' '));
  for ii=1:length(data_new)
	aedes_wbar(ii/length(data_new),h,...
	  sprintf('Writing File...\n%s',[dirname,fnames{ii}]))
	aedes_write_nifti(data_new{ii},[dirname,fnames{ii}]);
  end
  delete(h)
else
  return
end
