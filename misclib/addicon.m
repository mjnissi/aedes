function imdata=addicon(filename,fieldname)
% ADDICON
%
% Synopsis: 
%
% Description:
%
% Examples:
%
% See also:
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

if ~iscell(filename) && ~iscell(fieldname)
  filename = {filename};
  fieldname={fieldname};
end

% Load icon cdata
S=load('aedes_cdata.mat');
cdata = S.cdata;
  
for ii=1:length(filename)  
  % Read image
  imdata = imread(filename{ii});
  if length(size(imdata))<3
    continue
  end
  imdata=[];
  imdata = imread(filename{ii},'Backgroundcolor',[236 233 216]./255);
  
  % Scale image data and convert to double
  imdata = double(imdata);
  sz=size(imdata);
  
  % Add NaN:s
  for jj=1:sz(1)
    for kk=1:sz(2)
      if all(squeeze(imdata(jj,kk,:))'==[236 233 216])
        imdata(jj,kk,1)=NaN;
        imdata(jj,kk,2)=NaN;
        imdata(jj,kk,3)=NaN;
      end
    end
  end
  imdata = imdata/max(max(max(imdata)));
  %imdata(imdata==0)=NaN;
  
  cdata.(fieldname{ii}) = imdata;
end


% Show in figure menubar
fh=figure('menubar','none');
fl_names = fieldnames(cdata);
for ii=1:length(fl_names)
  if isnumeric(cdata.(fl_names{ii}))
    h=uipushtool('cdata',cdata.(fl_names{ii}));
  end
end
