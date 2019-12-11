function motion_estimate(DATA,ROI,AddInfo)
% MOTION_ESTIMATE - Estimate motion using center of mass
%   (this is an Aedes plugin)
%
% Synopsis: 
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


% Calculate index vector for center of mass
data_sz = [size(DATA{1}.FTDATA,1) ...
  size(DATA{1}.FTDATA,2) ...
  size(DATA{1}.FTDATA,3) ...
  size(DATA{1}.FTDATA,4)];
n1 = 1:data_sz(1);
n2 = 1:data_sz(2);
n3 = 1:data_sz(3);
nSlices = data_sz(3);
nVols = data_sz(4);
N2 = repmat(n2,data_sz(1),1);
N2 = N2(:);
N1 = n1.';
N1 = repmat(N1,data_sz(2),1);
N3 = n3;
N3 = repmat(N3,length(N2),1);
N3 = N3(:);
tmp = [N1 N2];
C = [repmat(tmp,nSlices,1) N3];

% Calculate center of mass to all volumes
if AddInfo.isDataMixed
  COM = zeros(length(DATA),3);
  for ii=1:length(DATA)
    data=double(DATA{ii}.FTDATA);
    COM(ii,:) = data(:).'*C/sum(data(:),'double');
  end
else
  COM = zeros(nVols,3);
  for ii=1:nVols
    data=double(DATA{1}.FTDATA(:,:,:,ii));
    COM(ii,:) = data(:).'*C/sum(data(:),'double');
  end
end

% Plot the center of mass indices relative to first volume
fh=figure;
ax=axes;
plot(ax,1:nVols,COM(:,1)-COM(1,1),'b');
hold on
plot(ax,1:nVols,COM(:,2)-COM(1,2),'g');
plot(ax,1:nVols,COM(:,3)-COM(1,3),'r');
hold off

xlabel(ax,'Volume number');
ylabel(ax,'Relative movement (pixels)');
title(ax,'Relative center of mass displacement')
legend(ax,'Z-dir','Y-dir','X-dir')


