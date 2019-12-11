function plot_profile(DATA,ROI,AddInfo)
% PLOT_PROFILE - Plot profile
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


if AddInfo.isDataMixed
  fh=figure;
  imagesc(DATA{AddInfo.CurrentSlice}.FTDATA,AddInfo.Clim)
  colormap gray
  axis image
  [cx,cy,c]=improfile;
  line(cx,cy,'color','c','linewidth',2)
  fh2=figure;
  plot(c)
else
  fh=figure;
  imagesc(DATA{1}.FTDATA(:,:,AddInfo.CurrentSlice(1)),AddInfo.Clim)
  colormap gray
  axis image
  [cx,cy,c]=improfile;
  line(cx,cy,'color','c','linewidth',2)
  fh2=figure;
  plot(c)
end
