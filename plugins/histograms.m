function histograms(DATA,ROI,AddInfo)
% HISTAGRAMS - Calculate histograms for image and ROIs (Aedes plugin)
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


%% Check if the data to be analyzed is valid
if AddInfo.isDataMixed
  errordlg('Sorry. This plugin does not work with mixed data.','Error',...
    'modal');
  return
end

% Default number of bins
nBins = 20;

%% Prompt for number of bins
resp = aedes_inputdlg('Please, input the number of bins',...
  'Number of bins',num2str(nBins));
if isempty(resp)
  % Canceled
  return
end
nBins = str2num(resp{1});

% Scale data using CLim values
data = double(DATA{1}.FTDATA);
data(data<AddInfo.Clim(1))=AddInfo.Clim(1);
data(data>AddInfo.Clim(2))=AddInfo.Clim(2);

%% Plot histogram of the whole data
[n_data,x_data]=hist(data(:),nBins);
fh=figure;
ax=axes('parent',fh);
bar(ax,x_data,n_data);
title(ax,'Whole data (scaled to CLim)');


%% Plot histograms of ROIs
if ~isempty(ROI)
  fh=figure;
  if length(ROI)<=3
    nCols = length(ROI);
    nRows = 1;
  elseif length(ROI)<=9
    nCols = 3;
    nRows = ceil(length(ROI)/nCols);
  elseif length(ROI)>9
    nCols = 4;
    nRows = ceil(length(ROI)/nCols);
  end
  
  for ii=1:length(ROI)
    data = DATA{1}.FTDATA(ROI(ii).voxels{1});
    [n_data,x_data]=hist(data,nBins);
    ax=subplot(nRows,nCols,ii,'align','parent',fh);
    %ax=axes('parent',fh);
    bar(ax,x_data,n_data);
    title(ax,['ROI ',num2str(ii),' (',ROI(ii).label,')']);
  end
end







