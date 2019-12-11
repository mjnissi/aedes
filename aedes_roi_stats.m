function Res = aedes_roi_stats(DATA,ROI)
% AEDES_ROI_STATS - Calculate ROI statistics
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


Res = [];

%% Check input arguments
if nargin<2
  error('Two input arguments are required')
elseif nargin>2
  error('Too many input arguments')
end

% Ensure that DATA is a cell array
if ~iscell(DATA)
  DATA = {DATA};
end

% if ROI is empty, return empty
if isempty(ROI)
  return
end

% Check if data is in mixed form
if length(DATA)>1 || ndims(DATA{1}.FTDATA)==2
  isDataMixed = true;
else
  isDataMixed = false;
end

% Generate timestamp
DateTime = datestr(now);

% Generate file info
for ii=1:length(DATA)
  FileInfo.DataFileName{ii} = DATA{ii}.HDR.fname;
  FileInfo.DataPathName{ii} = DATA{ii}.HDR.fpath;
end
Res.DateTime = DateTime;
Res.FileInfo = FileInfo;
Res.Stat = [];

for ii=1:length(ROI)
  
  if isDataMixed
    Res.Stat(ii).FileName = {};
  end 
  
  for kk=1:length(DATA)
    
    % Statistics for mixed data type
    if isDataMixed
      % Get ROI data
      data = DATA{kk}.FTDATA(ROI(ii).voxels{kk});
      data = double(data);
      Res.Stat(ii).FileName{kk,1} = DATA{kk}.HDR.fpath;
      Res.Stat(ii).FileName{kk,2} = DATA{kk}.HDR.fname;
      if kk==1
        Res.Stat(ii).isMixed = true;
        Res.Stat(ii).Label = ROI(ii).label;
      end
      
      if isempty(data)
        Res.Stat(ii).Mean(kk) = NaN;
        Res.Stat(ii).Std(kk) = NaN;
        Res.Stat(ii).Sum(kk) = NaN;  
        Res.Stat(ii).Max(kk) = NaN;
        Res.Stat(ii).Min(kk) = NaN;
        Res.Stat(ii).PixelCount(kk) = 0;
      else
        Res.Stat(ii).Mean(kk) = mean(data);
        Res.Stat(ii).Std(kk) = std(data);
        Res.Stat(ii).Sum(kk) = sum(data);
        Res.Stat(ii).Max(kk) = max(data);
        Res.Stat(ii).Min(kk) = min(data);
        Res.Stat(ii).PixelCount(kk) = length(data);
      end
        
    else % Statistics for normal type data
      
      %% Overall results for current ROI
      data=DATA{kk}.FTDATA(ROI(ii).voxels{kk});
      data = double(data);
      Res.Stat(ii).isMixed = false;
      Res.Stat(ii).Label = ROI(ii).label;
      Res.Stat(ii).Mean = mean(data);
      Res.Stat(ii).Std = std(data);
      Res.Stat(ii).Sum = sum(data);
      Res.Stat(ii).Max = max(data);
      Res.Stat(ii).Min = min(data);
      Res.Stat(ii).PixelCount = length(data);
      if isempty(data)
        Res.Stat(ii).Mean = NaN;
        Res.Stat(ii).Std = NaN;
        Res.Stat(ii).Sum = NaN;
        Res.Stat(ii).Max = NaN;
        Res.Stat(ii).Min = NaN;
        Res.Stat(ii).PixelCount = 0;
      end
      
      %% Calculate results in X direction
      Res.Stat(ii).XD.Mean = [];
      Res.Stat(ii).XD.Std = [];
      Res.Stat(ii).XD.Sum = [];
      Res.Stat(ii).XD.Max = [];
      Res.Stat(ii).XD.Min = [];
      Res.Stat(ii).XD.PixelCount = [];

			for jj=1:size(ROI(ii).voxels{kk},1)
        roix=ROI(ii).voxels{kk}(jj,:,:,:);
        datax = DATA{kk}.FTDATA(jj,:,:,:);
        tmpx=datax(roix);
        tmpx=double(tmpx);
        if isempty(tmpx)
          Res.Stat(ii).XD.Mean(end+1) = NaN;
          Res.Stat(ii).XD.Std(end+1) = NaN;
          Res.Stat(ii).XD.Sum(end+1) = NaN;
          Res.Stat(ii).XD.Max(end+1) = NaN;
          Res.Stat(ii).XD.Min(end+1) = NaN;
          Res.Stat(ii).XD.PixelCount(end+1) = 0;
        else
          Res.Stat(ii).XD.Mean(end+1) = mean(tmpx);
          Res.Stat(ii).XD.Std(end+1) = std(tmpx);
          Res.Stat(ii).XD.Sum(end+1) = sum(tmpx);
          Res.Stat(ii).XD.Max(end+1) = max(tmpx);
          Res.Stat(ii).XD.Min(end+1) = min(tmpx);
          Res.Stat(ii).XD.PixelCount(end+1) = length(tmpx);
        end
      end
      
      %% Calculate results in Y direction
      Res.Stat(ii).YD.Mean = [];
      Res.Stat(ii).YD.Std = [];
      Res.Stat(ii).YD.Sum = [];
      Res.Stat(ii).YD.Max = [];
      Res.Stat(ii).YD.Min = [];
      Res.Stat(ii).YD.PixelCount = [];
      
      for jj=1:size(ROI(ii).voxels{kk},2)
        roiy=ROI(ii).voxels{kk}(:,jj,:,:);
        datay = DATA{kk}.FTDATA(:,jj,:,:);
        tmpy=datay(roiy);
        tmpy=double(tmpy);
        if isempty(tmpy)
          Res.Stat(ii).YD.Mean(end+1) = NaN;
          Res.Stat(ii).YD.Std(end+1) = NaN;
          Res.Stat(ii).YD.Sum(end+1) = NaN;
          Res.Stat(ii).YD.Max(end+1) = NaN;
          Res.Stat(ii).YD.Min(end+1) = NaN;
          Res.Stat(ii).YD.PixelCount(end+1) = 0;
        else
          Res.Stat(ii).YD.Mean(end+1) = mean(tmpy);
          Res.Stat(ii).YD.Std(end+1) = std(tmpy);
          Res.Stat(ii).YD.Sum(end+1) = sum(tmpy);
          Res.Stat(ii).YD.Max(end+1) = max(tmpy);
          Res.Stat(ii).YD.Min(end+1) = min(tmpy);
          Res.Stat(ii).YD.PixelCount(end+1) = length(tmpy);
        end
      end
      
      %% Calculate results in Z direction
      Res.Stat(ii).ZD.Mean = [];
      Res.Stat(ii).ZD.Std = [];
      Res.Stat(ii).ZD.Sum = [];
      Res.Stat(ii).ZD.Max = [];
      Res.Stat(ii).ZD.Min = [];
      Res.Stat(ii).ZD.PixelCount = [];
      
			for jj=1:size(ROI(ii).voxels{kk},3)
        roiz=ROI(ii).voxels{kk}(:,:,jj,:);
        dataz = DATA{kk}.FTDATA(:,:,jj,:);
        tmpz=dataz(roiz);
        tmpz = double(tmpz);
        if isempty(tmpz)
          Res.Stat(ii).ZD.Mean(end+1) = NaN;
          Res.Stat(ii).ZD.Std(end+1) = NaN;
          Res.Stat(ii).ZD.Sum(end+1) = NaN;
          Res.Stat(ii).ZD.Max(end+1) = NaN;
          Res.Stat(ii).ZD.Min(end+1) = NaN;
          Res.Stat(ii).ZD.PixelCount(end+1) = 0;
        else
          Res.Stat(ii).ZD.Mean(end+1) = mean(tmpz);
          Res.Stat(ii).ZD.Std(end+1) = std(tmpz);
          Res.Stat(ii).ZD.Sum(end+1) = sum(tmpz);
          Res.Stat(ii).ZD.Max(end+1) = max(tmpz);
          Res.Stat(ii).ZD.Min(end+1) = min(tmpz);
          Res.Stat(ii).ZD.PixelCount(end+1) = length(tmpz);
        end
      end
      
      %% Calculate results in V direction
      Res.Stat(ii).VD.Mean = [];
      Res.Stat(ii).VD.Std = [];
      Res.Stat(ii).VD.Sum = [];
      Res.Stat(ii).VD.Max = [];
      Res.Stat(ii).VD.Min = [];
      Res.Stat(ii).VD.PixelCount = [];
      
      for jj=1:size(ROI(ii).voxels{kk},4)
        roiv=ROI(ii).voxels{kk}(:,:,:,jj);
        datav = DATA{kk}.FTDATA(:,:,:,jj);
        tmpv=datav(roiv);
        tmpv=double(tmpv);
        if isempty(tmpv)
          Res.Stat(ii).VD.Mean(end+1) = NaN;
          Res.Stat(ii).VD.Std(end+1) = NaN;
          Res.Stat(ii).VD.Sum(end+1) = NaN;
          Res.Stat(ii).VD.Max(end+1) = NaN;
          Res.Stat(ii).VD.Min(end+1) = NaN;
          Res.Stat(ii).VD.PixelCount(end+1) = 0;
        else
          Res.Stat(ii).VD.Mean(end+1) = mean(tmpv);
          Res.Stat(ii).VD.Std(end+1) = std(tmpv);
          Res.Stat(ii).VD.Sum(end+1) = sum(tmpv);
          Res.Stat(ii).VD.Max(end+1) = max(tmpv);
          Res.Stat(ii).VD.Min(end+1) = min(tmpv);
          Res.Stat(ii).VD.PixelCount(end+1) = length(tmpv);
        end
      end
    end
  end
end

% - EOF -
