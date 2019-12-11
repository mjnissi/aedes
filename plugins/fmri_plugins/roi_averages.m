function roi_averages(DATA,ROI,AddInfo)
%
% This Aedes plugin calculates and plots average time series from ROIs
%

if AddInfo.isDataMixed
  % Make data a 4D-matrix if it is loaded in aedes as a stack of 2D or 3D images
  DATA2{1}=DATA{1};
  DATA2{1}.FTDATA = zeros([size(DATA{1}.FTDATA,1),...
    size(DATA{1}.FTDATA,2),size(DATA{1}.FTDATA,3),length(DATA)],...
    class(DATA{1}.FTDATA));
  for ii=1:length(DATA)
    DATA2{1}.FTDATA(:,:,:,ii)=DATA{ii}.FTDATA;
  end
  DATA=DATA2;
  CurrentVol=AddInfo.CurrentSlice;
else
  CurrentVol=AddInfo.CurrentVol;
end

CurrentVol=AddInfo.CurrentVol;

% Check that there are ROIs defined at all
if isempty(ROI)
  % No ROIs defined, nothing to calculate
  errordlg('No ROIs defined. At least one ROI must be drawn to the data to allow average ROI time series calculation.',...
    'ROIs not defined.','modal');
  return
end

% Check that there are ROIs defined in this volume
RoiInds = [];
for ii=1:length(ROI)
  if any(any(any(ROI(ii).voxels{1}(:,:,:,CurrentVol))))
    RoiInds(end+1)=ii;
  end
end

% Return if there were no ROIs in the current volume
if isempty(RoiInds)
  errordlg('No ROIs defined in the current volume. At least one ROI must be drawn to the current volume to allow average ROI time series calculation.',...
    'ROIs not defined in current volume.','modal');
  return
end

% Plot time series from ROIs
nScans = size(DATA{1}.FTDATA,4);
fh=figure;
if length(RoiInds)<=3
  nRows = length(RoiInds);
  nCols=1;
elseif length(RoiInds)==4
  nRows = 2;
  nCols=2;
elseif length(RoiInds)<10
  nRows = 3;
  nCols = ceil(length(RoiInds)/nRows);
else
  nRows = 4;
  nCols = ceil(length(RoiInds)/nRows);
end
for kk=RoiInds
  
  % Get time-series indices from ROI
  ind = repmat(ROI(kk).voxels{1}(:,:,:,CurrentVol),[1 1 1 size(DATA{1}.FTDATA,4)]);
  
  % Mean EPI time series
  ts_data = reshape(double(DATA{1}.FTDATA(ind)),[],size(DATA{1}.FTDATA,4));
  ts_data = ts_data.';
  ts_data = detrend(ts_data)+repmat(mean(ts_data),size(ts_data,1),1);
  
  % Normalize mean to bold-%
  mean_ts_data = mean(ts_data,2);
  mean_ts_data = (mean_ts_data./mean(mean_ts_data)-1)*100;
  
  data_trend=aedes_trendest(double(mean_ts_data),10);
  
  % Plot results
  ax=subplot(nRows,nCols,kk,'align','parent',fh);
  line(1:length(mean_ts_data),mean_ts_data,'color','k',...
    'parent',ax);
  line(1:length(mean_ts_data),data_trend,...
    'color','r','linewidth',2,'parent',ax);
  title(['Time series for ROI: ',ROI(kk).label])
  ylabel(ax,'BOLD-%');
  set(ax,'xlim',[0 length(mean_ts_data)],...
    'ylim',[min(mean_ts_data)-min(mean_ts_data)*0.05 ...
    max(mean_ts_data)+max(mean_ts_data)*0.05]);
end
