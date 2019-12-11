function resting_state_fc(DATA,ROI,AddInfo)
%
% An Aedes plugin that calculates resting state funtional connectivity from
% seed ROIs.


% Only accept 3D or 4D data
if AddInfo.isDataMixed || ndims(DATA{1}.FTDATA)==2
  errordlg('Only 3D or 4D can be used with this plugin.','Error',...
    'modal');
  return
end

data=DATA{1}.FTDATA;


% Get the ROIs from current slice/volume
if ndims(data)==3
  % Data volume upon which to show the maps
  vol = data(:,:,1);
  for ii=1:length(ROI)
    ROI(ii).voxels{1} = ROI(ii).voxels{1}(:,:,AddInfo.CurrentSlice(3));
  end
else
  vol = data(:,:,:,1);
  for ii=1:length(ROI)
    ROI(ii).voxels{1} = ROI(ii).voxels{1}(:,:,:,AddInfo.CurrentVol);
  end
end

% Do global normalization
data = fmri_global_norm(data);

% Spatially smooth data
data = aedes_fmri_smooth(data,[2 2 1]);

% Do fMRI filtering
data = fmri_filter(data,2,'detrending','on',...
  'lowpass',0.08);

% Correlate seeds
corrmap = fmri_corr(data,ROI);

% Display results in a new Aedes window
dt = repmat(vol,[1 1 1 length(corrmap)]);
rs_maps = zeros(size(corrmap(1).ccc,1),size(corrmap(1).ccc,2),...
  size(corrmap(1).ccc,3),length(corrmap));
for kk=1:length(corrmap)
  rs_maps(:,:,:,kk) = corrmap(kk).ccc;
end
aedes(dt,[],rs_maps)









