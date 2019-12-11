function fmri_analysis(DATA,ROI,AddInfo)
%
% This Aedes plugin does a basic fMRI analysis.
%

if AddInfo.isDataMixed || ndims(DATA{1}.FTDATA)<4
	error('Input data has to have at least 4 dimensions.')
end

% Get brain mask if available
if ~isempty(ROI) && any(strcmpi({ROI(:).label},'mask'))
	roi_ind = find(strcmpi({ROI(:).label},'mask'));
	mask = ROI(roi_ind).voxels{1}(:,:,:,1);
else
	mask = [];
end

% Launch fMRI GUI
aedes_fmri_gui(DATA{1},mask)

