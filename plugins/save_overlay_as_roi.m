function save_overlay_as_roi(DATA,ROI,AddInfo)

if ~isfield(AddInfo,'ImageOverlay')
	errordlg('Image Overlay not defined.','Error','modal')
	return
end

% Get Threshold and clim
clim = round(((AddInfo.ImageOverlayClim-AddInfo.ImOverlayMin)*256)./...
	(AddInfo.ImOverlayMax-AddInfo.ImOverlayMin));
if clim(1)==clim(2)
	if clim(1)==0
		clim(2)=1;
	elseif clim(1)==256
		clim(1)=255;
	end
end
thold = round(((AddInfo.ImageOverlayThold-AddInfo.ImOverlayMin)*256)/...
	(AddInfo.ImOverlayMax-AddInfo.ImOverlayMin));

% Convert indexed image to RGB image
overlay=AddInfo.ImageOverlay(:,:,:,AddInfo.CurrentVol);
overlay=double(overlay);

% Get thresholded alpha indices
if AddInfo.ImageOverlayTholdDirPos==1
	overlay_th = overlay<thold;
else
	overlay_th = overlay>thold;
end

% Get clim alpha indices
overlay_clim = ( overlay>=clim(1) & overlay<=clim(2) );

roi_mask = false(size(overlay));
roi_mask(overlay_clim) = true;
roi_mask(overlay_th) = false;

ROI = [];
ROI.voxels{1} = roi_mask;
ROI.fpath{1} = '';
ROI.fname{1} = '';
ROI.label = 'overlay';
ROI.color = [255 0 0];
DateTime = datestr(now);
RotateFlip = {};

% Get default directory
if ispref('Aedes','GetOverlayFileDir')
	default_dir = [getpref('Aedes','GetOverlayFileDir'),'untitled.roi'];
else
	default_dir = [pwd,filesep,'untitled.roi'];
end

% Save ROI
[filename,filepath] = uiputfile({'*.roi','Aedes ROI-files (*.roi)'},...
	'Save ROI as',default_dir);
if isequal(filename,0)
	% Canceled
	return
end

[fp,fn,fe]=fileparts(filename);

ROI.label = fn;

save([filepath,filename],'ROI','DateTime','RotateFlip');