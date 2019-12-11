function open_in_imtool(DATA,ROI,AddInfo)
%
% This Aedes plugin opens the current slice in the imtool that comes with 
% Matlab Image Processing Toolbox.
%

if AddInfo.isDataMixed
	imdata = DATA{AddInfo.CurrentSlice}.FTDATA;
else
	dt_ind = AddInfo.CurrentSlice;
	if AddInfo.AxView==0 || AddInfo.AxView==1
		imdata = DATA{1}.FTDATA(:,:,dt_ind(1));
	elseif AddInfo.AxView==3
		imdata = squeeze(DATA{1}.FTDATA(dt_ind(3),:,:));
	else
		imdata = squeeze(DATA{1}.FTDATA(:,dt_ind(2),:));
	end
end

% Colormap and clim
cmap = get(AddInfo.hFigure,'colormap');
clim = AddInfo.Clim;

% Open imtool
imtool(imdata,'Colormap',cmap,...
	'DisplayRange',clim,...
	'InitialMagnification','fit');
