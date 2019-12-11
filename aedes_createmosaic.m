function aedes_createmosaic(DATA,ROI,overlay,varargin)
% AEDES_CREATEMOSAIC - Create a figure containing mosaic of MR images
%
%
% Synopsis:
%
% Description:
%
% Examples:
%
% See also:
%       AEDES, AEDES_EXPORT_GUI

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

% Defaults
sliceDir = 3; % Z-direction
roi_ind = 1:length(ROI);
drawSliceNbr = false;
drawColorbar = true;
drawFileName = false;
drawOverlay = true;
drawRois = true;
colmap = gray(256);
sliceInd = [];
rows = [];
cols = [];
Clim_in = [];
roiTransp = 1;
showRoiEdges = false;
Dat.currentVol = 1;

% Go through varargin
for ii=1:2:length(varargin)
	switch lower(varargin{ii})
		case 'slicedir'
			sliceDir = varargin{ii+1};
		case 'drawslicenbr'
			drawSliceNbr = varargin{ii+1};
		case 'drawrois'
			drawRois = varargin{ii+1};
		case 'drawoverlay'
			drawOverlay = varargin{ii+1};
		case 'drawcolorbar'
			drawColorbar = varargin{ii+1};
		case 'drawfilename'
			drawFileName = varargin{ii+1};
		case 'roiind'
			roi_ind = varargin{ii+1};
		case 'colormap'
			colmap = varargin{ii+1};
		case 'cols'
			cols = varargin{ii+1};
		case 'rows'
			rows = varargin{ii+1};
		case 'sliceind'
			sliceInd = varargin{ii+1};
		case 'vol'
			Dat.currentVol = varargin{ii+1};
		case 'clim'
			Clim_in = varargin{ii+1};
		case 'roitransp'
			roiTransp = varargin{ii+1};
		case 'showRoiEdges'
			showRoiEdges = varargin{ii+1};
	end
end

if isempty(ROI)
	drawRois = false;
	roi_ind = [];
end

if nargin < 3 || isempty(overlay)
	drawOverlay = false;
end

if ~iscell(DATA)
	DATA = {DATA};
end

if length(DATA)>1
	Dat.isDataMixed = true;
else
	Dat.isDataMixed = false;
end

% Store data sizes in all directions
if Dat.isDataMixed
  Dat.XSize = 0;
  Dat.YSize = 0;
  Dat.ZSize = length(DATA);
else
  Dat.XSize = size(DATA{1}.FTDATA,1);
  Dat.YSize = size(DATA{1}.FTDATA,2);
  Dat.ZSize = size(DATA{1}.FTDATA,3);
  Dat.VSize = size(DATA{1}.FTDATA,4);
end

if isempty(sliceInd)
	if sliceDir==1
		Dat.SliceInd = 1:Dat.XSize;
	elseif sliceDir==2
		Dat.SliceInd = 1:Dat.YSize;
	else
		Dat.SliceInd = 1:Dat.ZSize;
	end
else
	Dat.SliceInd = sliceInd;
end

if isempty(rows) && isempty(cols)
	rows = floor(sqrt(length(Dat.SliceInd)));
	cols = ceil(length(Dat.SliceInd)/rows);
elseif isempty(rows) 
	rows = ceil(length(Dat.SliceInd)/cols);
elseif isempty(cols)
	cols = ceil(length(Dat.SliceInd)/rows);
end

% Estimate space for the mosaic
if Dat.isDataMixed
	data_sz=size(DATA{Dat.SliceInd(1)}.FTDATA);
else
	if sliceDir==3
		data_sz=size(DATA{1}.FTDATA(:,:,Dat.SliceInd(1),Dat.currentVol));
	elseif sliceDir==2
		data_sz=size(squeeze(DATA{1}.FTDATA(:,Dat.SliceInd(1),:,Dat.currentVol)));
	else
		data_sz=size(squeeze(DATA{1}.FTDATA(Dat.SliceInd(1),:,:,Dat.currentVol)));
	end
end

% Determine figure size
scrsz = get(0,'screensize'); % Get current display resolution
scrsz(4) = scrsz(4)-45-70; % The 47px is for taskbar,
scrsz(3) = scrsz(3)-16;
gap=3;
fig_w = cols*data_sz(2)+(cols+1)*gap;
fig_h = rows*data_sz(1)+(rows+1)*gap;

if fig_w>scrsz(3) || fig_h>scrsz(4)
	
	aratio_w=[fig_w fig_h]./fig_w;
	aratio_h=[fig_w fig_h]./fig_h;
	for ii=1:2
		if fig_w>scrsz(3)
			tmp=aratio_w*scrsz(3);
			fig_w = tmp(1);
			fig_h = tmp(2);
		elseif fig_h>scrsz(4)
			tmp=aratio_h*scrsz(4);
			fig_w = tmp(1);
			fig_h = tmp(2);
		end
	end
end

%% Determine paper orientation
if fig_w>fig_h
	paperorient = 'landscape';
	
	%% Determine paperposition
	papersize = [29.7 21.0];
	tmp=[fig_w,fig_h]./fig_h;
	pap_h = tmp(2)*(papersize(2)-2);
	pap_w = tmp(1)*(papersize(2)-2);
	if pap_w>(papersize(1)-2)
		pap_h=((papersize(1)-2)/pap_w)*pap_h;
		pap_w = papersize(1)-2;
	end
	paperpos = [papersize(1)/2-pap_w/2 ...
		papersize(2)/2-pap_h/2 ...
		pap_w ...
		pap_h];
else
	paperorient = 'portrait';
	
	%% Determine paperposition
	papersize = [21.0 29.7];
	tmp=[fig_w,fig_h]./fig_w;
	pap_h = tmp(2)*(papersize(1)-2);
	pap_w = tmp(1)*(papersize(1)-2);
	if pap_h>(papersize(2)-2)
		pap_w=((papersize(2)-2)/pap_h)*pap_w;
		pap_h = papersize(2)-2;
	end
	paperpos = [papersize(1)/2-pap_w/2 ...
		papersize(2)/2-pap_h/2 ...
		pap_w ...
		pap_h];
end


%% Show wait dialog
[wait_h,txh]=aedes_calc_wait('Creating mosaic from slices...');

% Draw the mosaic figure
fh = figure('position',[scrsz(3)/2-fig_w/2+4 ...
	scrsz(4)/2-fig_h/2+45 fig_w fig_h],...
	'visible','on',...
	'inverthardcopy','off',...
	'renderer','opengl',...%'painters',...
	'numbertitle','off',...
	'name','Mosaic View',...
	'colormap',colmap,...
	'color','w',...
	'papertype','a4',...
	'paperpositionmode','manual',...
	'paperorientation',paperorient,...
	'paperunits','centimeters',...
	'paperposition',paperpos);
set(fh,'paperunits','inches')
if ( roiTransp==0 | not(drawRois)) && not(drawOverlay)
	set(fh,'renderer','painters');
end

% Set resolution information in the figure
if ispref('Aedes','ExportMosaicImageResolution')
	default_resolution = getpref('Aedes','ExportMosaicImageResolution');
else
	default_resolution = 300;
end
setappdata(fh,'Exportsetup',struct('Resolution',default_resolution))

%% Set header text
if ~Dat.isDataMixed
	hs = getappdata(fh,'PrintHeaderHeaderSpec');
	if isempty(hs)
		hs = struct('dateformat','none',...
			'string',[date,',  ',...
			strrep(strrep(DATA{1}.HDR.fpath,'\','\\'),'_','\_'),...
			strrep(strrep(DATA{1}.HDR.fname,'\','\\'),'_','\_')],...
			'fontname','Times',...
			'fontsize',12,...          % in points
			'fontweight','normal',...
			'fontangle','normal',...
			'margin',15);            % in points
	end
	setappdata(fh,'PrintHeaderHeaderSpec',hs)
end


% Look for custom aspect ratio
if ispref('Aedes','ExportMosaicUseCustomAspectRatio')
	data_aspectratio = getpref('Aedes','ExportMosaicUseCustomAspectRatio');
else
	data_aspectratio = [1 1 1];
end

%% images and axes
count=1;
gap_w=gap/fig_w;
gap_h=gap/fig_h;
for ii=1:rows
	for kk=1:cols
		ax=axes('parent',fh,...
			'units','normal',...
			'position',[(kk-1)*((1-(cols+1)*gap_w)/cols)+kk*gap_w ...
			(rows-ii)*((1-(rows+1)*gap_h)/rows)+(rows+1-ii)*gap_h ...
			(1-(cols+1)*gap_w)/cols (1-(rows+1)*gap_h)/rows],...
			'visible','off',...
			'ydir','reverse',...
			'xtick',[],...
			'ytick',[],...
			'xticklabel',[],...
			'yticklabel',[],...
			'DataAspectRatio',data_aspectratio,...
			'PlotBoxAspectRatio',[data_sz(2) data_sz(1) 1],...
			'PlotBoxAspectRatioMode','manual');
		if strcmpi(get(fh,'renderer'),'opengl')
			roi_ax = axes('parent',fh,...
				'units','normal',...
				'position',[(kk-1)*((1-(cols+1)*gap_w)/cols)+kk*gap_w ...
				(rows-ii)*((1-(rows+1)*gap_h)/rows)+(rows+1-ii)*gap_h ...
				(1-(cols+1)*gap_w)/cols (1-(rows+1)*gap_h)/rows],...
				'visible','off',...
				'ydir','reverse',...
				'xtick',[],...
				'ytick',[],...
				'xticklabel',[],...
				'yticklabel',[],...
				'DataAspectRatio',data_aspectratio,...
				'PlotBoxAspectRatio',[data_sz(2) data_sz(1) 1],...
				'PlotBoxAspectRatioMode','manual');
			overlay_ax = axes('parent',fh,...
				'units','normal',...
				'position',[(kk-1)*((1-(cols+1)*gap_w)/cols)+kk*gap_w ...
				(rows-ii)*((1-(rows+1)*gap_h)/rows)+(rows+1-ii)*gap_h ...
				(1-(cols+1)*gap_w)/cols (1-(rows+1)*gap_h)/rows],...
				'visible','off',...
				'ydir','reverse',...
				'xtick',[],...
				'ytick',[],...
				'xticklabel',[],...
				'yticklabel',[],...
				'DataAspectRatio',data_aspectratio,...
				'PlotBoxAspectRatio',[data_sz(2) data_sz(1) 1],...
				'PlotBoxAspectRatioMode','manual');
		end
		if count<=length(Dat.SliceInd)
			
			if Dat.isDataMixed
				imdata=DATA{Dat.SliceInd(count)}.FTDATA;
			else
				if sliceDir==3
					imdata=DATA{1}.FTDATA(:,:,Dat.SliceInd(count),Dat.currentVol);
					if drawOverlay
						if overlay.isOverlayRGB
							overlay_data = squeeze(overlay.ImageOverlay(:,:,Dat.SliceInd(count),:));
						else
							overlay_data = overlay.ImageOverlay(:,:,Dat.SliceInd(count),Dat.currentVol);
						end
					end
				elseif sliceDir==2
					imdata=squeeze(DATA{1}.FTDATA(:,Dat.SliceInd(count),:,Dat.currentVol));
					if drawOverlay
						if overlay.isOverlayRGB
							overlay_data = squeeze(overlay.ImageOverlay(:,Dat.SliceInd(count),:,:));
						else
							overlay_data = squeeze(overlay.ImageOverlay(:,Dat.SliceInd(count),:,Dat.currentVol));
						end
					end
				else
					imdata=squeeze(DATA{1}.FTDATA(Dat.SliceInd(count),:,:,Dat.currentVol));
					if drawOverlay
						if overlay.isOverlayRGB
							overlay_data = squeeze(overlay.ImageOverlay(Dat.SliceInd(count),:,:,:));
						else
							overlay_data = squeeze(overlay.ImageOverlay(Dat.SliceInd(count),:,:,Dat.currentVol));
						end
					end
				end
			end
			
			
			
			%% Plot image
			h=image('parent',ax,...
				'cdata',imdata,...
				'cdatamapping','scaled');
			axis(ax,'image')
			% Plot overlay
			if drawOverlay
				if ~overlay.isOverlayRGB
					ov_cmap = overlay.ImageOverlayCmap;
					ov_clim = round(((overlay.ImageOverlayClim-overlay.ImOverlayMin)*256)./...
						(overlay.ImOverlayMax-overlay.ImOverlayMin));
					if ov_clim(1)==ov_clim(2)
						if ov_clim(1)==0
							ov_clim(2)=1;
						elseif ov_clim(1)==256
							ov_clim(1)=255;
						end
					end
					ov_thold = round(((overlay.ImageOverlayThold-overlay.ImOverlayMin)*256)/...
						(overlay.ImOverlayMax-overlay.ImOverlayMin));
				end
				ov_alpha_val = overlay.ImageOverlayAlpha;
				
				if ~overlay.isOverlayRGB
					% Convert indexed image to RGB image
					%slice1_ind=overlay.ImageOverlay(:,:,Dat.Slices(1),Dat.CurrentVol);
					overlay_ind=double(overlay_data);
					
					% Get thresholded alpha indices
					if overlay.ImageOverlayTholdDirPos==1
						overlay_alpha_th = overlay_ind<ov_thold;
					else
						overlay_alpha_th = overlay_ind>ov_thold;
					end
					
					% Get clim alpha indices
					overlay_alpha_clim = ( overlay_ind>=ov_clim(1) & overlay_ind<=ov_clim(2) );
					
					overlay_ind(overlay_ind<ov_clim(1))=ov_clim(1);
					overlay_ind(overlay_ind>ov_clim(2))=ov_clim(2);
					
					overlay_ind=ceil((overlay_ind-ov_clim(1))./diff(ov_clim)*255+1);
					
					sz = size(overlay_ind);
					overlay_im = zeros([sz(1) sz(2) 3],'single');
					overlay_im(:,:,1)= reshape(ov_cmap(overlay_ind,1),sz);
					overlay_im(:,:,2)= reshape(ov_cmap(overlay_ind,2),sz);
					overlay_im(:,:,3)= reshape(ov_cmap(overlay_ind,3),sz);
					
					overlay_alpha = zeros(size(overlay_ind));
					overlay_alpha(overlay_alpha_clim) = ov_alpha_val;
					overlay_alpha(overlay_alpha_th) = 0;
				else
					overlay_im = overlay_data;
					overlay_alpha = ov_alpha_val;
				end
				
				
				
				h=image('parent',overlay_ax,...
					'cdata',overlay_im,...
					'cdatamapping','scaled',...
					'AlphaDataMapping','none',...
					'AlphaData',overlay_alpha,...
					'visible','on');
				
			end
			
			%% Check Clim
			if isempty(Clim_in)
				clim = [min(min(imdata)) max(max(imdata))];
				set(ax,'clim',clim)
			elseif size(Clim_in,1)>1
				set(ax,'clim',Clim_in(Dat.SliceInd(count),:))
			else
				clim = Clim_in;
				set(ax,'clim',clim)
			end
			
			% Draw ROIs
			if drawRois
				for tt=1:length(roi_ind)
					if Dat.isDataMixed
						roidata = ROI(roi_ind(tt)).voxels{Dat.SliceInd(count)};
					else
						if sliceDir==3
							roidata = ROI(roi_ind(tt)).voxels{1}(:,:,Dat.SliceInd(count),Dat.currentVol);
						elseif sliceDir==2
							roidata = squeeze(ROI(roi_ind(tt)).voxels{1}(:,Dat.SliceInd(count),:,Dat.currentVol));
						else
							roidata = squeeze(ROI(roi_ind(tt)).voxels{1}(Dat.SliceInd(count),:,:,Dat.currentVol));
						end
						if isempty(find(roidata))
							continue;
						end
					end
					
					if strcmpi(get(fh,'renderer'),'opengl')
						roidata=uint8(roidata);
						cdata = zeros([size(roidata) 3],'uint8');
						if ROI(roi_ind(tt)).color(1)~=0
							cdata(:,:,1) = roidata*ROI(roi_ind(tt)).color(1);
						end
						if ROI(roi_ind(tt)).color(2)~=0
							cdata(:,:,2) = roidata*ROI(roi_ind(tt)).color(2);
						end
						if ROI(roi_ind(tt)).color(3)~=0
							cdata(:,:,3) = roidata*ROI(roi_ind(tt)).color(3);
						end
						alphadata=double(roidata)*roiTransp;
						
						h=image('parent',roi_ax,'cdata',cdata,...
							'AlphaDataMapping','none',...
							'cdatamapping','scaled',...
							'AlphaData',alphadata,...
							'visible','on');
					end
					
					if showRoiEdges
						B=bwboundaries(roidata,4,'holes');
						for jj=1:length(B)
							boundary = B{jj};
							line('parent',ax,...
								'xdata',boundary(:,2),'ydata',boundary(:,1),...
								'color',ROI(roi_ind(tt)).color./255,...
								'tag','roiedge',...
								'linewidth',0.5,...
								'linestyle','-',...
								'hittest','off');
						end % for jj=1:length(B)
					end
				end % for tt=1:length(roi_ind)
			end % if drawRois
			
			% Draw slice number and filename
			if drawSliceNbr || drawFileName
				if (rows*cols)<=50
					fontsz=8;
				elseif (rows*cols)>50 & (rows*cols)<100
					fontsz=7;
				else
					fontsz=6;
				end
				if drawSliceNbr && drawFileName
					if Dat.isDataMixed
						fname = DATA{Dat.SliceInd(count)}.HDR.fname;
						fpath = DATA{Dat.SliceInd(count)}.HDR.fpath;
					else
						fname=DATA{1}.HDR.fname;
						fpath = DATA{1}.HDR.fpath;
					end
					if strcmp(fname,'fid')
						[fp,fn,fe]=fileparts(fpath(1:end-1));
						fname = [fn,fe];
					end
					slicetxt = sprintf('%03d: %s',Dat.SliceInd(count),fname);
				elseif drawFileName
					if Dat.isDataMixed
						fname = DATA{Dat.SliceInd(count)}.HDR.fname;
						fpath = DATA{Dat.SliceInd(count)}.HDR.fpath;
					else
						fname=DATA{1}.HDR.fname;
						fpath = DATA{1}.HDR.fpath;
					end
					if strcmp(fname,'fid')
						[fp,fn,fe]=fileparts(fpath(1:end-1));
						fname = [fn,fe];
					end
					slicetxt = fname;
				elseif drawSliceNbr
					slicetxt = sprintf('%03d',Dat.SliceInd(count));
				end
				tx_slice=text('parent',ax,...
					'units','normal',...
					'position',[0.01 0.99],...
					'verticalalign','top',...
					'horizontalalign','left',...
					'interpreter','none',...
					'string',slicetxt,...
					'backgroundcolor','w',...
					'clipping','off',...
					'fontsize',fontsz);
			end
			
			count=count+1;
		end % if count<=length(Dat.
	end % for kk=1:cols
end % for ii=1:rows


delete(wait_h)
