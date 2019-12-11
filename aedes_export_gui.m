function aedes_export_gui(DATA,ROI,Clim_in,colmap_in,roiTransp,showRoiEdges,overlay,currentVol)
% AEDES_EXPORT_GUI - Graphical user interface for exporting mri images
%   
%
% Synopsis:
%
% Description:
%
% Examples:
%
% See also:
%       AEDES

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


% Public variables
H=[];
Dat=[];

if nargin==0
  error('Too few input arguments!')
elseif nargin==1
  ROI = [];
else
  % Check ROI structure
  if ~isempty(ROI) && ( ~isstruct(ROI) || ~isfield(ROI,'voxels') )
    error('Invalid ROI structure!')
  end
end

if ~exist('Clim_in','var')
  Clim_in = [];
end
if ~exist('colmap_in','var')
  colmap_in = [];
end
if ~exist('overlay','var')
  overlay = [];
end
if ~exist('showRoiEdges','var')
  showRoiEdges = false;
end
if ~exist('roiTransp','var')
  roiTransp = 1;
end
if ~exist('currentVol','var')
  currentVol = 1;
end
% $$$ if ~exist('flip_in','var')
% $$$   flip_in = [];
% $$$ end
% $$$ if ~exist('rotation_in','var')
% $$$   rotation_in = [];
% $$$ end

% Check DATA structure
Dat.currentVol = currentVol;
showError = false;
Dat.isDataMixed = false;
if isstruct(DATA)
  if ~isfield(DATA,'FTDATA')
    showError = true;
  else
    DATA = {DATA};
  end
elseif iscell(DATA)
  if ~isstruct(DATA{1})
    showError = true;
  else
    if length(DATA)>1
      Dat.isDataMixed = true;
    end
  end
else
  showError = true;
end

if showError
  error('First input argument must be a valid DATA-form structure!')
  return
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

% Draw GUI figure and uicontrols
H=l_DrawGui;


%%%%%%%%%%%%%%%%%%%%%%%%
% Draw GUI
%%%%%%%%%%%%%%%%%%%%%%%%
function H=l_DrawGui

  %% Load default font and colors
  GD=aedes_gui_defaults;
  
	
  fig_w = 350;
  fig_h = 585;
	fig_location = aedes_dialoglocation([fig_w,fig_h]);
	fig_pos = [fig_location(1) fig_location(2) fig_w fig_h];
	
  %% Main Figure ----------------------------
  H.MAINFIG = figure('Units','Pixel', ...
                     'position',fig_pos,...
                     'Name','Export Dialog', ...
                     'Numbertitle','off', ...
                     'Tag','aedes_export_fig', ...
                     'Color',GD.col.mainfig, ...
                     'Toolbar','none', ...
                     'Menubar','none', ...
                     'DockControls','off',...
                     'renderer','painters',...
                     'KeyPressFcn','',...
                     'Closereq',@l_Close,...
                     'Handlevisibility','off',...
                     'windowstyle','normal');
	if ~GD.HG2graphics
		set(H.MAINFIG,'DoubleBuffer','on')
	end
  set(H.MAINFIG,'resize','off')
  
  %% Uipanels --------------------------------
  
  % Export options
  H.EXPORT_TYPE = uipanel('parent',H.MAINFIG,...
                          'units','pixel',...
                          'position',[10 fig_h-10-210 330 210],...
                          'title','Export Type and Path',...
						  'backgroundcolor',GD.col.frame);
  H.FILE_FORMAT_TX = uicontrol('parent',H.EXPORT_TYPE,...
                               'units','pixel',...
                               'position',[10 170 80 15],...
                               'string','File Format:',...
                               'horizontalalign','left',...
                               'style','text',...
							   'backgroundcolor',GD.col.frame);
  tmp=get(H.FILE_FORMAT_TX,'position');
  H.FILE_FORMAT_POPUP = uicontrol('parent',H.EXPORT_TYPE,...
                                  'units','pixel',...
                                  'position',[tmp(1)+tmp(3) tmp(2) 225 18],...
                                  'style','popup',...
                                  'string',...
                                  {'JPG - JPEG Files',...
                      'TIFF - Tagged Image File Format',...
                      'PNG - Portable Network Graphics',...
                      'BMP - Windows Bitmap',...
                      'EPS - Encapsulated PostScript'},...
                                  'backgroundcolor','w',...
                                  'callback',@l_ChangeFileFormat);
  H.FILETYPE_OPTIONS = uipanel('parent',H.EXPORT_TYPE,...
                               'units','pixel',...
                               'position',[tmp(1)+tmp(3) tmp(2)-80 225 75],...
                               'title','File Format Options',...
							   'backgroundcolor',GD.col.frame);
  H.RESOLUTION_TX = uicontrol('parent',H.FILETYPE_OPTIONS,...
                              'units','pixel',...
                              'position',[30 35 80 15],...
                              'string','Resolution (dpi):',...
                              'horizontalalign','left',...
                              'style','text',...
							  'backgroundcolor',GD.col.frame);
  tmp=get(H.RESOLUTION_TX,'position');
  H.RESOLUTION_EDIT = uicontrol('parent',H.FILETYPE_OPTIONS,...
                                'units','pixel',...
                                'position',[tmp(1)+tmp(3)+10 tmp(2) 50 18],...
                                'backgroundcolor','w',...
                                'horizontalalign','center',...
                                'style','edit',...
                                'string','96',...
                                'userdata','96',...
                                'callback',@l_ChangeResolutionQuality);
  H.COMPRESSION_CHBOX = uicontrol('parent',H.FILETYPE_OPTIONS,...
                                  'units','pixel',...
                                  'position',[tmp(1) tmp(2)-25 150 18],...
                                  'horizontalalign','center',...
                                  'style','checkbox',...
                                  'string','Use Compression',...
                                  'value',1,...
                                  'visible','off',...
								  'backgroundcolor',GD.col.frame);
  H.QUALITY_TX = uicontrol('parent',H.FILETYPE_OPTIONS,...
                           'units','pixel',...
                           'position',[tmp(1) tmp(2)-25 75 18],...
                           'string','Quality (%):',...
                           'horizontalalign','left',...
                           'style','text',...
						   'backgroundcolor',GD.col.frame);
  tmp2=get(H.QUALITY_TX,'position');
  tmp3=get(H.RESOLUTION_EDIT,'position');
  H.QUALITY_EDIT = uicontrol('parent',H.FILETYPE_OPTIONS,...
                             'units','pixel',...
                             'position',[tmp3(1) tmp2(2) 50 18],...
                             'horizontalalign','center',...
                             'style','edit',...
                             'backgroundcolor','w',...
                             'string','80',...
                             'value',1,...
                             'visible','on',...
                             'userdata','80',...
                             'callback',@l_ChangeResolutionQuality);
  tmp2=get(H.FILETYPE_OPTIONS,'position');
  tmp=get(H.FILE_FORMAT_TX,'position');
  H.FILEPREFIX_TX = uicontrol('parent',H.EXPORT_TYPE,...
                              'units','pixel',...
                              'position',[tmp(1) tmp2(2)-tmp(4)-20 tmp(3:4)],...
                              'string','File Prefix:',...
                              'horizontalalign','left',...
                              'style','text',...
							  'backgroundcolor',GD.col.frame);
  tmp=get(H.FILEPREFIX_TX,'position');
  H.FILEPREFIX_EDIT = uicontrol('parent',H.EXPORT_TYPE,...
                                'units','pixel',...
                                'position',[tmp(1)+tmp(3) tmp(2) 225 18],...
                                'backgroundcolor','w',...
                                'horizontalalign','left',...
                                'style','edit',...
                                'string','mrimage',...
                                'keypressfcn',@l_FilePrefixEditKeyPress);
  % This dummy uicontrol is used to safely shift focus
  % between uicontrols. Focus shift is used as a workaround to a Matlab
  % keypressfcn bug with editboxes...
  H.DUMMYUICH = uicontrol('parent',H.EXPORT_TYPE,...
                          'units','pixel',...
                          'position',[-20 0 2 2],...
                          'style','pushbutton',...
                          'visible','on');
  H.OUTDIR_TX = uicontrol('parent',H.EXPORT_TYPE,...
                          'units','pixel',...
                          'position',[tmp(1) tmp(2)-25 90 15],...
                          'string','Output directory:',...
                          'horizontalalign','left',...
                          'style','text',...
						  'backgroundcolor',GD.col.frame);
  tmp=get(H.OUTDIR_TX,'position');
  H.OUTDIR_EDIT = uicontrol('parent',H.EXPORT_TYPE,...
                            'units','pixel',...
                            'position',[tmp(1) tmp(2)-18 280 18],...
                            'string','',...
                            'backgroundcolor','w',...
                            'horizontalalign','left',...
                            'style','edit',...
                            'callback',@l_CheckOutputDir);
  % Set default output directory
  try
    filepath=getpref('Aedes','ExportFileDir');
    set(H.OUTDIR_EDIT,'string',filepath)
  catch
    set(H.OUTDIR_EDIT,'string',[pwd,filesep])
  end
  tmp=get(H.OUTDIR_EDIT,'position');
  H.BROWSE_BTN = uicontrol('parent',H.EXPORT_TYPE,...
                           'units','pixel',...
                           'position',[tmp(1)+tmp(3)+2 tmp(2) 25 18],...
                           'string','...',...
                           'style','pushbutton',...
                           'callback',@l_CheckOutputDir);
  
  % Separate files/mosaic uipanel
  tmp=get(H.EXPORT_TYPE,'position');
  H.FILES_UIPANEL = uipanel('parent',H.MAINFIG,...
                            'units','pixel',...
                            'position',[tmp(1) tmp(2)-10-120 160 120],...
                            'title','Files',...
							'backgroundcolor',GD.col.frame);
  set(H.MAINFIG,'color',get(H.FILES_UIPANEL,'backgroundcolor'));
  H.RADIOGRP_FILES = uibuttongroup('parent',H.FILES_UIPANEL,...
                                   'units','pixel',...
                                   'position',[5 5 150 100],...
                                   'bordertype','none',...
                                   'selectionchangefcn',@l_FilesSelectionCB,...
								   'backgroundcolor',GD.col.frame);
  H.SEP_FILES = uicontrol('parent',H.RADIOGRP_FILES,...
                          'units','pixel',...                        
                          'position',[10 75 130 25],...
                          'style','radio',...
                          'string','Separate Files',...
                          'value',1,...
						  'backgroundcolor',GD.col.frame);
  tmp=get(H.SEP_FILES,'position');
  H.MOSAIC_FILES = uicontrol('parent',H.RADIOGRP_FILES,...
                             'units','pixel',...                        
                             'position',[tmp(1) tmp(2)-tmp(4) tmp(3:4)],...
                             'style','radio',...
                             'string','Mosaic',...
                             'value',0,...
							 'backgroundcolor',GD.col.frame);
  tmp=get(H.MOSAIC_FILES,'position');
  H.ROWS_TX = uicontrol('parent',H.RADIOGRP_FILES,...
                        'units','pixel',...
                        'position',[40 tmp(2)-20 50 15],...
                        'style','text',...
                        'horizontalalign','left',...
                        'string','Rows:',...
                        'enable','off',...
						'backgroundcolor',GD.col.frame);
  H.ROWS_EDIT = uicontrol('parent',H.RADIOGRP_FILES,...
                          'units','pixel',...
                          'position',[90 tmp(2)-20 50 18],...
                          'style','edit',...
                          'backgroundcolor','w',...
                          'string','',...
                          'enable','off',...
                          'callback',@l_CheckMosaic);
  tmp=get(H.ROWS_TX,'position');
  H.COLS_TX = uicontrol('parent',H.RADIOGRP_FILES,...
                        'units','pixel',...
                        'position',[tmp(1) tmp(2)-25 tmp(3:4)],...
                        'style','text',...
                        'horizontalalign','left',...
                        'string','Columns:',...
                        'enable','off',...
						'backgroundcolor',GD.col.frame);
  H.COLS_EDIT = uicontrol('parent',H.RADIOGRP_FILES,...
                          'units','pixel',...
                          'position',[90 tmp(2)-25 50 18],...
                          'style','edit',...
                          'backgroundcolor','w',...
                          'string','',...
                          'enable','off',...
                          'callback',@l_CheckMosaic);
  
  % Set default values for mosaic rows/columns
  rows = floor(sqrt(Dat.ZSize));
  cols = ceil(Dat.ZSize/rows);
  set(H.ROWS_EDIT,'string',num2str(rows),...
                  'userdata',num2str(rows))
  set(H.COLS_EDIT,'string',num2str(cols),...
                  'userdata',num2str(cols))
  
  % X-dir, Y-dir, Z-dir
  tmp=get(H.FILES_UIPANEL,'position');
  H.DIR_UIPANEL = uipanel('parent',H.MAINFIG,...
                          'units','pixel',...
                          'position',[180 tmp(2) tmp(3) tmp(4)],...
                          'title','Slice Direction',...
						  'backgroundcolor',GD.col.frame);
  H.RADIOGRP_DIR = uibuttongroup('parent',H.DIR_UIPANEL,...
                                 'units','pixel',...
                                 'position',[5 5 150 100],...
                                 'bordertype','none',...
                                 'selectionchangefcn',@l_DirChangeFcn,...
								 'backgroundcolor',GD.col.frame);
  H.XDIR_RADIO = uicontrol('parent',H.RADIOGRP_DIR,...
                           'units','pixel',...                        
                           'position',[10 65 130 25],...
                           'style','radio',...
                           'string','X-Direction',...
                           'value',0,...
						   'backgroundcolor',GD.col.frame);
  tmp=get(H.XDIR_RADIO,'position');
  H.YDIR_RADIO = uicontrol('parent',H.RADIOGRP_DIR,...
                           'units','pixel',...                        
                           'position',[tmp(1) tmp(2)-tmp(4) tmp(3:4)],...
                           'style','radio',...
                           'string','Y-Direction',...
                           'value',0,...
						   'backgroundcolor',GD.col.frame);
  tmp=get(H.YDIR_RADIO,'position');
  H.ZDIR_RADIO = uicontrol('parent',H.RADIOGRP_DIR,...
                           'units','pixel',...                        
                           'position',[tmp(1) tmp(2)-tmp(4) tmp(3:4)],...
                           'style','radio',...
                           'string','Z-Direction',...
                           'value',1,...
						   'backgroundcolor',GD.col.frame);
  if Dat.isDataMixed
    % Disable Y- and Z-radio buttons for mixed data
    set([H.YDIR_RADIO,H.XDIR_RADIO],'enable','off')
  end
  
  % Slice selection
  tmp=get(H.FILES_UIPANEL,'position');
  H.SLICE_UIPANEL = uipanel('parent',H.MAINFIG,...
                            'units','pixel',...
                            'position',[10 tmp(2)-180-10 tmp(3) 180],...
                            'title','Slice Selection',...
							'backgroundcolor',GD.col.frame);
  H.RADIOGRP_SLICE = uibuttongroup('parent',H.SLICE_UIPANEL,...
                                   'units','pixel',...
                                   'position',[5 5 150 170],...
                                   'bordertype','none',...
                                   'selectionchangefcn',@l_SliceSelection,...
								   'backgroundcolor',GD.col.frame);
  H.ALL_SLICES = uicontrol('parent',H.RADIOGRP_SLICE,...
                           'units','pixel',...                        
                           'position',[10 135 130 25],...
                           'style','radio',...
                           'string','All Slices',...
                           'value',1,...
						   'backgroundcolor',GD.col.frame);
  tmp=get(H.ALL_SLICES,'position');
  H.CUSTOM_SLICES = uicontrol('parent',H.RADIOGRP_SLICE,...
                              'units','pixel',...                        
                              'position',[tmp(1) tmp(2)-tmp(4) tmp(3:4)],...
                              'style','radio',...
                              'string','Custom Slices',...
                              'value',0,...
							  'backgroundcolor',GD.col.frame);
  tmp=get(H.CUSTOM_SLICES,'position');
  H.CUSTOM_EDIT = uicontrol('parent',H.RADIOGRP_SLICE,...
                            'units','pixel',...                        
                            'position',[tmp(1)+20 tmp(2)-tmp(4)+5 110 18],...
                            'style','edit',...
                            'backgroundcolor','w',...
                            'horizontalalign','left',...
                            'string','1:end',...
                            'userdata','1:end',...
                            'enable','off',...
                            'callback',@l_CheckCustomSliceSelection);
  
  % Additional Options
  tmp=get(H.FILES_UIPANEL,'position');
  H.OPTIONS_UIPANEL = uipanel('parent',H.MAINFIG,...
                              'units','pixel',...
                              'position',[180 tmp(2)-10-180 160 180],...
                              'title','Additional Options',...
							  'backgroundcolor',GD.col.frame);
  if isempty(ROI)
    roi_enable = 'off';
  else
    roi_enable = 'on';
  end
  H.DRAW_ROIS = uicontrol('parent',H.OPTIONS_UIPANEL,...
                          'units','pixel',...
                          'position',[10 140 130 25],...
                          'style','checkbox',...
                          'string','Draw ROIs',...
                          'value',0,...
                          'enable',roi_enable,...
                          'callback',@l_DrawRoisCheckBoxCB,...
						  'backgroundcolor',GD.col.frame);
  if isempty(ROI)
    % Disable "Draw ROIs" checkbox if ROIs are not available
    set(H.DRAW_ROIS,'value',0,'enable','off')
  end
  tmp=get(H.DRAW_ROIS,'position');
  H.ROI_LBOX = uicontrol('parent',H.OPTIONS_UIPANEL,...
                         'units','pixel',...
                         'position',[tmp(1)+20 tmp(2)-40 110 40],...
                         'style','listbox',...
                         'backgroundcolor','w',...
                         'min',1,'max',3,...
                         'string',{''},...
                         'enable','off');
  if ~isempty(ROI)
    Dat.RoiLabels = {ROI(:).label};
    set(H.ROI_LBOX,'string',Dat.RoiLabels,...
			'enable','on','value',[1:length(ROI)])
		set(H.DRAW_ROIS,'value',1)
  end
  tmp=get(H.DRAW_ROIS,'position');
  H.DRAW_SLICENBR = uicontrol('parent',H.OPTIONS_UIPANEL,...
                              'units','pixel',...
                              'position',[tmp(1) tmp(2)-60 130 20],...
                              'style','checkbox',...
                              'string','Draw Slice Nbr.',...
                              'value',0,...
							  'backgroundcolor',GD.col.frame);
  tmp=get(H.DRAW_SLICENBR,'position');
  H.DRAW_FILENAME = uicontrol('parent',H.OPTIONS_UIPANEL,...
                              'units','pixel',...
                              'position',[tmp(1) tmp(2)-23 130 20],...
                              'style','checkbox',...
                              'string','Draw File Name',...
                              'value',0,...
							  'backgroundcolor',GD.col.frame);
  tmp=get(H.DRAW_FILENAME,'position');
  H.DRAW_COLORBAR = uicontrol('parent',H.OPTIONS_UIPANEL,...
                              'units','pixel',...
                              'position',[tmp(1) tmp(2)-23 130 20],...
                              'style','checkbox',...
                              'string','Draw Colorbar',...
                              'value',0,...
							  'backgroundcolor',GD.col.frame);
  tmp=get(H.DRAW_COLORBAR,'position');
  H.DRAW_OVERLAY = uicontrol('parent',H.OPTIONS_UIPANEL,...
                              'units','pixel',...
                              'position',[tmp(1) tmp(2)-23 130 20],...
                              'style','checkbox',...
                              'string','Draw Overlay',...
                              'value',0,...
							  'backgroundcolor',GD.col.frame);
  if isempty(overlay)
		set(H.DRAW_OVERLAY,'enable','off')
	else
		set(H.DRAW_OVERLAY,'value',1)
  end
  tmp=get(H.OPTIONS_UIPANEL,'position');
  
  % Export and Cancel buttons
  H.EXPORT_BTN = uicontrol('parent',H.MAINFIG,...
                       'units','pixel',...
                       'position',[tmp(1) tmp(2)-30-5 77.5 30],...
                       'string','Export...',...
                       'style','pushbutton',...
                       'callback',@l_CheckOptions);
  tmp = get(H.EXPORT_BTN,'position');
  H.CANCEL_BTN = uicontrol('parent',H.MAINFIG,...
                           'units','pixel',...
                           'position',[tmp(1)+tmp(3)+5 tmp(2:4)],...
                           'string','Cancel',...
                           'style','pushbutton',...
                           'callback','close(gcbf)');
  
  % Initialize some internal variables
  Dat.SliceInd = 1:Dat.ZSize;
  
  
end % function l_DrawGui(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GUI callbacks
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_ChangeFileFormat(h,evd)

popupval = get(H.FILE_FORMAT_POPUP,'value');

if popupval==1
  set([H.QUALITY_TX,H.QUALITY_EDIT],'visible','on')
  set(H.COMPRESSION_CHBOX,'visible','off')
elseif popupval==2
  set([H.QUALITY_TX,H.QUALITY_EDIT],'visible','off')
  set(H.COMPRESSION_CHBOX,'visible','on')
else
  set([H.QUALITY_TX,H.QUALITY_EDIT,H.COMPRESSION_CHBOX],...
      'visible','off')
end
  

end % function l_ChangeFileFormat(h,

function l_FilePrefixEditKeyPress(h,evd)

% This ugly stabbing of Matlab doesn't work anymore correctly with newer
% versions so return immediately...
return  

% Switch focus between uicontrols (a Matlab bug workaround)
uicontrol(H.DUMMYUICH)
uicontrol(H.FILEPREFIX_EDIT)

% Call drawnow twice to refresh the editbox (yet another Matlab bug
% workaround)
drawnow
drawnow

% Query the editbox string value
str=get(H.FILEPREFIX_EDIT,'string');
if isempty(str)
  return
end

% If the string contains forbidden characters, display error and remove
% the last character
ind=ismember(str,'\/:*?"<>|');
if any(ind)
  warndlg({'A file name cannot contain any of the following characters:',...
           '','\/:*?"<>|'},'Forbidden character entered.','modal')
  str(ind)=[];
  set(H.FILEPREFIX_EDIT,'string',str)
  uicontrol(H.DUMMYUICH)
  uicontrol(H.FILEPREFIX_EDIT)
  drawnow
  drawnow
end


end % function l_FilePrefixEditKeyPress(h,

function l_DrawRoisCheckBoxCB(h,evd)

if get(H.DRAW_ROIS,'value')
  set(H.ROI_LBOX,'enable','on')
else
  set(H.ROI_LBOX,'enable','off')
end

end % function l_DrawRoisCheckBoxCB(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_FilesSelectionCB(h,evd)

if ~isempty(evd)
  if evd.NewValue==H.MOSAIC_FILES
    child_h = findobj(H.EXPORT_TYPE,'type','uicontrol');
    %child_h=get(H.EXPORT_TYPE,'children');
    %child_h = child_h(child_h~=H.FILETYPE_OPTIONS);
    set(child_h,'enable','off')
    set([H.ROWS_TX,H.ROWS_EDIT,H.COLS_TX,H.COLS_EDIT],'enable','on')
    %set([H.FILE_FORMAT_TX,H.FILE_FORMAT_POPUP,...
    %    H.OUTDIR_TX,H.OUTDIR_EDIT,H.BROWSE_BTN],'enable','off')
  else
    child_h = findobj(H.EXPORT_TYPE,'type','uicontrol');
    %child_h=get(H.EXPORT_TYPE,'children');
    %child_h = child_h(child_h~=H.FILETYPE_OPTIONS);
    set(child_h,'enable','on')
    %set(get(H.EXPORT_TYPE,'children'),'enable','on')
    set([H.ROWS_TX,H.ROWS_EDIT,H.COLS_TX,H.COLS_EDIT],'enable','off')
    %set([H.FILE_FORMAT_TX,H.FILE_FORMAT_POPUP,...
    %    H.OUTDIR_TX,H.OUTDIR_EDIT,H.BROWSE_BTN],'enable','on')
  end
end

end % function l_FilesSelectionCB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_ChangeResolutionQuality(h,evd)

if h==H.RESOLUTION_EDIT
  res_val=floor(str2num(get(H.RESOLUTION_EDIT,'string')));
  if isempty(res_val) || ~isreal(res_val) || res_val<0 || res_val>1200
    hh=errordlg(['Invalid value for resolution! The value must be between ' ...
                 '0-1200 dpi.'],'Error!','modal');
    set(H.RESOLUTION_EDIT,'string',get(H.RESOLUTION_EDIT,'userdata'))
    return
  else
    set(H.RESOLUTION_EDIT,'string',num2str(res_val),...
                      'userdata',num2str(res_val))
  end
elseif h==H.QUALITY_EDIT
  qval = floor(str2num(get(H.QUALITY_EDIT,'string')));
  if isempty(qval) || ~isreal(qval) || qval<1 || qval>100
    hh=errordlg('Invalid value for quality! The value must be between 1-100%.',...
                'Error!','modal');
    set(H.QUALITY_EDIT,'string',get(H.QUALITY_EDIT,'userdata'))
    return
  else
    set(H.QUALITY_EDIT,'string',num2str(qval),...
                      'userdata',num2str(qval))
  end
end

end % function l_ChangeResolutionQuality(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_CheckMosaic(h,evd)

if ~isempty(h)

  % Check that the value is numeric
  val = str2num(get(h,'string'));
  if isempty(val) || ~isreal(val) || val<=0
    errordlg('Invalid row/column value!','Error. Invalid value!','modal');
    set(H.COLS_EDIT,'string',get(H.COLS_EDIT,'userdata'))
    set(H.ROWS_EDIT,'string',get(H.ROWS_EDIT,'userdata'))
    return
  end

  total_sz = length(Dat.SliceInd);
  
  if val>total_sz
    val = total_sz;
  end

  % Set values for mosaic rows/columns
  if h==H.ROWS_EDIT
    rows = val;
    cols = ceil(total_sz/rows);
    set(H.COLS_EDIT,'string',num2str(cols),...
                    'userdata',num2str(cols))
    set(H.ROWS_EDIT,'string',num2str(val))
  elseif h==H.COLS_EDIT
    cols = val;
    rows = ceil(total_sz/cols);
    set(H.ROWS_EDIT,'string',num2str(rows),...
                    'userdata',num2str(rows))
    set(H.COLS_EDIT,'string',num2str(val))
  end
  
else
  
% $$$   % Keep number of rows constant and change the number of columns
% $$$   total_sz = length(Dat.SliceInd);
% $$$   val = str2num(get(H.ROWS_EDIT,'string'));
% $$$   if val>total_sz
% $$$     val = total_sz;
% $$$   end
% $$$   rows = val;
% $$$   cols = ceil(total_sz/rows);
% $$$   set(H.COLS_EDIT,'string',num2str(cols),...
% $$$                   'userdata',num2str(cols))
% $$$   set(H.ROWS_EDIT,'string',num2str(val))
  
  % Keep number of columns constant and change the number of rows
  total_sz = length(Dat.SliceInd);
  val = str2num(get(H.COLS_EDIT,'string'));
  if val>total_sz
    val = total_sz;
  end
  cols = val;
  rows = ceil(total_sz/cols);
  set(H.ROWS_EDIT,'string',num2str(rows),...
                  'userdata',num2str(rows))
  set(H.COLS_EDIT,'string',num2str(val))
  
end

end % function l_CheckMosaic(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_CheckOutputDir(h,evd)

if h==H.BROWSE_BTN % Called from browse btn
  
  try
    filepath=getpref('Aedes','ExportFileDir');
  catch
    filepath='';
  end
  filepath = uigetdir(filepath);
  if isequal(filepath,0)
    return
  end
  filepath = [filepath,filesep];
  setpref('Aedes','ExportFileDir',filepath)
  set(H.OUTDIR_EDIT,'string',filepath)
  
elseif h==H.OUTDIR_EDIT % Called from outdir edit
  
  filepath=get(H.OUTDIR_EDIT,'string');
  
  % Check if the inputted directory exists
  if ~isdir(filepath)
    errordlg({'The export path','',...
              ['"',filepath,'"'],'','does not exist!'},...
             'Error! File path does not exist!','modal')
    return
  end
  
  if filepath(end)~=filesep
    filepath(end+1)=filesep;
    set(H.OUTDIR_EDIT,'string',filepath)
  end
  
end

end % function l_CheckOutputDir(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_CheckCustomSliceSelection(h,evd)

% Get selected direction
selobj=get(H.RADIOGRP_DIR,'selectedobject');
if selobj==H.XDIR_RADIO
  sz = Dat.XSize;
elseif selobj==H.YDIR_RADIO
  sz = Dat.YSize;
else
  sz = Dat.ZSize;
end

%% Get custom string
str = get(H.CUSTOM_EDIT,'string');
if isempty(str)
  hh=warndlg('The custom text field cannot be empty!',...
            'Invalid expression','modal');
  set(H.CUSTOM_EDIT,'string',get(H.CUSTOM_EDIT,'userdata'))
  return
end

%% Check that the string contains only valid characters
tmp_ind=ismember(str,'0123456789end:,*+-/');
if ~all(tmp_ind)
  hh=warndlg(['The custom text "' str '" contains invalid characters'],...
            'Invalid expression','modal');
  set(H.CUSTOM_EDIT,'string',get(H.CUSTOM_EDIT,'userdata'))
  return
end

% replace 'end' statement with ROI size
custom_str=strrep(str,'end','sz');
slice_ind=[];

% Try to evaluate the index string
try
  eval(['slice_ind=[' custom_str '];']);
catch
  hh=warndlg({['Could not evaluate expression "' custom_str '".'],...
             ['Custom string has to be a valid Matlab vector expression.']},...
            'Could not evaluate expression','modal');
  set(H.CUSTOM_EDIT,'string',get(H.CUSTOM_EDIT,'userdata'))
  return
end

% Sort indices and make sure that all indices are unique
%slice_ind=unique(slice_ind);

% Make sure that min>=1 and max<=sz
if max(slice_ind)>sz
  hh=warndlg({['Could not evaluate expression "' custom_str '".'],...
             'Maximum slice index exceeds the number of slices.'},...
            'Could not evaluate expression','modal');
  set(H.CUSTOM_EDIT,'string',get(H.CUSTOM_EDIT,'userdata'))
  return
end
if min(slice_ind)<1
  hh=warndlg({['Could not evaluate expression "' custom_str '".'],...
             ['Slice indices have to be larger than zero.']},...
            'Could not evaluate expression','modal');
  set(H.CUSTOM_EDIT,'string',get(H.CUSTOM_EDIT,'userdata'))
  return
end

Dat.SliceInd = slice_ind;
set(H.CUSTOM_EDIT,'userdata',str)

l_CheckMosaic([],[])

end % function l_CheckCustomSliceSelection(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_DirChangeFcn(h,evd)

set(H.CUSTOM_EDIT,'string','1:end')
if get(H.XDIR_RADIO,'value')==1
  Dat.SliceInd = 1:Dat.XSize;
elseif get(H.YDIR_RADIO,'value')==1
  Dat.SliceInd = 1:Dat.YSize;
else
  Dat.SliceInd = 1:Dat.ZSize;
end

l_CheckMosaic([],[])

end % function l_DirChangeFcn(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_SliceSelection(h,evd)

if h==H.RADIOGRP_SLICE
  if evd.NewValue==H.ALL_SLICES
    % Disable Custom slice editbox if All Slices is selected
    set(H.CUSTOM_EDIT,'enable','off')
    if get(H.XDIR_RADIO,'value')==1
      Dat.SliceInd = 1:Dat.XSize;
    elseif get(H.YDIR_RADIO,'value')==1
      Dat.SliceInd = 1:Dat.YSize;
    else
      Dat.SliceInd = 1:Dat.ZSize;
    end
  else
    set(H.CUSTOM_EDIT,'enable','on')
  end
end

l_CheckMosaic([],[])

end % function l_SliceSelection(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check Options and continue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_CheckOptions(h,evd)

% Check output dir
if get(H.SEP_FILES,'value')==1
  createMosaic = false;
  outdir = get(H.OUTDIR_EDIT,'string');
  if ~isdir(outdir)
    hh=errordlg({'The output directory',...
               ['"',outdir,'"'],'does not exist!'},'Error!','modal');
    return
  end
else
  outdir = get(H.OUTDIR_EDIT,'string');
  createMosaic = true;
end

if outdir(end)~=filesep
  outdir(end+1)=filesep;
end

% Check File Format Options
useCompression = logical(get(H.COMPRESSION_CHBOX,'value'));
resolution = get(H.RESOLUTION_EDIT,'string');
quality = str2num(get(H.QUALITY_EDIT,'string'));

% Get file prefix (for separate files)
filePrefix = get(H.FILEPREFIX_EDIT,'string');

% If the string contains forbidden characters, display error and return
ind=ismember(filePrefix,'\/:*?"<>|');
if any(ind)
  hh=warndlg({'A file name cannot contain any of the following characters:',...
	'','\/:*?"<>|'},'Forbidden character entered.','modal');
  uiwait(hh)
  uicontrol(H.DUMMYUICH)
  uicontrol(H.FILEPREFIX_EDIT)
  drawnow
  drawnow
  return
end

% Check direction
if get(H.XDIR_RADIO,'value')==1
  sliceDir = 1;
elseif get(H.YDIR_RADIO,'value')==1
  sliceDir = 2;
elseif get(H.ZDIR_RADIO,'value')==1
  sliceDir = 3;
end

% Check additional options
if get(H.DRAW_ROIS,'value')==1
  drawRois = true;
  roi_ind=get(H.ROI_LBOX,'value');
  if isempty(roi_ind)
    drawRois = false;
  end
else
  drawRois = false;
  roi_ind = [];
end
if get(H.DRAW_SLICENBR,'value')==1
  drawSliceNbr = true;
else
  drawSliceNbr = false;
end
if get(H.DRAW_FILENAME,'value')==1
  drawFileName = true;
else
  drawFileName = false;
end
if get(H.DRAW_COLORBAR,'value')==1
  drawColorbar = true;
else
  drawColorbar = false;
end
if get(H.DRAW_OVERLAY,'value')==1
  drawOverlay = true;
else
  drawOverlay = false;
end

% Check Colormap
if isempty(colmap_in)
  colmap = gray(256);
else
  colmap = colmap_in;
end

% % Check Clim
% if isempty(Clim_in)
%   colmap = gray(256);
% else
%   colmap = colmap_in;
% end



% Draw Mosaic figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if createMosaic
  
	%% Show wait dialog
	[wait_h,txh]=aedes_calc_wait('Creating mosaic from slices...');
	
  % Estimate space for the mosaic
  cols=str2num(get(H.COLS_EDIT,'string'));
  rows=str2num(get(H.ROWS_EDIT,'string'));
	
	aedes_createmosaic(DATA,ROI,overlay,...
		'slicedir',sliceDir,...
		'drawslicenbr',drawSliceNbr,...
		'drawrois',drawRois,...
		'drawfilename',drawFileName,...
		'drawcolorbar',drawColorbar,...
		'drawoverlay',drawOverlay,...
		'colormap',colmap,...
		'cols',cols,...
		'rows',rows,...
		'roiind',roi_ind,...
		'sliceind',Dat.SliceInd,...
		'clim',Clim_in,...
		'vol',Dat.currentVol,...
		'roitransp',roiTransp,...
		'showRoiEdges',showRoiEdges);
  
  
  delete(wait_h)
  
% Draw separate files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  
  % Get file extension
  fext_val = get(H.FILE_FORMAT_POPUP,'value');
  if fext_val==1
    fext = 'jpg';
    print_type = 'jpeg';
    quality = sprintf('%02d',quality);
  elseif fext_val==2
    fext = 'tif';
    quality = '';
    if useCompression
      print_type = 'tiff';
    else
      print_type = 'tiffnocompression';
    end
  elseif fext_val==3
    fext = 'png';
    print_type = 'png';
    quality = '';
  elseif fext_val==4
    fext = 'bmp';
    print_type = 'bmp';
    quality = '';
  elseif fext_val==5
    fext = 'eps';
    print_type = 'epsc2';
    quality = '';
  end

  % Create file names
  fnames={};
  for ii=1:length(Dat.SliceInd)
    fnames{ii}=sprintf('%s_%03d.%s',filePrefix,Dat.SliceInd(ii),fext);
  end

  % Get file names from the output directory
  tmp=dir(outdir);
  outdir_files={tmp(~[tmp(:).isdir]).name};

  if ~isunix % Check file names, in windows ignore case
    ind=ismember(lower(fnames),lower(outdir_files));
  else
    ind=ismember(fnames,outdir_files);
  end
  if any(ind)
    % Warn if some files are about to be overwritten
    overwrited_files = {fnames{ind}};

    % Limit the file list to 20 files
    if length(overwrited_files)>20
      overwrited_files = {overwrited_files{1:20}};
      overwrited_files{end+1}='...';
    end

    resp=questdlg({['The following files already exist in the output directory' ...
      ' "',outdir,'"'],'','(NOTE: max. 20 files shown here)','',...
      overwrited_files{:},'',...
      'Do you want to overwrite these files?'},...
      'Overwrite Existing Files?','Overwrite','Abort','Abort');
    if strcmpi(resp,'Abort')
      return
    end
  end

  % Initialize waitbar
  wbar_h = aedes_wbar(0,sprintf('Creating %s-image %d/%d to folder\n"%s"',...
    upper(fext),0,length(fnames),outdir));

  % Print images to files
  for ii=1:length(Dat.SliceInd)

    % Initialize waitbar
    aedes_wbar(ii/length(fnames),wbar_h,sprintf('Creating %s-image %d/%d to folder\n"%s"',...
      upper(fext),ii,length(fnames),...
      outdir));


    % Image data
    if Dat.isDataMixed
      imdata = DATA{Dat.SliceInd(ii)}.FTDATA;
    else
      if sliceDir==3
        imdata = DATA{1}.FTDATA(:,:,Dat.SliceInd(ii),Dat.currentVol);
        if drawOverlay
          if overlay.isOverlayRGB
            overlay_data = squeeze(overlay.ImageOverlay(:,:,Dat.SliceInd(ii),:));
          else
            overlay_data = overlay.ImageOverlay(:,:,Dat.SliceInd(ii),Dat.currentVol);
          end
        end
      elseif sliceDir==2
        imdata = squeeze(DATA{1}.FTDATA(:,Dat.SliceInd(ii),:,Dat.currentVol));
        if drawOverlay
          if overlay.isOverlayRGB
            overlay_data = squeeze(overlay.ImageOverlay(:,Dat.SliceInd(ii),:,:));
          else
            overlay_data = squeeze(overlay.ImageOverlay(:,Dat.SliceInd(ii),:,Dat.currentVol));
          end
        end
      else
        imdata = squeeze(DATA{1}.FTDATA(Dat.SliceInd(ii),:,:,Dat.currentVol));
        if drawOverlay
          if overlay.isOverlayRGB
            overlay_data = squeeze(overlay.ImageOverlay(Dat.SliceInd(ii),:,:,:));
          else
            overlay_data = squeeze(overlay.ImageOverlay(Dat.SliceInd(ii),:,:,Dat.currentVol));
          end
        end
      end
    end

    sz = size(imdata);
		
		% A fix for Windows 7. It appears that in Windows 7 one cannot create a
		% figure window that is narrower than the buttons in the top right
		% corner...
		if ispc && any(sz<128)
			sz = ceil((128/max(sz))*sz);
		end

    % Create invisible figure
    if drawColorbar
      sz(2)=sz(2)+65;
    end
    fh = figure('position',[300 200 sz(2) sz(1)],...
      'visible','off',...
      'inverthardcopy','off',...
      'renderer','opengl',...
      'paperpositionmode','auto',...
      'color','w');
    if ( roiTransp==0 | not(drawRois)) && not(drawOverlay)
      set(fh,'renderer','painters');
    end

    % Create axes
    if isempty(Clim_in)
      clim = [min(min(imdata)) max(max(imdata))];
    elseif size(Clim_in,1)>1
      clim = Clim_in(Dat.SliceInd(ii),:);
    else
      clim = Clim_in;
    end
    ax = axes('parent',fh,...
      'units','normal',...
      'position',[0 0 1 1],...
      'clim',clim,...
      'xlim',[0.5 sz(2)+0.5],...
      'ylim',[0.5 sz(1)+0.5],...
      'ydir','reverse',...
      'visible','off');
    if strcmpi(get(fh,'renderer'),'opengl')
      roi_ax = axes('parent',fh,...
        'units','normal',...
        'position',[0 0 1 1],...
        'visible','off',...
        'ydir','reverse',...
        'xlim',[0.5 sz(2)+0.5],...
        'ylim',[0.5 sz(1)+0.5]);
      overlay_ax = axes('parent',fh,...
        'units','normal',...
        'units','normal',...
        'position',[0 0 1 1],...
        'visible','off',...
        'ydir','reverse',...
        'xlim',[0.5 sz(2)+0.5],...
        'ylim',[0.5 sz(1)+0.5]);
    end

    % Create image
    im = image('parent',ax,...
      'cdatamapping','scaled',...
      'cdata',imdata);
    set(fh,'colormap',colmap)

    % Draw Colorbar
    if drawColorbar
      colorbar('peer',ax,'East',...
        'fontsize',7);
    end

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

    % Draw ROIs
    if drawRois
      for kk=1:length(roi_ind)
        if Dat.isDataMixed
          roidata = ROI(roi_ind(kk)).voxels{Dat.SliceInd(ii)};
        else
          if sliceDir==3
            roidata = ROI(roi_ind(kk)).voxels{1}(:,:,Dat.SliceInd(ii),Dat.currentVol);
          elseif sliceDir==2
            roidata = squeeze(ROI(roi_ind(kk)).voxels{1}(:,Dat.SliceInd(ii),:,Dat.currentVol));
          else
            roidata = squeeze(ROI(roi_ind(kk)).voxels{1}(Dat.SliceInd(ii),:,:,Dat.currentVol));
          end
          if isempty(find(roidata))
            continue;
          end
        end

        if strcmpi(get(fh,'renderer'),'opengl')
          roidata=uint8(roidata);
          cdata = zeros([size(roidata) 3],'uint8');
          if ROI(roi_ind(kk)).color(1)~=0
            cdata(:,:,1) = roidata*ROI(roi_ind(kk)).color(1);
          end
          if ROI(roi_ind(kk)).color(2)~=0
            cdata(:,:,2) = roidata*ROI(roi_ind(kk)).color(2);
          end
          if ROI(roi_ind(kk)).color(3)~=0
            cdata(:,:,3) = roidata*ROI(roi_ind(kk)).color(3);
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
              'color',ROI(roi_ind(kk)).color./255,...
              'tag','roiedge',...
              'linewidth',0.5,...
              'linestyle','-',...
              'hittest','off');
          end
        end
      end
    end

    % Draw slice number
    if drawSliceNbr
      if drawFileName
        slicenbr = sprintf('%03d:',Dat.SliceInd(ii));
      else
        slicenbr = sprintf('%03d',Dat.SliceInd(ii));
      end
      tx_nbr=text('parent',ax,...
        'units','normal',...
        'position',[0.015 0.965],...
        'interpreter','none',...
        'string',slicenbr,...
        'backgroundcolor','w',...
        'clipping','off',...
        'fontsize',7);
    end

    % Draw file name text
    if drawFileName
      if drawSliceNbr
        tmp_ext=get(tx_nbr,'extent');
      else
        tmp_ext=[0.015 0.965 0 0];
      end
      if Dat.isDataMixed
        fname = DATA{Dat.SliceInd(ii)}.HDR.fname;
        fpath = DATA{Dat.SliceInd(ii)}.HDR.fpath;
      else
        fname=DATA{1}.HDR.fname;
        fpath = DATA{1}.HDR.fpath;
      end
      if strcmp(fname,'fid')
        [fp,fn,fe]=fileparts(fpath(1:end-1));
        fname = [fn,fe];
      end
      tx_fname=text('parent',ax,...
        'units','normal',...
        'position',[tmp_ext(1)+tmp_ext(3) 0.965],...
        'interpreter','none',...
        'string',fname,...
        'backgroundcolor','w',...
        'clipping','off',...
        'fontsize',7);
    end

    % $$$     % Draw slice number
    % $$$     if drawSliceNbr
    % $$$       tx=text('parent',ax,...
    % $$$               'units','normal',...
    % $$$               'position',[0.015 0.965],...
    % $$$               'string',sprintf('%03d',Dat.SliceInd(ii)),...
    % $$$               'backgroundcolor','w',...
    % $$$               'fontsize',7);
    % $$$     end

    % Print figure
    print(fh,['-d',print_type,quality],['-r',resolution],[outdir,fnames{ii}]);

    delete(fh);
  end

  % Delete waitbar
  delete(wbar_h)
  
end

end % function l_CheckOptions(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_Close(h,evd)

% Store figure handle before clearing variables
fig_h = H.MAINFIG;

% Clear public variables
clear H Dat DATA ROI

% Close figure
delete(fig_h);

end % function l_Close(h,

end
