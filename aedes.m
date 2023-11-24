function aedes(DATA,ROI,inputOverlay)
% AEDES - Graphical user interface for analyzing MRI images
%   
%
% Synopsis: 
%       aedes(InArg)
%       
%       or
%
%       aedes
%
% Description:
%       AEDES is a graphic user interface (GUI) based tool for region
%       of interest (ROI) analysis of (mainly) MRI images.
%
%       The input argument InArg is optional. If the function is called
%       without input arguments, an empty GUI window is initialized. InArg
%       can be a DATA structure, containing at least the FTDATA and HDR
%       fields, and a string or a cell array of strings containing the full
%       path(s) to supported data file(s). Input argument can also be a 2D,
%       3D, or 4D-matrix containing image slices and possibly volumes in the
%       4th dimension.
%
% Examples:
%      %% Example 1
%      aedes  % Initialize an empty GUI window
%
%      %% Example 2
%      DataMtx = rand(100); % Generate a 100x100 random matrix
%      aedes(DataMtx)    % Open the data in Aedes
%
%      %% Example 3
%      DATA=aedes_data_read;      % Read various image data to a DATA-structure
%      aedes(DATA)       % Open data in Aedes
%
%      %% Example 4
%      % Read data from file with path given as a string
%      aedes('C:\MyDataDirectory\MyDataFiles\MyData.fid')
%
%      %% Example 5 
%      % Read data from multiple files given as a cell
%      % array. Each file has to contain only one slice.
%      files={'C:\MyData\MyDataFile1.nii','C:\MyData\MyDataFile1.nii'};
%      aedes(files)
%
% See also:
%      AEDES_DATA_READ, AEDES_READFID, AEDES_RESVIEWER, AEDES_JUIGETFILES

% Aedes - A graphical tool for analyzing medical images
%
% Copyright (C) 2006 Juha-Pekka Niskanen <Juha-Pekka.Niskanen@uku.fi>
% 
% Department of Physics, Department of Neurobiology
% University of Eastern Finland, Kuopio, FINLAND
%
% This program may be used under the terms of the GNU General Public
% License version 2.0 as published by the Free Software Foundation
% and appearing in the file LICENSE.TXT included in the packaging of
% this program.
%
% This program is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
% WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.


% Define "public" variables. Since Aedes utilizes nested functions,
% these variables will be visible to all the subfunctions.
H = [];   % Handles structure. Contains all handles to figure,
          % uicontrols, axes, etc...
Dat = []; % Internal data structure for exchanging information between
          % functions
switch nargin
  case 0;
    DATA = []; % The data structure containing image data and header
    ROI = []; % The ROI structure
    inputOverlay = [];
  case 1
    ROI = [];
    inputOverlay = [];
  case 2
    inputOverlay = [];
end


		  
% % Auto-check for updates ------------------------------
% try
%   if getpref('Aedes','AutoCheckUpdates')
% 	Dat=aedes_update('semiprompt');
% 	if Dat
% 	  return
% 	else
% 	  Dat = [];
% 	end
%   end
% catch
%   % pass
% end

% Detect Matlab version
[Dat.MatlabVersion,Dat.isImageProc] = aedes_getmatlabversion;


% Show license and warranty notification ----------------
%
% NOTE: You can disable this license notification using the following
% command: setpref('Aedes','ShowLicenseAtStartUp',false)
if ~ispref('Aedes','ShowLicenseAtStartUp') || ...
	getpref('Aedes','ShowLicenseAtStartUp')
  l_PrintLicense([],[]);
end


%% Parse input arguments -------------------------
if nargin>0
  if ischar(DATA) %% String input
    tmp_filename = DATA;
    DATA = [];
    
    % Draw blank GUI
    H=l_draw_gui();

    % Read data and initialize
    l_OpenFile([],[],'single',tmp_filename)
    
    % Initialize GUI
    %l_Initialize([],[])
    
  elseif iscell(DATA) %% Cell input
    if ischar(DATA{1})
      tmp_filename = DATA;
      DATA = [];
      
      % Draw blank GUI
      H=l_draw_gui();
      
      % Read data and initialize
      l_OpenFile([],[],'multi',tmp_filename)
      
      % Initialize GUI
      %l_Initialize([],[])
      
    elseif ~isstruct(DATA{1})
      error('Invalid type for input argument!')
    else
      
      % Draw blank GUI
      H=l_draw_gui();
      
      % Initialize GUI
      l_Initialize([],[])
	  
    end
  else
    if isstruct(DATA) %% Structure input
      DATA = {DATA};
	  
      % Check if kspace should be viewed
      if isfield(DATA{1},'KSPACE') && ...
		  (isempty(DATA{1}.FTDATA) & ~isempty(DATA{1}.KSPACE))
    
    resp=questdlg(['The inputted DATA structure contains only k-space',...
      ' information. Do you want to view real, imaginary, or absolute',...
      ' part of the complex k-space?'],...
      'Complex data inputted',...
      'Real','Imaginary','Absolute','Absolute');
        if strcmpi(resp,'Real')
          DATA{1}.FTDATA = real(DATA{1}.KSPACE);
        elseif strcmpi(resp,'Imaginary')
          DATA{1}.FTDATA = imag(DATA{1}.KSPACE);
        elseif strcmpi(resp,'Absolute')
          DATA{1}.FTDATA = abs(DATA{1}.KSPACE);
        else
          % Canceled
          return
        end
		
        %DATA{1}.FTDATA = abs(DATA{1}.KSPACE);
      end
    elseif isnumeric(DATA) || islogical(DATA)
      if isreal(DATA)
        DATA_tmp=DATA;
      else
        resp=questdlg(['The inputted data is complex. Do you want ' ...
          'to view real, imaginary, or absolute part of the data?'],...
          'Complex data inputted',...
          'Real','Imaginary','Absolute','Absolute');
        if strcmpi(resp,'Real')
          DATA_tmp = real(DATA);
        elseif strcmpi(resp,'Imaginary')
          DATA_tmp = imag(DATA);
        else
          DATA_tmp = abs(DATA);
        end
      end
      DATA=[];
      DATA.FTDATA = DATA_tmp;
      DATA.HDR.fpath = '';
      DATA.HDR.fname = '';
      DATA={DATA};
      clear DATA_tmp
    else
      error('Invalid type for input argument')
    end
    
    % Draw blank GUI
    H=l_draw_gui();
    
    % Initialize GUI
    l_Initialize([],[])
    
  end
  
  % Load ROI if given as input argument
  if nargin>1
    if ~isempty(ROI)
      l_RoiLoad([],[],'');
    end
    
    if nargin==3 && ~isempty(inputOverlay)
      % Load overlay
      l_LoadImageOverlay([],[],'inputarg')
    end
  end
  
else
  % Open an empty GUI window, if called without input arguments
  H=l_draw_gui;
end

% Build plugins menu
l_BuildPluginsMenu();

% Wait for quit in standalone
if isdeployed
  waitfor(H.FIG)
end

% Debug function - paste debug information into workspace...
function l_debug(h,evd)
  
assignin('base','DATA',DATA);
assignin('base','H',H);
assignin('base','ROI',ROI);
assignin('base','Dat',Dat);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK FOR UPDATES
%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_CheckUpdates(h,evd)
  done = aedes_update('prompt');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%
function H=l_draw_gui()
  
%% Toggle debug on/off
debug=true;

%% Determine paths
tmp=which('aedes');
[fp,fn,fe]=fileparts(tmp);
Dat.AedesFolder = [fp,filesep];
Dat.PluginsFolder = [Dat.AedesFolder,'plugins',filesep];
Dat.IconsFolder = [Dat.AedesFolder,'icons',filesep];
Dat.TablibFolder = [Dat.AedesFolder,'tablib',filesep];

%% Figure Size Definitions ----------------
taskbar_sz = 70; % The windows taskbar height in pixels
fig_w = 1000; % Screen width in pixels
fig_h = 661;%768-taskbar_sz; % Screen height in pixels
dx = 1/fig_w;          % Horizontal space
dy = 1/fig_h;         % Vertical space
dxx = 5*dx;                % Horizontal big space
dyy = 5*dy;                % Vertical big space
btn_h = 25/fig_h;     % Push button height
btn_w = 75/fig_w;      % Push button width
btns_w = 25/fig_w;
popup_h = 25/fig_h;   % Popup menu height
popup_w = 100/fig_w;   % Popup menu width
edit_w = 45/fig_w;     % Edit box width
edit_h = 20/fig_h;    % Edit box height


%% Load default font and colors
if isunix
  DefaultColor = [239 235 222]/255;
else
  DefaultColor=get(0,'DefaultUicontrolBackgroundcolor');
end
GD=aedes_gui_defaults;
GD.col.frame = DefaultColor;
Dat.HG2graphics = GD.HG2graphics;

% Calculate default position for Aedes window
try
	% If multiple monitors are connected, draw Aedes on the monitor where the
	% mouse cursor is located
	scrsz = get(0,'MonitorPositions');
	ind = aedes_getcurrentmonitor;
	fig_pos = [scrsz(ind,3)/2-fig_w/2 scrsz(ind,4)/2-fig_h/2-20 fig_w fig_h];
	fig_pos(1) = fig_pos(1)+scrsz(ind,1);
catch
	% Get screen size
	scrsz=get(0,'ScreenSize');
	
	% Calculate gui position on the center of the screen
	fig_pos = [scrsz(3)/2-fig_w/2 scrsz(4)/2-fig_h/2-20 fig_w fig_h];
end

%% Check if other Aedes windows exist
figs=findall(0,'tag','aedes_main_fig');
if ~isempty(figs)
  tmp_pos=get(figs(1),'position');
  fig_pos(1)=tmp_pos(1)+15;
  fig_pos(2)=tmp_pos(2)-15;
end

%% Draw main figure
H.FIG=figure('Position',fig_pos, ...
             'Units','Pixel', ...
             'Name','Aedes 1.0', ...
             'Numbertitle','off', ...
             'Tag','aedes_main_fig', ...
             'Color',GD.col.mainfig, ...
             'Toolbar','none', ...
             'Menubar','none', ...
             'DockControls','off',...
             'renderer','OpenGL',...
             'KeyPressFcn',@l_KeyPressFcn,...
             'CloseRequestFcn',@l_quit,...
             'Handlevisibility','off');


if ~Dat.HG2graphics
	set(H.FIG,'DoubleBuffer','on')
end

% File Uimenu ---------------------------
file_h = uimenu('Label','File','Accelerator','F', ...
                'Parent',H.FIG);
H.FILEMENU_NEW = uimenu(file_h,'Label','New window',...
                        'Callback','aedes');
H.FILEMENU_OPEN_SINGLE = uimenu(file_h,'Label','&Open File',...
  'Accelerator','O', ...
  'Callback',{@l_OpenFile,'single'},'Tag', ...
  'open_single_file',...
  'separator','on');						 
H.FILEMENU_OPEN_MULTI = uimenu(file_h,'Label','Open &Multiple Files',...
                                'Callback',{@l_OpenFile,'multi'},'Tag', ...
                                'open_multiple_file');
H.FILEMENU_OPEN_RECENT = uimenu(file_h,'Label','Open Recent',...
  'Tag','open_recent',...
  'separator','on');
l_BuildRecentFilesMenu([],[],H.FILEMENU_OPEN_RECENT);

H.FILEMENU_SAVE_IMAGE = uimenu(file_h,'Label','Save Image Data as', ...
                               'Callback',@l_SaveImageData, ...
                               'Separator','on','userdata','','enable','off');
H.FILEMENU_SAVERES = uimenu(file_h,'Label','Save results as','Accelerator','S', ...
                            'Callback',{@l_SaveResults,[]}, ...
                            'Separator','off','userdata','','enable','off');
H.FILEMENU_EXPORT = uimenu(file_h,'Label','Export Image Data...', ...
                           'Callback',@l_ExportImages, ...
                           'Separator','on','userdata','','enable','off');
% $$$ H.FILEMENU_PRINT = uimenu(file_h,'Label','Print...', ...
% $$$                           'Callback',@l_Print, ...
% $$$                           'Separator','on','userdata','','enable','off');
closefile_h =  uimenu(file_h,'Label','Close File', ...
                      'Callback',@l_CloseFile, ...
                      'Separator','on','userdata','','enable','off');
H.FILEMENU_QUIT = uimenu(file_h,'Label','Exit Aedes','Accelerator','Q', ...
                         'Callback',@l_quit,'Separator','on');


% Edit Uimenu ---------------------------------
edit_h = uimenu('Label','Edit', ...
                'Parent',H.FIG,...
                'enable','off');
H.UIEDIT_IMSTACK = uimenu(edit_h,'Label','&Image Stack',...
                          'Callback',@l_EditImageStack,...
                          'separator','off',...
                          'enable','off');

% Unfold data
H.UNFOLD_DATA = uimenu(edit_h,'Label','Unfold data',...
                       'Callback',@l_UnfoldData,...
                       'separator','on',...
                       'enable','off');

% Rotate/flip view
rotate_h=uimenu(edit_h,'Label','Rotate/Flip Images','separator','on',...
                'enable','off');
H.UIMENU_ROTATEFLIP = rotate_h;

rot_h=uimenu(rotate_h,'Label','Rotate Current Slice');
tmp_str={'90','180','270'};
for ii=1:3
  H.ROTATE(ii) = uimenu(rot_h,'Label',[tmp_str{ii},char(176)],...
                        'callback',...
                        {@l_RotateFlip,tmp_str{ii}});
end

flip_h=uimenu(rotate_h,'Label','Flip Current Slice');
tmp_str={'FlipLR','FlipUD'};
for ii=1:2
  H.FLIP(ii) = uimenu(flip_h,'Label',tmp_str{ii},...
                      'callback',...
                      {@l_RotateFlip,tmp_str{ii}});
end
rotcustom_h = uimenu(rotate_h,'Label','Rotate/Flip Custom...',...
                     'callback',...
                     {@l_RotateFlip,'custom'},...
                     'separator','on');
H.UIROTRESET = uimenu(rotate_h,'Label','Reset Rotation/Flip',...
                      'callback',...
                      {@l_RotateFlip,'reset'},...
                      'separator','on');

% Rotate/Flip volumes ----------------------------------------
rotate3d_h=uimenu(edit_h,'Label','Rotate/Flip Volumes','separator','off',...
                'enable','off');
H.UIMENU_ROTATEFLIP3D = rotate3d_h;
rot3d_h=uimenu(rotate3d_h,'Label','Rotate');
flip3d_h=uimenu(rotate3d_h,'Label','Flip');

tmpx_h = uimenu(rot3d_h,'Label','X-dir');
tmpy_h = uimenu(rot3d_h,'Label','Y-dir');
tmpz_h = uimenu(rot3d_h,'Label','Z-dir');

tmp_str={'90','180','270'};
tmp_h = [tmpx_h,tmpy_h,tmpz_h];
tmp_str2 = {'X','Y','Z'};
for kk=1:3
  for ii=1:3
	H.ROTATE3D(kk,ii) = uimenu(tmp_h(kk),'Label',[tmp_str{ii},char(176)],...
	  'callback',...
	  {@l_RotateFlip,['3d',tmp_str2{kk},'_',tmp_str{ii}]});
  end
end
H.FLIP3D_X = uimenu(flip3d_h,'Label','X-dir',...
  'callback',{@l_RotateFlip,'flipX'});
H.FLIP3D_Y = uimenu(flip3d_h,'Label','Y-dir',...
  'callback',{@l_RotateFlip,'flipY'});
H.FLIP3D_Z = uimenu(flip3d_h,'Label','Z-dir',...
  'callback',{@l_RotateFlip,'flipZ'});
H.FLIP3D_V = uimenu(flip3d_h,'Label','V-dir',...
  'callback',{@l_RotateFlip,'flipV'},...
  'enable','off');


% View Uimenu ------------------------------
view_h = uimenu('Label','View', ...
                'Parent',H.FIG,...
                'enable','off');



% View XYZ (3D) direction (default)
H.UIVIEW_3D = uimenu(view_h,'Label','View 3D',...
                     'callback',{@l_ChangeView,0},...
                     'checked','on');
% View only X direction
H.UIVIEW_X = uimenu(view_h,'Label','View X direction',...
                    'callback',{@l_ChangeView,1},...
                    'separator','on');
% View only X direction
H.UIVIEW_Y = uimenu(view_h,'Label','View Y direction',...
                    'callback',{@l_ChangeView,2});
% View only X direction
H.UIVIEW_Z = uimenu(view_h,'Label','View Z direction',...
                    'callback',{@l_ChangeView,3});

% View axes ticks and grid
H.UIVIEW_GRID = uimenu(view_h,'Label','Show Grid',...
  'callback',{@l_ViewAxesUnits,'pixel'},...
  'separator','on',...
  'checked','off');

% View Header in aedes_headerbrowser
H.viewheader_h = uimenu(view_h,'Label','File Header',...
  'callback',@l_LaunchHeaderBrowser,...
  'enable','on','Accelerator','H',...
  'separator','on');

% Show Info text background
H.UIVIEW_INFOTXTBACKGROUND = uimenu(view_h,'Label','Info/ROI Text Background',...
                                    'callback',@l_InfoTextBackground,...
                                    'separator','on',...
                                    'checked','off');

% View Voxel Time Series
H.UIVIEW_TIMESERIES = uimenu(view_h,'Label','Voxel time-series',...
                             'callback',{@l_ShowTimeSeries,'toggle'},...
                             'separator','on',...
                             'checked','off');

% % Volume menu
% H.UIVIEW_VOL = uimenu(view_h,'Label','Select volume',...
%                       'separator','on');
% 

% $$$ H.UIVIEW_FOV_UNITS = uimenu(view_h,'Label','View FOV units',...
% $$$                             'callback',{@l_ViewAxesUnits,'FOV'},...
% $$$                             'checked','off');


% $$$ flipcustom_h=uimenu(rotate_h,'Label','Flip Custom',...
% $$$                     'callback',...
% $$$                     {@l_RotateZoom,'flip','custom'},...
% $$$                     'separator','off');



% % Zoom view
% tmp_str={'normalized',...
%          '0.5x',...
%          '1x',...
%          '2x',...
%          '3x',...
%          '4x',...
%          '5x',...
%          '6x',...
%          '7x',...
%          '8x',...
%          '9x',...
%          '10x',...
% 		 '11x',...
%          '12x',...
%          '13x',...
%          '14x',...
%          '15x',...
%          '16x',...
%          '17x',...
%          '18x',...
%          '19x',...
%          '20x'};
% tmp_fact = [0 0.5 1 2 3 4 5 6 7 8 9 10 ...
% 11 12 13 14 15 16 17 18 19 20];
% zoom_h=uimenu(view_h,'Label','Zoom',...
%   'separator','on');
% for ii=1:length(tmp_str)
%   H.UIZOOM(ii)=uimenu(zoom_h,'Label',tmp_str{ii},...
%                       'callback',{@l_Zoom,tmp_fact(ii)},...
%                       'tag',['uizoom_' num2str(tmp_fact(ii))],...
%                       'userdata',tmp_fact(ii));
%   if ii==2
%     set(H.UIZOOM(ii),'separator','on')
%   end
% end
% set(H.UIZOOM(3),'checked','on')


% TOOLS UIMENU -------------------------------------
tools_h = uimenu('Label','Tools', ...
                 'Parent',H.FIG,...
                 'enable','on');
vnmredit_h = uimenu(tools_h,'Label','Edit VNMR defaults', ...
  'callback','aedes_readfidprefs',...
  'enable','on');
resview_h = uimenu(tools_h,'Label','Results Viewer', ...  
  'callback','aedes_resviewer',...
  'enable','on',...
  'separator','on');
% update_h = uimenu(tools_h,'Label','Check for updates', ...
%   'callback',@l_CheckUpdates,...
%   'enable','on',...
%   'separator','on');
  

				 
% Overlay menu ---------------------------------------
H.UIOVERLAY = uimenu('Label','Overlay', ...
                     'Parent',H.FIG,...
                     'enable','off');

load_overlay_h = uimenu(H.UIOVERLAY,'Label','Load Image Overlay', ...
                        'separator','off',...
                        'enable','on');
uimenu(load_overlay_h,'Label','From File', ...
                        'callback',{@l_LoadImageOverlay,'file'},...
                        'separator','off',...
                        'enable','on');
uimenu(load_overlay_h,'Label','From Variable', ...
                        'callback',{@l_LoadImageOverlay,'var'},...
                        'separator','on',...
                        'enable','on');
H.UIOVERLAY_SAVE = uimenu(H.UIOVERLAY,'Label','Save Overlay to file', ...
                              'callback',{@l_OverlayControls,'save'},...
                              'separator','on',...
                              'enable','off');

H.UIOVERLAY_CONTROLS = uimenu(H.UIOVERLAY,'Label','Show Overlay Controls', ...
                              'callback',{@l_OverlayControls,'show'},...
                              'separator','on',...
                              'enable','off');
H.UIOVERLAY_DELETE = uimenu(H.UIOVERLAY,'Label','Delete Image Overlay', ...
                            'callback',{@l_OverlayControls,'delete'},...
                            'separator','on',...
                            'enable','off');



% ROI Tools ------------------------------------------
roi_tools_h = uimenu('Label','ROI', ...
                     'Parent',H.FIG,...
                     'enable','off');

% $$$ H.UIROITOOLS_AUTOCLOSE = uimenu(roi_tools_h,'Label','Autoclose',...
% $$$                                 'callback',...
% $$$                                 ['if strcmp(get(gcbo,''checked''),''on''),' ...
% $$$                     'set(gcbo,''checked'',''off''),else,' ...
% $$$                     'set(gcbo,''checked'',''on''),end'],...
% $$$                                 'checked','off');

H.UIROITOOLS_LOAD = uimenu(roi_tools_h,'Label','Load ROI(s)',...
                           'callback',{@l_RoiLoad,'interactive'},...
                           'separator','off');
H.UIROITOOLS_SAVE = uimenu(roi_tools_h,'Label','Save ROI(s)',...
                           'callback',@l_RoiSave,...
                           'separator','off',...
                           'enable','off');
H.UIROITOOLS_SAVETEMPLATE = uimenu(roi_tools_h,'Label','Save ROI template',...
                                   'callback',{@l_RoiSave,'template'},...
                                   'separator','off',...
                                   'enable','off');
H.UIROISTATS = uimenu(roi_tools_h,'Label','View/Export ROI Statistics',...
                      'callback',@l_RoiViewStats,...
                      'separator','on',...
                      'enable','off');
H.UIROITOOLS_COPY = uimenu(roi_tools_h,'Label','Copy ROI',...
                           'separator','on',...
                           'enable','off');
H.UIROITOOLS_COPYSLICES = uimenu(H.UIROITOOLS_COPY,'Label','Copy ROI to slices/volumes',...
                                 'callback',@l_RoiCopy,...
                                 'separator','off',...
                                 'enable','on');
H.UIROITOOLS_COPYADD = uimenu(H.UIROITOOLS_COPY,'Label','Copy current ROI to new',...
                              'callback',{@l_RoiAdd,'copy'},...
                              'separator','off',...
                              'enable','on');
H.UIROITOOLS_FLIP = uimenu(roi_tools_h,'Label','Flip current ROI',...
                           'separator','off',...
                           'enable','off');
tmp_str={'FlipLR','FlipUD'};
for ii=1:2
  H.ROIFLIP(ii) = uimenu(H.UIROITOOLS_FLIP,'Label',tmp_str{ii},...
                         'callback',...
                         {@l_RoiFlip,tmp_str{ii}});
end
H.UIROITOOLS_COMP = uimenu(roi_tools_h,'Label','Boolean operations',...
                          'callback',@l_RoiBooleanOperations,...
                          'separator','off',...
                          'checked','off',...
                          'enable','off');
H.UIROISHOWEDGES = uimenu(roi_tools_h,'Label','Show ROI edges',...
                          'callback',@l_UiRoiEdgesCB,...
                          'separator','off',...
                          'checked','off',...
                          'enable','off');
H.UIROIRENAME = uimenu(roi_tools_h,'Label','Rename current ROI',...
                       'callback',@l_RoiRename,...
                       'separator','off',...
                       'checked','off',...
                       'enable','off');
H.UIROISETCOLOR = uimenu(roi_tools_h,'Label','Set current ROI color',...
                       'callback',@l_RoiSetColor,...
                       'separator','off',...
                       'checked','off',...
                       'enable','off');                     

H.UIROITOOLS_UNDO = uimenu(roi_tools_h,'Label','Undo draw (0)',...
                           'callback',@l_RoiUndo,...
                           'separator','on',...
                           'enable','off',...
                           'accelerator','z');

%% Plugins menu --------------------------------------------------------
%% DO NOT ADD/REMOVE PLUGINS FROM HERE.

% Construct plugins uimenu
H.UIPLUGINS = uimenu('Label','Plugins', ...
  'Parent',H.FIG,...
  'enable','off');


%% ----------------------------------------------------------------------

%% Help uimenu
help_h = uimenu('Label','Help', ...
                'Parent',H.FIG,...
                'enable','on');
help_about_h = uimenu(help_h,...
                      'Label','About Aedes', ...
                      'enable','on',...
                      'callback','eval(get(gcbo,''userdata''))');
set(help_about_h,'userdata',...
                 'aedes_helpabout');


if debug
  debug_h=uimenu('Label','Debug', ...
                 'Parent',H.FIG,...
                 'enable','on',...
                 'callback',@l_debug);
end

%% Uitoolbar controls --------------------------------
H.UITOOLBAR = uitoolbar('parent',H.FIG);

%% load CData for buttons
btn_cdata = load('aedes_cdata.mat');
H.btn_cdata = btn_cdata;
%tmp=NaN(18,20,3);
%tmp(:,3:end,:)=btn_cdata.cdata.open;

% Open file
H.UIPUSH_OPEN = uipushtool('parent',H.UITOOLBAR,...
                           'CData',btn_cdata.cdata.open2,...
                           'ClickedCallback',{@l_OpenFile,'single'},...
                           'TooltipString','Open a new data file');

% Open multiple files
H.UIPUSH_OPENMULTI = uipushtool('parent',H.UITOOLBAR,...
                                'CData',btn_cdata.cdata.openmulti,...
                                'ClickedCallback',{@l_OpenFile,'multi'},...
                                'TooltipString',...
                                'Open multiple single-slice data files');

% Save file
H.UIPUSH_SAVE = uipushtool('parent',H.UITOOLBAR,...
                           'CData',btn_cdata.cdata.save2_small,...
                           'ClickedCallback',{@l_SaveResults,[]},...
                           'TooltipString','Save',...
                           'enable','off');

% Print
% $$$ H.UIPUSH_PRINT = uipushtool('parent',H.UITOOLBAR,...
% $$$                             'CData',btn_cdata.cdata.print2,...
% $$$                             'ClickedCallback','',...
% $$$                             'TooltipString','Print',...
% $$$                             'enable','off');

% Zoom Normalized
H.UITOGGLE_ZOOMNORM = uitoggletool('parent',H.UITOOLBAR,...
                                   'CData',btn_cdata.cdata.fitheight,...
                                   'separator','on',...
                                   'ClickedCallback',{@l_Zoom,'normalize'},...
                                   'TooltipString','Fit to Window Height',...
                                   'enable','off');

% Zoom in
H.UIPUSH_ZOOMIN = uipushtool('parent',H.UITOOLBAR,...
                             'CData',btn_cdata.cdata.zoomin,...
                             'separator','off',...
                             'ClickedCallback',{@l_Zoom,'+'},...
                             'TooltipString','Zoom In',...
                             'enable','off');

% Zoom out
H.UIPUSH_ZOOMOUT = uipushtool('parent',H.UITOOLBAR,...
                              'CData',btn_cdata.cdata.zoomout,...
                              'ClickedCallback',{@l_Zoom,'-'},...
                              'TooltipString','Zoom Out',...
                              'enable','off');

% Show/hide crosshairs
H.UITOGGLE_CROSSHAIRS = uitoggletool('parent',H.UITOOLBAR,...
                               'CData',btn_cdata.cdata.crosshairs_small,...
                               'ClickedCallback',@l_ShowHideCrossbars,...
                               'TooltipString','Show/Hide Crosshairs',...
                               'enable','off',...
                               'state','off',...
                               'separator','on');
														
% Show Grid
H.UITOGGLE_GRID = uitoggletool('parent',H.UITOOLBAR,...
                               'CData',btn_cdata.cdata.grid,...
                               'ClickedCallback',{@l_ViewAxesUnits,'pixel'},...
                               'TooltipString','Show/Hide Grid',...
                               'enable','off',...
                               'state','off',...
                               'separator','on');

% Show/hide colorbar
H.UITOGGLE_COLORBAR = uitoggletool('parent',H.UITOOLBAR,...
                                   'CData',btn_cdata.cdata.colorbar,...
                                   'ClickedCallback',@l_ShowColorbar,...
                                   'TooltipString','Show/Hide Colorbar',...
                                   'enable','off',...
                                   'state','off',...
                                   'separator','off');

% Show Info Text
H.UITOGGLE_INFO = uitoggletool('parent',H.UITOOLBAR,...
                               'CData',btn_cdata.cdata.about,...
                               'ClickedCallback',@l_UpdateInfoText,...
                               'TooltipString','Show/Hide Info Text',...
                               'enable','off',...
                               'state','on',...
                               'separator','off');

% Show image stack
H.UIPUSH_IMSTACK = uipushtool('parent',H.UITOOLBAR,...
                              'CData',btn_cdata.cdata.editimagestack,...
                              'ClickedCallback',@l_EditImageStack,...
                              'TooltipString','Edit Image Stack',...
                              'enable','off',...
                              'separator','on');

% Show rotate/flip custom menu
H.UIPUSH_ROTATE = uipushtool('parent',H.UITOOLBAR,...
                             'CData',btn_cdata.cdata.rotateflip,...
                             'ClickedCallback',{@l_RotateFlip,'custom'},...
                             'TooltipString','Rotate/Flip Custom',...
                             'enable','off',...
                             'separator','off');


% View 3D
H.UITOGGLE_VIEW3D = uitoggletool('parent',H.UITOOLBAR,...
                               'CData',btn_cdata.cdata.view3d,...
                               'ClickedCallback',{@l_ChangeView,0},...
                               'TooltipString','View 3D',...
                               'enable','off',...
                               'state','off',...
                               'separator','on');

H.UITOGGLE_VIEWX = uitoggletool('parent',H.UITOOLBAR,...
                               'CData',btn_cdata.cdata.viewx,...
                               'ClickedCallback',{@l_ChangeView,1},...
                               'TooltipString','View X direction',...
                               'enable','off',...
                               'state','off',...
                               'separator','off');

H.UITOGGLE_VIEWY = uitoggletool('parent',H.UITOOLBAR,...
                               'CData',btn_cdata.cdata.viewy,...
                               'ClickedCallback',{@l_ChangeView,2},...
                               'TooltipString','View Y direction',...
                               'enable','off',...
                               'state','off',...
                               'separator','off');

H.UITOGGLE_VIEWZ = uitoggletool('parent',H.UITOOLBAR,...
                               'CData',btn_cdata.cdata.viewz,...
                               'ClickedCallback',{@l_ChangeView,3},...
                               'TooltipString','View Z direction',...
                               'enable','off',...
                               'state','off',...
                               'separator','off');

% Mouse uitoggletools
H.UITOGGLE_ARROW = uitoggletool('parent',H.UITOOLBAR,...
                               'CData',btn_cdata.cdata.arrow,...
                               'ClickedCallback',...
                                ['set(get(gcbo,''userdata''),''state'',''off''),' ...
                    'set(gcbo,''state'',''on'')'],...
                               'TooltipString','Mouse: Normal',...
                               'enable','off',...
                               'state','on',...
                               'separator','on');

H.UITOGGLE_PAN = uitoggletool('parent',H.UITOOLBAR,...
                               'CData',btn_cdata.cdata.pan,...
                               'ClickedCallback',...
                              ['set(get(gcbo,''userdata''),''state'',''off''),' ...
                    'set(gcbo,''state'',''on'')'],...
                              'TooltipString','Mouse: Pan image',...
                              'enable','off',...
                              'state','off',...
                              'separator','off');

% ROI buttons
H.UITOGGLE_DRAWROI = uitoggletool('parent',H.UITOOLBAR,...
                                  'CData',btn_cdata.cdata.freedraw,...
                                  'ClickedCallback',...
                                  ['set(get(gcbo,''userdata''),''state'',''off''),' ...
                    'set(gcbo,''state'',''on'')'],...
                                  'TooltipString','Mouse: Draw ROI',...
                                  'enable','off',...
                                  'state','off',...
                                  'separator','off');
H.UITOGGLE_ERASEROI = uitoggletool('parent',H.UITOOLBAR,...
                                  'CData',btn_cdata.cdata.eraser,...
                                   'ClickedCallback',...
                                   ['set(get(gcbo,''userdata''),''state'',''off''),' ...
                    'set(gcbo,''state'',''on'')'],...
                                   'TooltipString','Mouse: Erase ROI',...
                                   'enable','off',...
                                   'state','off',...
                                   'separator','off');

H.UITOGGLE_FILLROI = uitoggletool('parent',H.UITOOLBAR,...
                                  'CData',btn_cdata.cdata.bucketfill,...
                                  'ClickedCallback',...
                                  ['set(get(gcbo,''userdata''),''state'',''off''),' ...
                    'set(gcbo,''state'',''on'')'],...
                                   'TooltipString','Mouse: Fill ROI',...
                                   'enable','off',...
                                   'state','off',...
                                   'separator','off');

H.UITOGGLE_MOVEROI = uitoggletool('parent',H.UITOOLBAR,...
                                  'CData',btn_cdata.cdata.move,...
                                  'ClickedCallback',...
                                  ['set(get(gcbo,''userdata''),''state'',''off''),' ...
                    'set(gcbo,''state'',''on'')'],...
                                   'TooltipString','Mouse: Move ROI in slice',...
                                   'enable','off',...
                                   'state','off',...
                                   'separator','off');
H.UITOGGLE_MOVEROI3D = uitoggletool('parent',H.UITOOLBAR,...
                                    'CData',btn_cdata.cdata.move3d,...
                                    'ClickedCallback',...
                                    ['set(get(gcbo,''userdata''),''state'',''off''),' ...
                    'set(gcbo,''state'',''on'')'],...
                                    'TooltipString','Mouse: Move ROI in volume',...
                                    'enable','off',...
                                    'state','off',...
                                    'separator','off');

% Show ROI statistics
H.UIPUSH_STATS_ROI = uipushtool('parent',H.UITOOLBAR,...
                                'CData',btn_cdata.cdata.statistics,...
                                'ClickedCallback',@l_RoiViewStats,...
                                'TooltipString','View ROI statistics',...
                                'enable','off',...
                                'separator','on');

H.UIPUSH_NEWROI = uipushtool('parent',H.UITOOLBAR,...
                             'CData',btn_cdata.cdata.add_small,...
                             'ClickedCallback',@l_RoiAdd,...
                             'TooltipString','New ROI',...
                             'enable','off',...
                             'separator','on');

H.UIPUSH_COPYROI = uipushtool('parent',H.UITOOLBAR,...
                             'CData',btn_cdata.cdata.copy_small,...
                             'ClickedCallback',@l_RoiCopy,...
                              'TooltipString','Copy ROI to slices/volumes',...
                              'enable','off',...
                              'separator','off');

% $$$ H.UIPUSH_DELETEROI = uipushtool('parent',H.UITOOLBAR,...
% $$$                                 'CData',btn_cdata.cdata.delete_small,...
% $$$                                 'ClickedCallback',{@l_RoiFunctions,'delete'},...
% $$$                                 'TooltipString','Delete selected ROI(s)',...
% $$$                                 'enable','off',...
% $$$                                 'separator','off');
% $$$ 
% $$$ H.UIPUSH_DELETEALLROI = uipushtool('parent',H.UITOOLBAR,...
% $$$                                 'CData',btn_cdata.cdata.deleteall_small,...
% $$$                                 'ClickedCallback',{@l_RoiFunctions,'delete_all'},...
% $$$                                 'TooltipString','Delete ALL ROI(s)',...
% $$$                                 'enable','off',...
% $$$                                 'separator','off');

H.UIPUSH_UNDOROI = uipushtool('parent',H.UITOOLBAR,...
                              'CData',btn_cdata.cdata.undo_small,...
                              'ClickedCallback',@l_RoiUndo,...
                              'TooltipString','Undo last ROI action',...
                              'enable','off',...
                              'separator','off');






set([H.UITOGGLE_ARROW,H.UITOGGLE_PAN,...
     H.UITOGGLE_DRAWROI,H.UITOGGLE_ERASEROI,H.UITOGGLE_MOVEROI,...
     H.UITOGGLE_MOVEROI3D,H.UITOGGLE_FILLROI],...
    'userdata',[H.UITOGGLE_ARROW,H.UITOGGLE_PAN,...
                H.UITOGGLE_DRAWROI,H.UITOGGLE_ERASEROI,H.UITOGGLE_MOVEROI,...
                H.UITOGGLE_MOVEROI3D,H.UITOGGLE_FILLROI])


fig_pos = get(H.FIG,'position');
H.ORIG_FIG_POS = fig_pos;
%% Draw uicontrol frame behind the uicontrols
%  to make them draw faster in opengl mode
H.SIDEBAR_FRAME=uicontrol('parent',H.FIG,...
                          'units','pixel',...
                          'position',[0 0 232 fig_pos(4)+2],...
                          'backgroundcolor',GD.col.frame,...
                          'style','frame',...
                          'visible','on');


%% Draw sidebar -------------------------------
% H.SIDEBAR=uipanel('parent',H.FIG,...
%                   'units','pixel',...
%                   'position',[0 0 232 fig_pos(4)+1],...
%                   'backgroundcolor',GD.col.frame);

%% Draw image sliders
%tmp=get(H.SIDEBAR,'position');
tmp=get(H.SIDEBAR_FRAME,'position');
% Draw frame to hide the opening in the corner
H.IMSLIDER_FRAME = uicontrol('parent',H.FIG,...
                             'units','pixel',...
                             'position',[tmp(1)+tmp(3)+1 0 fig_pos(3)-tmp(1)-tmp(3) 17],...
                             'style','frame',...
                             'backgroundcolor',GD.col.mainfig,...
                             'visible','off');
% x-slider
H.IMSLIDER(1) = uicontrol('parent',H.FIG,...
                          'units','pixel',...
                          'position',[tmp(1)+tmp(3) 0 fig_pos(3)-tmp(1)-tmp(3)-17 17],...
                          'style','slider',...
                          'min',0,...
                          'max',1,...
                          'visible','off');%'callback',@l_AxesPositions);%@l_ImSliderCB);
                                           % y-slider
H.IMSLIDER(2) = uicontrol('parent',H.FIG,...
                          'units','pixel',...
                          'position',[fig_pos(3)-17 17 17 fig_pos(4)-17+1],...
                          'style','slider',...
                          'min',-1,...
                          'max',0,...
                          'visible','off');%'callback',@l_AxesPositions);%@l_ImSliderCB);

% If JavaFigures are enabled, set image sliders to work while moving
if ~Dat.HG2graphics
if Dat.MatlabVersion>7.03
  SliderListener1 = handle.listener(H.IMSLIDER(1),...
	'ActionEvent',...
	@l_AxesPositions);%@l_ImSliderCB);
  SliderListener2 = handle.listener(H.IMSLIDER(2),...
	'ActionEvent',...
	@l_AxesPositions);%@l_ImSliderCB);
  setappdata(H.IMSLIDER_FRAME,'ImSliderListener',...
	[SliderListener1,...
	SliderListener2])
else
  if  feature('javafigures')
	SliderListener1 = handle.listener(H.IMSLIDER(1),...
	  'ActionEvent',...
	  @l_AxesPositions);%@l_ImSliderCB);
	SliderListener2 = handle.listener(H.IMSLIDER(2),...
	  'ActionEvent',...
	  @l_AxesPositions);%@l_ImSliderCB);
	setappdata(H.IMSLIDER_FRAME,'ImSliderListener',...
	  [SliderListener1,...
	  SliderListener2])
  else
	set(H.IMSLIDER,'callback',@l_AxesPositions)
  end
end
else
	SliderListener1 = addlistener(H.IMSLIDER(1),...
		'ContinuousValueChange',@l_AxesPositions);
	SliderListener2 = addlistener(H.IMSLIDER(2),...
		'ContinuousValueChange',@l_AxesPositions);
	setappdata(H.IMSLIDER_FRAME,'ImSliderListener',[SliderListener1,...
		SliderListener2])
end

%% Slice viewing tools uipanel
% H.SLICE_TOOLS=uipanel('parent',H.FIG,...H.SIDEBAR,...
%                       'units','pixel',...
%                       'position',[round(0.01*232) round(0.79*fig_pos(4)) ...
%                     round(0.98*232)-2 round(0.18*fig_pos(4))],...
%                       'title','Slice viewing tools',...
%                       'fontsize',8,...
%                       'backgroundcolor',GD.col.frame,...%'highlightcolor',[0.5 0.5 0.5],...%'shadowcolor',[0 0 0],...
%                       'fontweight','bold',...
%                       'bordertype','etchedin');
fig_pos=get(H.FIG,'position');
H.SLICE_TOOLS=uicontrol('parent',H.FIG,...
                        'style','frame',...
                        'units','pixel',...
                        'position',[5 fig_pos(4)-105-2 ...
                    232-10 105],...%135
                        'backgroundcolor',GD.col.frame);


% Slice Tools text
tmp=get(H.SLICE_TOOLS,'position');
H.SLICE_TOOLS_TX = uicontrol('parent',H.FIG,...
                             'units','pixel',...
                             'position',[tmp(1)+5 tmp(2)+tmp(4)-20 150 15],...
                             'style','text',...
                             'horizontalalign','left',...
                             'fontweight','bold',...
                             'backgroundcolor',GD.col.frame,...
                             'string','Image Coordinates');

% Volume popup
tmp = get(H.SLICE_TOOLS,'position');
% $$$ H.VOL_POPUP = uicontrol('parent',H.FIG,...
% $$$                         'units','pixel',...%'normal',...
% $$$                         'position',[85 ...
% $$$                     80+tmp(2) ...
% $$$                     45 ...
% $$$                     34],...%[0.37 0.8 0.2 0.15],...
% $$$                         'style','popup',...
% $$$                         'backgroundcolor','w',...
% $$$                         'string',{'1'},...
% $$$                         'callback','');
% $$$ 
% $$$ % Volume text
% $$$ h1 = uicontrol('parent',H.FIG,...
% $$$                'units','pixel',...
% $$$                'position',[tmp(1)+5 tmp(2)+95 72 15],...
% $$$                'style','text',...
% $$$                'horizontalalign','left',...
% $$$                'backgroundcolor',GD.col.frame,...
% $$$                'string','Select volume');

% Slice slider
H.SL_SLIDER = uicontrol('parent',H.FIG,...
                        'units','pixel',...
                        'position',[tmp(1)+5 tmp(2)+60 215 20],...
                        'style','slider',...
                        'backgroundcolor',GD.col.slider);

% If JavaFigures are enabled, set slider to work while moving
if ~Dat.HG2graphics
if Dat.MatlabVersion>7.03
  SliderListener = handle.listener(H.SL_SLIDER,...
	'ActionEvent',...
	  {@l_ChangeSlice,'slider'});
	setappdata(H.SL_SLIDER,'SliderListener',SliderListener)
else
  if feature('javafigures')
	SliderListener = handle.listener(H.SL_SLIDER,...
	  'ActionEvent',...
	  {@l_ChangeSlice,'slider'});
	setappdata(H.SL_SLIDER,'SliderListener',SliderListener)
  else
	set(H.SL_SLIDER,'callback',{@l_ChangeSlice,'slider'})
  end
end
else
	SliderListener = addlistener(H.SL_SLIDER,...
		'ContinuousValueChange',@l_ChangeSlice);
	setappdata(H.SL_SLIDER,'SliderListener',SliderListener)
end


% X --------------------------
H.X_BTN = uicontrol('parent',H.FIG,...
                    'units','pixel',...
                    'position',[tmp(1)+5 tmp(2)+35 30 22],...
                    'style','togglebutton',...
                    'fontweight','bold',...
                    'string','X:',...
                    'value',1,...
                    'callback',{@l_ChangeSliderOrient,1});
tmp=get(H.X_BTN,'pos');

H.X_EDIT = uicontrol('parent',H.FIG,...
                     'units','pixel',...
                     'position',[5+tmp(3)+5 tmp(2) 35 tmp(4)],...
                     'style','edit',...
                     'backgroundcolor','w',...
                     'string','1',...
                     'callback',{@l_ChangeSlice,'editbox'});

% Y -----------------------------
tmp_gap = 8;
tmp=get(H.X_EDIT,'pos');
H.Y_BTN = uicontrol('parent',H.FIG,...
                    'units','pixel',...
                    'position',[tmp(1)+tmp(3)+tmp_gap tmp(2) 30 tmp(4)],...
                    'style','togglebutton',...
                    'string','Y:',...
                    'fontweight','normal',...
                    'value',0,...
                    'callback',{@l_ChangeSliderOrient,2});
tmp=get(H.Y_BTN,'pos');

H.Y_EDIT = uicontrol('parent',H.FIG,...
                     'units','pixel',...
                     'position',[tmp(1)+tmp(3) tmp(2) 35 tmp(4)],...
                     'style','edit',...
                     'backgroundcolor','w',...
                     'string','1',...
                     'callback',{@l_ChangeSlice,'editbox'});

% Z -----------------------------
tmp=get(H.Y_EDIT,'pos');
H.Z_BTN = uicontrol('parent',H.FIG,...
                    'units','pixel',...
                    'position',[tmp(1)+tmp(3)+tmp_gap tmp(2) 30 tmp(4)],...
                    'style','togglebutton',...
                    'fontweight','normal',...
                    'string','Z:',...
                    'value',0,...
                    'callback',{@l_ChangeSliderOrient,3});
tmp=get(H.Z_BTN,'pos');

H.Z_EDIT = uicontrol('parent',H.FIG,...
                     'units','pixel',...
                     'position',[tmp(1)+tmp(3) tmp(2) 35 tmp(4)],...
                     'style','edit',...
                     'backgroundcolor','w',...
                     'string','1',...
                     'callback',{@l_ChangeSlice,'editbox'});

% V ------------------------------
tmp=get(H.Z_BTN,'pos');
H.V_BTN = uicontrol('parent',H.FIG,...
                    'units','pixel',...
                    'position',[tmp(1) tmp(2)-tmp(4)-tmp_gap 30 tmp(4)],...
                    'style','togglebutton',...
                    'fontweight','normal',...
                    'string','V:',...
                    'value',0,...
                    'callback',{@l_ChangeSliderOrient,4},...
                    'enable','off');
tmp=get(H.V_BTN,'pos');

H.V_EDIT = uicontrol('parent',H.FIG,...
                     'units','pixel',...
                     'position',[tmp(1)+tmp(3) tmp(2) 35 tmp(4)],...
                     'style','edit',...
                     'backgroundcolor','w',...
                     'string','1',...
                     'callback',{@l_ChangeVolume,'editbox'},...
                     'enable','off');


% Crossbars
tmp=get(H.SLICE_TOOLS,'pos');
% H.SHOW_CROSSBARS = uicontrol('parent',H.FIG,...
%                              'units','pixel',...
%                              'position',[tmp(1)+5 tmp(2)+25 130 20],...
%                              'style','checkbox',...
%                              'fontweight','normal',...
%                              'string','Show Crossbars',...
%                              'backgroundcolor',GD.col.frame,...
%                              'value',0,...
%                              'callback',@l_ShowHideCrossbars);

H.VOXEL_VALUE_TX = uicontrol('parent',H.FIG,...
                             'units','pixel',...
                             'position',[tmp(1)+5 tmp(2)+20 130 13],...
                             'style','text',...
                             'fontweight','normal',...
                             'string','Current pixel value:',...
														 'horizontalalign','left',...
                             'backgroundcolor',GD.col.frame);													 

H.VOXEL_VALUE = uicontrol('parent',H.FIG,...
                             'units','pixel',...
                             'position',[tmp(1)+5 tmp(2)+5 130 13],...
                             'style','text',...
                             'fontweight','normal',...
                             'string','-',...
														 'horizontalalign','left',...
                             'backgroundcolor',GD.col.frame);	
													 
tmp=get(H.SLICE_TOOLS,'position');
fig_pos = get(H.FIG,'position');
%% Window/Level tools frame
H.WINDOWLEVEL_TOOLS=uicontrol('parent',H.FIG,...
                           'style','frame',...
                           'units','pixel',...
                           'position',[tmp(1) tmp(2)-5-220 ...
                    tmp(3) 220],...
                           'backgroundcolor',GD.col.frame);


% Window/Level text
tmp=get(H.WINDOWLEVEL_TOOLS,'position');
H.WINDOWLEVEL_TOOLS_TX = uicontrol('parent',H.FIG,...
                                'units','pixel',...
                                'position',[tmp(1)+5 tmp(2)+tmp(4)-18 150 13],...
                                'style','text',...
                                'horizontalalign','left',...
                                'fontweight','bold',...
                                'backgroundcolor',GD.col.frame,...
                                'string','Window/Level');


% Window slider ------------------------------
tmp=get(H.WINDOWLEVEL_TOOLS,'position');
H.WINDOW_SLIDER = uicontrol('parent',H.FIG,...
                              'units','pixel',...
                              'position',[tmp(1)+5 tmp(2)+165 173 18],...
                              'style','slider',...
                              'min',0,...
                              'max',100,...
                              'value',100,...
                              'sliderstep',[0.01 0.1],...
                              'backgroundcolor',GD.col.slider);
tmp=get(H.WINDOW_SLIDER,'position');
% Window Edit box
H.WINDOW_EDIT = uicontrol('parent',H.FIG,...
                            'units','pixel',...
                            'position',[tmp(1)+tmp(3)+3 tmp(2) 37 tmp(4)],...
                            'style','edit',...
                            'backgroundcolor','w',...
                            'string','',...
                            'callback',@l_SetWindowLevel);
% Window text
h1 = uicontrol('parent',H.FIG,...
               'units','pixel',...
               'position',[tmp(1) tmp(2)+tmp(4)+2 150 12],...
               'style','text',...
               'horizontalalign','left',...
               'backgroundcolor',GD.col.frame,...
               'string','Window (0% - 100%)');

% Level slider ---------------------------------
tmp=get(H.WINDOW_SLIDER,'position');
tmp_gap=20;
H.LEVEL_SLIDER = uicontrol('parent',H.FIG,...
                                'units','pixel',...
                                'position',[tmp(1) tmp(2)-tmp(4)-tmp_gap tmp(3) tmp(4)],...
                                'style','slider',...
                                'min',0,...
                                'max',100,...
                                'value',50,...
                                'sliderstep',[0.01 0.1],...
                                'backgroundcolor',GD.col.slider);
tmp=get(H.LEVEL_SLIDER,'position');
% Level Edit box
H.LEVEL_EDIT = uicontrol('parent',H.FIG,...
                              'units','pixel',...
                              'position',[tmp(1)+tmp(3)+3 tmp(2) 37 tmp(4)],...
                              'style','edit',...
                              'backgroundcolor','w',...
                              'string','',...
                              'callback',@l_SetWindowLevel);
% Level text
h1 = uicontrol('parent',H.FIG,...
               'units','pixel',...
               'position',[tmp(1) tmp(2)+tmp(4)+2 170 12],...
               'style','text',...
               'horizontalalign','left',...
               'backgroundcolor',GD.col.frame,...
               'string','Level (0% - 100%)');

if ~Dat.HG2graphics
if Dat.MatlabVersion>7.03
  WindowSliderListener = handle.listener(H.WINDOW_SLIDER,...
	  'ActionEvent',@l_SetWindowLevel);
	LevelSliderListener = handle.listener(H.LEVEL_SLIDER,...
	  'ActionEvent',@l_SetWindowLevel);
	setappdata(H.WINDOW_SLIDER,'HandleListeners',...
	  [WindowSliderListener,LevelSliderListener])
else
  if feature('javafigures')
	WindowSliderListener = handle.listener(H.WINDOW_SLIDER,...
	  'ActionEvent',@l_SetWindowLevel);
	LevelSliderListener = handle.listener(H.LEVEL_SLIDER,...
	  'ActionEvent',@l_SetWindowLevel);
	setappdata(H.WINDOW_SLIDER,'HandleListeners',...
	  [WindowSliderListener,LevelSliderListener])
  else
	set(H.WINDOW_SLIDER,'callback',@l_SetWindowLevel)
	set(H.LEVEL_SLIDER,'callback',@l_SetWindowLevel)
  end
end
else
	WindowSliderListener = addlistener(H.WINDOW_SLIDER,...
	  'ContinuousValueChange',@l_SetWindowLevel);
	LevelSliderListener = addlistener(H.LEVEL_SLIDER,...
	  'ContinuousValueChange',@l_SetWindowLevel);
	setappdata(H.WINDOW_SLIDER,'HandleListeners',...
	  [WindowSliderListener,LevelSliderListener])
end
tmp=get(H.LEVEL_SLIDER,'position');
% Clim Min text
H.CLIMRANGETEXT = uicontrol('parent',H.FIG,...
                          'units','pixel',...
                          'position',...
                          [10 tmp(2)-25 50 15],...
                          'string','Range:',...
                          'horizontalalign','left',...
                          'style','text',...
                          'backgroundcolor',GD.col.frame);
tmp=get(H.CLIMRANGETEXT,'position');
H.CLIMMINEDIT = uicontrol('parent',H.FIG,...
                          'units','pixel',...
                          'position',...
                          [tmp(1)+tmp(3)+5 tmp(2) 60 18],...
                          'string','',...
                          'style','edit',...
                          'backgroundcolor','w',...
                          'callback',@l_ChangeClim);
tmp=get(H.CLIMMINEDIT,'position');
% % Clim Max text
H.CLIMMAXTEXT = uicontrol('parent',H.FIG,...
                          'units','pixel',...
                          'position',...
                          [tmp(1)+tmp(3)+5 tmp(2) 25 15],...
                          'string','-',...
                          'horizontalalign','center',...
                          'style','text',...
                          'backgroundcolor',GD.col.frame);
% tmp=get(H.CLIMMAXTEXT,'position');
% Max (Clim(2))
H.CLIMMAXEDIT = uicontrol('parent',H.FIG,...
                          'units','pixel',...
                          'position',...
                          [tmp(1)+tmp(3)+35 tmp(2) 60 18],...
                          'string','',...
                          'style','edit',...
                          'backgroundcolor','w',...
                          'callback',@l_ChangeClim);
tmp=get(H.CLIMMINEDIT,'position');
% Clim range behavior text
H.CLIMRANGEBEHTEXT = uicontrol('parent',H.FIG,...
	'units','pixel',...
	'position',...
	[10 tmp(2)-25 200 15],...
	'string','Range in image stack:',...
	'horizontalalign','left',...
	'style','text',...
	'backgroundcolor',GD.col.frame);
tmp=get(H.CLIMMINEDIT,'position');
% Clim range popup
H.CLIMRANGE_POPUP = uicontrol('parent',H.FIG,...
                              'units','pixel',...
                              'position',...
                              [10 tmp(2)-45 210 20],...
                              'string',{'Range fixed (all images)',...
                    'Range fixed (image)','Range image min-max','Range auto-balanced W/L'},...
                              'horizontalalign','left',...
                              'style','popup',...
                              'enable','on',...
                              'value',1,...
                              'callback',@l_ChangeClimRange,...
                              'backgroundcolor','w');

% Window/level autobalance
tmp=get(H.CLIMRANGE_POPUP,'position');
H.WINDOWLEVEL_AUTO = uicontrol('parent',H.FIG,...
                            'units','pixel',...
                            'position',...
                            [tmp(1) tmp(2)-tmp(4)-30 35 35],...
                            'string','',...
                            'CData',btn_cdata.cdata.contrast,...
                            'tooltip','Window/level auto-balance',...
                            'style','pushbutton',...
                            'callback',@l_SetAutoWindowLevel);
tmp=get(H.WINDOWLEVEL_AUTO,'position');
H.CONTRAST_INVERT = uicontrol('parent',H.FIG,...
                              'units','pixel',...
                              'position',...
                              [tmp(1)+tmp(3)+5 tmp(2) 35 35],...
                              'string','',...
                              'CData',btn_cdata.cdata.invertcontrast,...
                              'style','togglebutton',...
                              'tooltip','Invert contrast',...
                              'callback',@l_ColormapInvert);
tmp=get(H.CONTRAST_INVERT,'position');
% Colormap popup
H.COLMAP_POPUP = uicontrol('parent',H.FIG,...
                           'style','popup',...
                           'units','pixel',...
                           'position',[tmp(1)+tmp(3)+5 tmp(2) 130 25],...
                           'string',{'Gray','Jet','HSV','Cool','Hot','Spring',...
                    'Summer','Autumn','Winter','Bone','Copper','Pink','Custom'},...
                           'enable','on',...
                           'value',1,...
                           'backgroundcolor','w',...
                           'callback',{@l_ChangeColormap,[]});
tmp=get(H.COLMAP_POPUP,'position');
h1=uicontrol('parent',H.FIG,...
             'units','pixel',...
             'position',[tmp(1) tmp(2)+tmp(4)+2 80 12],...
             'style','text',...
             'horizontalalign','left',...
             'backgroundcolor',GD.col.frame,...
             'string','Colormap');







%% Analyze tools ------------------------
tmp = get(H.WINDOWLEVEL_TOOLS,'position');
H.ANALYZE_TOOLS = uicontrol('parent',H.FIG,...H.SIDEBAR,...
'style','frame',...
    'units','pixel',...
    'position',[tmp(1) tmp(2)-318-5 tmp(3) 318],...%[tmp(1) tmp(2)-348-5+40 tmp(3) 308],...
    'backgroundcolor',GD.col.frame);

% Tab-group for uitabs
% $$$ H.ANALYZE_TABGROUP = uitabgroup('parent',H.ANALYZE_TOOLS,...
% $$$                                 'units','normalized',...
% $$$                                 'position',[0 0 1 1],...
% $$$                                 'backgroundcolor',GD.col.frame);
% H.ANALYZE_TABGROUP = uipanel('parent',H.ANALYZE_TOOLS,...
%                              'units','normalized',...
%                              'position',[0 0 1 1],...
%                              'backgroundcolor',GD.col.frame);

% Create uitabs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ H.ROI_UITAB = uitab('parent',H.ANALYZE_TABGROUP,...
% $$$                     'title','ROI tools');
% H.ROI_ANALYSIS = uicontrol('parent',H.FIG,...H.ROI_UITAB,...
%                          'units','normal',...
%                          'position',[0 0 1 1],...
%                          'backgroundcolor',GD.col.frame,...
%                          'bordertype','etchedin');
tmp = get(H.ANALYZE_TOOLS,'position');
h1 = uicontrol('parent',H.FIG,...
               'units','pixel',...
               'position',[10 tmp(2)+tmp(4)-20 200 14],...%[10 332-40 200 14],...
               'style','text',...
               'string','ROI Analysis',...
               'horizontalalign','left',...
               'fontweight','bold',...
               'backgroundcolor',GD.col.frame);

% Radiobutton group for ROI buttons
tmp=get(h1,'position');
H.ROITOOLS_FRAME=uicontrol('parent',H.FIG,...
                           'style','frame',...
                           'units','pixel',...
                           'position',[tmp(1) tmp(2)-55 212 50]);%[10 259-40 212 70]);
tmp=get(H.ROITOOLS_FRAME,'position');



% ROI buttons
tmp=get(H.ROITOOLS_FRAME,'position');
H.ROIBTN_NEW = uicontrol('parent',H.FIG,...
                         'style','pushbutton',...
                         'units','pixel',...
                         'position',[14 tmp(2)+5 40 40],...
                         'string','',...
                         'CData',btn_cdata.cdata.add,...
                         'callback',@l_RoiAdd,...
                         'tooltip','Create a new ROI');

tmp=get(H.ROIBTN_NEW,'position');
tmp_gap=1;
H.ROIBTN_COPY = uicontrol('parent',H.FIG,...
                          'style','pushbutton',...
                          'units','pixel',...
                          'position',[tmp(1)+tmp(3)+tmp_gap tmp(2:4)],...
                          'string','',...
                          'CData',btn_cdata.cdata.copy,...
                          'callback',@l_RoiCopy,...
                          'enable','off',...
                          'tooltip','Copy ROI to slices/volumes');
tmp=get(H.ROIBTN_COPY,'position');
H.ROIBTN_DEL = uicontrol('parent',H.FIG,...
                         'style','pushbutton',...
                         'units','pixel',...
                         'position',[tmp(1)+tmp(3)+tmp_gap tmp(2:4)],...
                         'string','',...
                         'enable','off',...
                         'CData',btn_cdata.cdata.delete,...
                         'callback',{@l_RoiDelete,'selected'},...
                         'tooltip','Delete selected ROI(s)');
tmp=get(H.ROIBTN_DEL,'position');
H.ROIBTN_DELALL = uicontrol('parent',H.FIG,...
                            'style','pushbutton',...
                            'units','pixel',...
                            'position',[tmp(1)+tmp(3)+tmp_gap tmp(2:4)],...
                            'string','',...
                            'enable','off',...
                            'CData',btn_cdata.cdata.deleteall,...
                            'callback',{@l_RoiDelete,'all'},...
                            'tooltip','Delete ALL ROIs');
tmp=get(H.ROIBTN_DELALL,'position');
H.ROIBTN_UNDO = uicontrol('parent',H.FIG,...
                          'style','pushbutton',...
                          'units','pixel',...
                          'position',[tmp(1)+tmp(3)+tmp_gap tmp(2:4)],...
                          'string','',...
                          'fontweight','bold',...
                          'enable','off',...
                          'CData',btn_cdata.cdata.undo,...
                          'tooltip','Undo last draw action',...
                          'callback',@l_RoiUndo);
tmp=get(H.ROITOOLS_FRAME,'position');

H.ROIBTN_LOAD = uicontrol('parent',H.FIG,...
                          'style','pushbutton',...
                          'units','pixel',...
                          'position',[tmp(1) tmp(2)-70 40 40],...
                          'string','',...
                          'fontweight','bold',...
                          'enable','on',...
                          'CData',btn_cdata.cdata.load,...
                          'tooltip','Load ROI(s)',...
                          'callback',{@l_RoiLoad,'interactive'});
tmp=get(H.ROIBTN_LOAD,'position');
H.ROIBTN_SAVE = uicontrol('parent',H.FIG,...
                          'style','pushbutton',...
                          'units','pixel',...
                          'position',[tmp(1) tmp(2)-tmp(4)-5 tmp(3:4)],...
                          'string','',...
                          'fontweight','bold',...
                          'enable','off',...
                          'CData',btn_cdata.cdata.save2,...
                          'tooltip','Save ROI(s)',...
                          'callback',@l_RoiSave);
tmp=get(H.ROIBTN_SAVE,'position');
H.ROIBTN_STATS = uicontrol('parent',H.FIG,...
                          'style','pushbutton',...
                          'units','pixel',...
                          'position',[tmp(1) tmp(2)-tmp(4)-5 tmp(3:4)],...
                          'string','',...
                          'fontweight','bold',...
                          'enable','off',...
                          'CData',btn_cdata.cdata.statistics_large,...
                          'tooltip','View ROI statistics',...
                          'callback',@l_RoiViewStats);
tmp=get(H.ROIBTN_LOAD,'position');
% ROI selection listbox
%uicm=uicontextmenu;
%uimenu(uicm,'Label','Add')
%uimenu(uicm,'Label','testi2')
H.ROIVIEW_LBOX = uicontrol('parent',H.FIG,...
                           'style','listbox',...
                           'units','pixel',...
                           'position',[tmp(1)+tmp(3)+5 tmp(2)-85+40 165 85],...
                           'string','',...
                           'min',0,...
                           'max',2,...
                           'value',[],...
                           'backgroundcolor','w',...
                           'fontweight','bold',...
                           'string','',...
                           'callback',@l_RoiView);
tmp=get(H.ROIVIEW_LBOX,'position');
h1=uicontrol('parent',H.FIG,...
             'units','pixel',...
             'position',[tmp(1) tmp(2)+tmp(4)+2 86 14],...
             'style','text',...
             'string','View ROI(s)',...
             'horizontalalign','left',...
             'fontweight','normal',...
             'backgroundcolor',GD.col.frame);
h1=uicontrol('parent',H.FIG,...
             'units','pixel',...
             'position',[tmp(1) tmp(2)-21 150 14],...
             'style','text',...
             'string','Edit/Current ROI',...
             'horizontalalign','left',...
             'fontweight','normal',...
             'backgroundcolor',GD.col.frame);
tmp=get(h1,'position');

H.ROI_EDIT = uicontrol('parent',H.FIG,...
                       'style','popup',...
                       'units','pixel',...
                       'position',[tmp(1) tmp(2)-24 165 24],...
                       'string',{''},...
                       'enable','off',...
                       'backgroundcolor','w',...
                       'callback',@l_UpdateRoiInfoText);
tmp=get(H.ROI_EDIT,'position');
h1 = uicontrol('parent',H.FIG,...
               'units','pixel',...
               'position',[10 tmp(2)-35 133 14],...
               'style','text',...
               'string','ROI Transparency',...
               'horizontalalign','left',...
               'fontweight','normal',...
               'backgroundcolor',GD.col.frame);
tmp=get(h1,'position');
H.ROI_TRANSP_SLIDER = uicontrol('parent',H.FIG,...
                                'units','pixel',...
                                'position',[10 tmp(2)-18 170 18],...
                                'style','slider',...
                                'min',0,...
                                'max',1,...
                                'value',0.5,...
                                'sliderstep',[0.1 0.3],...
                                'backgroundcolor',GD.col.slider);
if ~Dat.HG2graphics
RoiTransparencyListener = handle.listener(H.ROI_TRANSP_SLIDER,...
                                          'ActionEvent',...
                                          {@l_SetRoiTransparency, ...
                    'slider'});
									set(H.ROI_TRANSP_SLIDER,'userdata',RoiTransparencyListener)
else
	RoiTransparencyListener=addlistener(H.ROI_TRANSP_SLIDER,...
		'ContinuousValueChange',@l_SetRoiTransparency);
	set(H.ROI_TRANSP_SLIDER,'userdata',RoiTransparencyListener)
end
tmp=get(H.ROI_TRANSP_SLIDER,'position');
h1 = uicontrol('parent',H.FIG,...
               'units','pixel',...
               'position',[10 tmp(2)-14 88 14],...
               'style','text',...
               'string','Transparent',...
               'horizontalalign','left',...
               'fontweight','normal',...
               'fontsize',8,...
               'backgroundcolor',GD.col.frame);
h1 = uicontrol('parent',H.FIG,...
               'units','pixel',...
               'position',[tmp(1)+tmp(3)-88 tmp(2)-14 88 14],...
               'style','text',...
               'string','Opaque',...
               'horizontalalign','right',...
               'fontweight','normal',...
               'fontsize',8,...
               'backgroundcolor',GD.col.frame);
H.ROI_TRANSP_EDIT = uicontrol('parent',H.FIG,...
                              'units','pixel',...
                              'position',[tmp(1)+tmp(3)+3 tmp(2) 38 tmp(4)],...
                              'style','edit',...
                              'backgroundcolor','w',...
                              'string','0.5',...
                              'callback',{@l_SetRoiTransparency,'edit'},...
                              'userdata',0.5);


%% Draw invisible uipanel for axes
%tmp=get(H.SIDEBAR,'position');
tmp=get(H.SIDEBAR_FRAME,'position');
% $$$ H.AXES_UIPANEL = uipanel('parent',H.FIG,...
% $$$                          'units','pixel',...
% $$$                          'position',[tmp(1)+tmp(3) ...
% $$$                     0 fig_pos(3)-(tmp(1)+tmp(3)) fig_pos(4)],...
% $$$                          'backgroundcolor',GD.col.mainfig,...
% $$$                          'bordertype','none',...
% $$$                          'tag','ax_uipanel');

%% Draw image axes
gap=10*dxx;
%pos=get(H.SIDEBAR,'position');
pos=get(H.SIDEBAR_FRAME,'position');
% Calculate axis positions
%ax_pos=l_AxesPositions([],[],0,[256,256;256,256;256,256]);

ax_pos = ones(3,4);
H.IMAX1=axes('Parent',H.FIG,...
             'units','pixel',...
             'Position', ...
             ax_pos(1,:),'Ylim',[0 1], ...
             'xlim',[0 1],...
			 'ydir','reverse',...%'Yticklabel',{},... 'Xticklabel',{},...
			 'Box','on', ...
             'xgrid','on',...
             'ygrid','on',...
             'Fontsize',GD.ax_fs,...
             'visible','off');
H.IMOVERLAYAX(1)=axes('Parent',H.FIG,...
                   'units','pixel',...
                   'Position', ...
                   ax_pos(1,:),'Ylim',[0 1], ...
                   'xlim',[0 1],...
				   'ydir','reverse',...
                   'Box','on', ...
                   'xgrid','off',...
                   'ygrid','off',...
                   'Fontsize',GD.ax_fs,...
                   'visible','off',...
                   'hittest','off');

H.IMAX2=axes('Parent',H.FIG,...
             'units','pixel',...
             'Position', ...
             ax_pos(2,:),'Ylim',[0 1], ...
             'xlim',[0 1],...
			 'ydir','reverse',...%'Yticklabel',{},... 'Xticklabel',{},...
			 'Box','on', ...
             'Fontsize',GD.ax_fs,...
             'visible','off');
H.IMOVERLAYAX(2)=axes('Parent',H.FIG,...
                   'units','pixel',...
                   'Position', ...
                   ax_pos(2,:),'Ylim',[0 1], ...
                   'xlim',[0 1],...
				   'ydir','reverse',...
				   'Box','on', ...
                   'xgrid','off',...
                   'ygrid','off',...
                   'Fontsize',GD.ax_fs,...
                   'visible','off',...
                   'hittest','off');

H.IMAX3=axes('Parent',H.FIG,...
             'units','pixel',...
             'Position', ...
             ax_pos(3,:),'Ylim',[0 1], ...
             'xlim',[0 1],...
			 'ydir','reverse',...%'Yticklabel',{},... 'Xticklabel',{},...
			 'Box','on', ...
             'Fontsize',GD.ax_fs,...
             'visible','off');
H.IMOVERLAYAX(3)=axes('Parent',H.FIG,...
                   'units','pixel',...
                   'Position', ...
                   ax_pos(3,:),'Ylim',[0 1], ...
                   'xlim',[0 1],...
				   'ydir','reverse',...
                   'Box','on', ...
                   'xgrid','off',...
                   'ygrid','off',...
                   'Fontsize',GD.ax_fs,...
                   'visible','off',...
                   'hittest','off');

H.COLORBAR_AX = axes('Parent',H.FIG,...
                     'units','pixel',...
                     'Position', ...
                     [1 1 1 1],'Ylim',[0 1], ...
                     'xlim',[0 1],...
                     'yaxislocation','right',...
                     'ydir','reverse',...%'Yticklabel',{},... 'Xticklabel',{},...
                     'Box','on', ...%'fontunits','normal',...'Fontsize',0.023,...
                     'fontsize',8,...
                     'visible','on');

H.COLORBAR_IM = [];
H.IMROI=[];
H.ROIAX=[];
H.IMOVERLAY=[];

% Info text

H.INFOTEXT = text('parent',H.IMAX1,...
  'units','normalized',...
  'position',[0.01 0.99],...
  'verticalalign','top',...
  'horizontalalign','left',...
  'string',{'File: ','Slice: ','Image size: ','Zoom: '},...
  'visible','off',...
	'color','w',...
  'interpreter','none',...
  'fontsize',9,...
  'fontname','fixedwidth',...
  'backgroundcolor','none');
if ~Dat.HG2graphics
	set(H.INFOTEXT,'erasemode','normal')
end
if isunix
  set(H.INFOTEXT,'fontname','bitstream vera sans mono')
end
H.ROIINFOTEXT = text('parent',H.IMAX1,...
  'units','normalized',...
  'position',[0.01 0.01],...
  'verticalalign','bottom',...
  'horizontalalign','left',...
	'color','w',...
  'string',...
  {'Current ROI:',...
  'Mean:',...
  'STD:',...
  'Sum:',...
  'Min:',...
  'Max:',...
  'Pixel Count:'},...
  'visible','off',...
  'interpreter','none',...
  'fontsize',9,...
  'fontname','fixedwidth',...
  'backgroundcolor','none');
if ~Dat.HG2graphics
	set(H.ROIINFOTEXT,'erasemode','normal')
end
if isunix
  set(H.ROIINFOTEXT,'fontname','bitstream vera sans mono')
end

% Set axes to the bottom of the uistack
fig_childs = get(H.FIG,'children');
ind=ismember(double(fig_childs),double([H.IMAX1,H.IMAX2,H.IMAX3])');
fig_childs(end+1:end+3)=fig_childs(ind);
fig_childs(ind)=[];
set(H.FIG,'children',fig_childs)

% Set KeyPressFcn for uicontrol objects
h=findobj(H.FIG,'type','uicontrol');
set(h,'KeyPressFcn',@l_KeyPressFcn);
set(h(strcmpi(get(h,'style'),'edit')),'KeyPressFcn','')



%% Draw frame
H.BLANK_FRAME = uicontrol('Parent',H.FIG,'Units','normalized', ...
			  'Backgroundcol',[0.5 0.5 0.5], ...
			  'Tag','BLANK_FRAME', ...
			  'Position',[0 0 1 1],'Style','frame',...
                          'visible','on');

% Get uicontrol original positions for ResizeFcn -------------
uich=findobj(H.FIG,'type','uicontrol');
uich(find(uich==H.IMSLIDER_FRAME))=[];
uich(find(uich==H.SIDEBAR_FRAME))=[];
uich(find(uich==H.IMSLIDER(1)))=[];
uich(find(uich==H.IMSLIDER(2)))=[];
uich(find(uich==H.BLANK_FRAME))=[];
tmp=get(uich,'position');
H.UICH_ORIG_POS = reshape([tmp{:}],4,length(tmp))';

% Get handles to uicontrols that need to be enabled after initialising
H.UICH_ENABLED=uich(strcmp(get(uich,'enable'),'on'));
H.UICH_ENABLED(end+1)=view_h;
H.UICH_ENABLED(end+1)=edit_h;
H.UICH_ENABLED(end+1)=roi_tools_h;
H.UICH_ENABLED(end+1)=closefile_h;
%H.UICH_ENABLED(end+1)=H.viewprocpar_h;
H.UICH_ENABLED(end+1)=H.UIOVERLAY;
H.UICH_ENABLED(end+1)=H.UIPUSH_SAVE;
%H.UICH_ENABLED(end+1)=H.UIPUSH_PRINT;
H.UICH_ENABLED(end+1)=H.FILEMENU_SAVE_IMAGE;
H.UICH_ENABLED(end+1)=H.FILEMENU_SAVERES;
H.UICH_ENABLED(end+1)=H.UIPUSH_IMSTACK;
H.UICH_ENABLED(end+1)=H.UIPUSH_ROTATE;
H.UICH_ENABLED(end+1)=H.UITOGGLE_VIEW3D;
H.UICH_ENABLED(end+1)=H.UITOGGLE_VIEWX;
H.UICH_ENABLED(end+1)=H.UITOGGLE_VIEWY;
H.UICH_ENABLED(end+1)=H.UITOGGLE_VIEWZ;
H.UICH_ENABLED(end+1)=H.UITOGGLE_GRID;
H.UICH_ENABLED(end+1)=H.UIPLUGINS;
H.UICH_ENABLED(end+1)=H.UITOGGLE_INFO;
H.UICH_ENABLED(end+1)=H.UIPUSH_NEWROI;
H.UICH_ENABLED(end+1)=H.ROIBTN_LOAD;
H.UICH_ENABLED(end+1)=H.FILEMENU_EXPORT;
H.UICH_ENABLED(end+1)=H.UITOGGLE_COLORBAR;
H.UICH_ENABLED(end+1)=H.UITOGGLE_CROSSHAIRS;
H.UICH_ENABLED(end+1)=H.UNFOLD_DATA;

% Get handles to the uicontrols that need to enabled when ROIs are
% defined
H.UICH_ROIENABLED = [H.UIROITOOLS_SAVE,...
                    H.UIROITOOLS_SAVETEMPLATE,...
                    H.UIROITOOLS_COPY,...
                    H.UIROITOOLS_COMP,...
                    H.UIROITOOLS_FLIP,...
                    H.UIROISTATS,...
                    H.UITOGGLE_DRAWROI,...
                    H.UITOGGLE_ERASEROI,...
										H.UITOGGLE_FILLROI,...
                    H.UITOGGLE_MOVEROI,...
                    H.UITOGGLE_MOVEROI3D,...
                    H.UIPUSH_COPYROI,...
                    H.ROIBTN_COPY,...
                    H.ROIBTN_DEL,...
                    H.ROIBTN_DELALL,...
                    H.ROIBTN_SAVE,...
										H.ROIBTN_STATS,...
                    H.ROI_EDIT,...
                    H.UIROISHOWEDGES,...
                    H.UIROIRENAME,...
                    H.UIROISETCOLOR,...
                    H.UIPUSH_STATS_ROI];


% Set enable off for all the uicontrols, when there is no data
set(uich,'enable','off')


end % function H=l_draw_gui()

%%%%%%%%%%%%%%%%%%%%%%
%%
%% Open file(s)
%%
%%%%%%%%%%%%%%%%%%%%%%
function l_OpenFile(h,evd,opt,FileName)
try

% Ask if ROI(s) should be saved
cancel=l_CheckRoiSaved;
if cancel
  return
end

if exist('FileName','var')
  PromptFiles = false;
else
  PromptFiles = true;
end

%Dat.LoadRoiAtStartUp = false;
  
if PromptFiles
  
  % Default filepath
  try
    filepath=getpref('Aedes','GetDataFileDir');
  catch
    filepath='';
  end
  [filefilt,dataformats] = aedes_getfilefilter;
  
  % Ask for a file
  if strcmpi(opt,'single')
    [f_name, f_path, f_index] = uigetfile( ...
      filefilt, ...
      'Select data file',filepath);
    if ( all(f_name==0) | all(f_path==0) ) % Cancel is pressed
      return
		end
		
		% There is a bug in uigetfile in OSX version of Matlab and the
		% fid-directory may be returned instead of the fid-file.
		if ismac
			if length(f_name)>3 && strcmpi(f_name(end-3:end),'.fid')
				f_path = [fullfile(f_path,f_name),filesep];
				f_name = 'fid';
			end
		end
  else
    [f_name, f_path, f_index] = aedes_juigetfiles( ...
      filefilt, ...
      'Select data file(s)',filepath);
    if ~iscell(f_name) && ~iscell(f_path)
      if ( all(f_name==0) | all(f_path==0) ) % Cancel is pressed
        return
      end
    end
  end
  
  % Put strings in a cell array so we don't have to check it all the
  % time...
  if ~iscell(f_name)
    f_path = {f_path};
    f_name = {f_name};
  end
  
  % Set current directory to preferences
  setpref('Aedes','GetDataFileDir',f_path{1})
  
else
  if strcmpi(opt,'single')
	[f_path,f_name,f_ext]=fileparts(FileName);
	if isempty(f_path)
	  f_path = {[pwd,filesep]};
	else
	  if ispc && length(f_path)<3
		f_path = {[pwd,filesep,f_path,filesep]};
	  elseif isunix && f_path(1)~=filesep
		f_path = {[pwd,filesep,f_path,filesep]};
	  else
		f_path = {[f_path,filesep]};
	  end
	end
	f_name = {[f_name,f_ext]};
% 	f_index = 0;
% 	dataformat = '';
	
	% Check that the file exists...
	if exist([f_path{1},f_name{1}],'file')~=2
	  showError = true;
	  
	  % If the file cannot be found, it might be a VNMR file with only
	  % folder reference...
	  if strcmpi(f_ext,'.fid')
		if exist([f_path{1},f_name{1},filesep,'fid'],'file')~=2
		  showError = true;
		else
		  showError = false;
		  f_path = {[f_path{1},f_name{1},filesep]};
		  f_name = {'fid'};
		end
	  end
	  
	  if showError
		h=errordlg({'Could not find file',...
		  ['"',f_path{1},f_name{1},'"'],...
		  '','Check that the file has not been removed and can be accessed.'},...
		  'Cannot not find file', ...
		  'modal');
		uiwait(h);
		return
	  end
	end
  else
	f_name = {};
	f_path = {};
	%f_index = 0;
	
	%% Extract file and path names to separate cell arrays
	for ii=1:length(FileName)
	  [fp,fn,fe]=fileparts(FileName{ii});
	  f_path{ii} = [fp,filesep];
	  f_name{ii} = [fn,fe];
	end
	%dataformat = '';
  end
end

% Close file if necessary
if ~isempty(DATA)
  l_CloseFile([],[]);
end


%% Read the data from other data files
if strcmpi(opt,'single') || length(f_name)==1
  try
    DATA=aedes_data_read([f_path{1},f_name{1}]);
  catch
    h=errordlg({['Could not read file "',[f_path{1},f_name{1}],'"'],'',lasterr},...
      'Error reading file', ...
      'modal');
    uiwait(h);
    return
  end
  % Set to Recent Files menu
  l_BuildRecentFilesMenu([],[],[],[f_path{1},f_name{1}]);
else
  DATA={};
  count=0;
  SkippedFiles={};
  SkippedInd = [];
  [h,txh]=aedes_calc_wait({['Reading file 1/',...
	num2str(length(f_name))],'""'});
  for ii=1:length(f_name)
    set(txh,'String',...
      sprintf('%s\n%s',['Reading file ' num2str(ii) '/' num2str(length(f_name))],...
      ['"',f_path{ii},f_name{ii},'"']))
    drawnow
    try
      tmp=aedes_data_read([f_path{ii},f_name{ii}],'wbar','off');
      
      % Determine from the first file if we are reading one slice data into
      % "mixed" cell arrayed format or 3D volumes into a 4D file
      if ii==1
        if ndims(tmp.FTDATA)==3
          % Read 3D data into 4D array
          Read2Mixed = false;
          
          % Reserve space for data
          DATA = tmp;
          DATA.FTDATA = zeros(size(tmp.FTDATA,1),...
            size(tmp.FTDATA,2),size(tmp.FTDATA,3),...
            length(f_name),class(tmp.FTDATA));
          DATA.FTDATA(:,:,:,1)=tmp.FTDATA;
        elseif ndims(tmp.FTDATA)==2
          % Read 2D slices into mixed cell array...
          Read2Mixed = true;
          DATA = {tmp};
        else
          % Throw error...
          Read2Mixed = false;
          error('Reading of multiple 4D files is not supported.')
        end
      end
      if ~Read2Mixed
        if isequal(size(DATA.FTDATA(:,:,:,1)),size(tmp.FTDATA))
          DATA.FTDATA(:,:,:,ii) = tmp.FTDATA;
        else
          count=count+1;
          SkippedFiles{count} = [f_path{ii},f_name{ii}];
          disp(['Aedes: Warning: Data size in file "',f_path{ii},...
            f_name{ii},'" does not match with the first data file! Skipping file...'])
        end
      elseif Read2Mixed
        DATA{ii}=tmp;
        if length(size(tmp.FTDATA))>2
          count=count+1;
          SkippedInd(end+1) = ii;
          SkippedFiles{count} = [f_path{ii},f_name{ii}];
          DATA{ii}=[];
          disp(['Aedes: Warning: File "',f_path{ii},...
            f_name{ii},'" contains multiple slices! Skipping file...'])
        end
      end
    catch
      if ii==1
        h=errordlg({'Could not read file',...
          ['"',f_path{ii},f_name{ii},'"'],...
          '','Returned error:',lasterr},...
          'Could not read file','modal');
        DATA = [];
        break
      elseif Read2Mixed
        count=count+1;
        SkippedFiles{count} = [f_path{ii},f_name{ii}];
        DATA{ii}=[];
        disp(['Aedes: Warning: Could not read file "',...
          f_path{ii},f_name{ii},'". Skipping file...'])
      else
        % Throw the actual error and break
        h=errordlg({'Could not read file',...
          ['"',f_path{ii},f_name{ii},'"'],...
          '','Returned error:',lasterr},...
          'Could not read file','modal');
        DATA = [];
        break
      end
    end
  end
  delete(h)
  
  % Ensure that there are no empty values in DATA
  if iscell(DATA)
    empty_val=cellfun('isempty',DATA);
    if all(empty_val)
      DATA=[];
    else
      DATA(empty_val)=[];
    end
  end

	% If the data consists of multiple 2D files that are all of the same
	% size, ask if the user wants to view them as 3D data.
	if Read2Mixed
		nFiles = length(DATA);
		X = zeros(1,nFiles);
		Y = zeros(1,nFiles);
		for ii = 1:nFiles
			X(ii)=size(DATA{ii}.FTDATA,1);
			Y(ii)=size(DATA{ii}.FTDATA,2);
		end
		if all(X==X(1)) && all(Y==Y(1))
			resp = questdlg('Do you want to view data as a 3D data set or a 2D image stack?',...
				'View data 3D or 2D?','View as 3D data','View as 2D image stack','View as 2D image stack');
			if isempty(resp)
				viewAs3D = false;
			elseif strcmpi(resp,'View as 3D data')
				viewAs3D = true;
			else
				viewAs3D = false;
			end
			if viewAs3D
				tmp = DATA{1};
				tmp.FTDATA = zeros(X(1),Y(1),nFiles,class(DATA{1}.FTDATA));
				for ii = 1:nFiles
					tmp.FTDATA(:,:,ii) = DATA{ii}.FTDATA;
					DATA{ii}.FTDATA=[];
				end
				DATA=[];
				DATA=tmp;
			end
		end
	end
  
  % Remove skipped volumes
  if ~Read2Mixed && ~isempty(SkippedInd)
    DATA.FTDATA(:,:,:,SkippedInd) = [];
  end
  
  % Warn about skipped files
  if count~=0
    if isempty(DATA)
      h=errordlg(['All files were skipped either due to errors ',...
        'or because they contained multiple slices.'],...
        'Could not read files','modal');
      return
    else
      h=warndlg({['The following file(s) were ',...
        'skipped either due to errors or ',...
        'because they contained multiple slices:'],'',SkippedFiles{:}},...
        'Some Files Were Skipped','modal');
      uiwait(h)
    end
  end
end

if isempty(DATA)
  % Canceled
  return
end

if isstruct(DATA)
  DATA = {DATA};
end

%% Load ROI also if they can be found in the DATA structure
if isfield(DATA{1},'ROI')
  ROI = DATA{1}.ROI;
  DATA{1}=rmfield(DATA{1},'ROI');
end

%% If the DATA was loaded from MAT-File, it might also contain rotation
%% and flipping and clim information
if isfield(DATA{1},'DataRotation') && isfield(DATA{1},'DataFlip')
  Dat.DataRotation = DATA{1}.DataRotation;
  Dat.DataFlip = DATA{1}.DataFlip;
  DATA{1}=rmfield(DATA{1},'DataRotation');
  DATA{1}=rmfield(DATA{1},'DataFlip');
end
if isfield(DATA{1},'RotateFlip3d')
  Dat.RotateFlip3d = DATA{1}.RotateFlip3d;
  DATA{1}=rmfield(DATA{1},'RotateFlip3d');
end

% Check if data file contain only kspace information
if length(DATA)==1 && isfield(DATA{1},'KSPACE') && ...
		(isempty(DATA{1}.FTDATA) & ~isempty(DATA{1}.KSPACE) )
	drawnow;
	resp=questdlg(['The inputted DATA structure contains only k-space',...
		' information. Do you want to view real, imaginary, or absolute',...
		' part of the complex k-space?'],...
		'Complex data inputted',...
		'Real','Imaginary','Absolute','Absolute');
	if strcmpi(resp,'Real')
		DATA{1}.FTDATA = real(DATA{1}.KSPACE);
	elseif strcmpi(resp,'Imaginary')
		DATA{1}.FTDATA = imag(DATA{1}.KSPACE);
	elseif strcmpi(resp,'Absolute')
		DATA{1}.FTDATA = abs(DATA{1}.KSPACE);
	else
		DATA=[];
		return
	end
end

% Initialize GUI
l_Initialize([],[]);

if isfield(DATA{1},'SliceClim')
  Dat.SliceClim = DATA{1}.SliceClim;
  DATA{1}=rmfield(DATA{1},'SliceClim');
	set(H.CLIMRANGE_POPUP,'value',2)
	l_ChangeClimRange([],[])
end

%% If the ROI-structure has been given in file load time, refresh ROIs
if ~isempty(ROI)
  l_RoiLoad([],[],'')
end


catch
  aedes_errordump(lasterror);
end
end % function l_OpenFile(h,


%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Save Image Data
%%
%%%%%%%%%%%%%%%%%%%%%%%%%
function l_SaveImageData(h,evd)
try
%% Get default save path
try
  filepath=getpref('Aedes','PutDataFileDir');
catch
  filepath='';
end


%% Save 3D data
if length(DATA)==1

  %% Get default save file name
  try
    if isfield(DATA{Dat.DataInd},'DataFormat') && ...
        strcmpi(DATA{Dat.DataInd}.DataFormat,'vnmr')
      [fp,fn,fe]=fileparts([DATA{Dat.DataInd}.HDR.fpath,DATA{Dat.DataInd}.HDR.fname]);
      [fp,fn,fe]=fileparts(fp);
    else
      [fp,fn,fe]=fileparts([DATA{Dat.DataInd}.HDR.fpath,DATA{Dat.DataInd}.HDR.fname]);
    end
  catch
    fp=[pwd,filesep];
    fn='';
    fe='';
  end
  
  %% Get saved file name and type
  [fname, fpath, findex] = uiputfile( ...
    {'*.nii;*.NII','NIfTI-files (*.nii)'; ...
     '*.hdr;*.HDR','Analyze 7.5 files (*.hdr)';...
     '*.mat:*.MAT','Matlab MAT-File (*.mat)';...
     '*.mri;*.MRI','MRI-files (*.mri)'}, ...
    'Save Image Data as', [filepath,fn]);
  if ( all(fname==0) | all(fpath==0) ) % Cancel is pressed
    return
  end
  
  %% Set fpath to preferences
  setpref('Aedes','PutDataFileDir',fpath);
  
  % Write NIfTI-files -------------------------
  if findex==1
    
    %% Prompt save type
    resp=questdlg('Save slices in separate NIfTI-files or in a single file?',...
                  'Separate files or single file?',...
                  'Separate files','Single file','Cancel',...
                  'Separate files');
    if isempty(resp) || strcmpi(resp,'Cancel')
      return
    elseif strcmpi(resp,'Separate files')
      singleFile = false;
    else
      singleFile = true;
    end
    fext = '.nii';
    fileType = 'NIfTI';
    
    
  elseif findex==2 % Write Analyze 7.5 file(s) ----------------
    
     %% Prompt save type
    resp=questdlg('Save slices in separate Analyze75-files or in a single file?',...
                  'Separate files or single file?',...
                  'Separate files','Single file','Cancel',...
                  'Separate files');
    if isempty(resp) || strcmpi(resp,'Cancel')
      return
    elseif strcmpi(resp,'Separate files')
      singleFile = false;
    else
      singleFile = true;
    end
    fext = '.hdr';
    fileType = 'Analyze 7.5';
    
  elseif findex==3 % Save data in MATLAB MAT-File ------------------
    
    [fp,fn,fe]=fileparts([fpath,fname]);
    try
      SliceClim = Dat.SliceClim;
      if strcmpi(fe,'.mat')
        save([fpath,fname],'DATA','SliceClim','-mat')
      else
        save([fp,filesep,fn,'.mat'],'DATA','SliceClim','-mat')
      end
    catch
      errordlg({'Could not save MAT-file. The following error was returned:',...
                '',lasterr},'Error saving MAT-file!','modal')
    end
    return
    
  else % Save data in MRI-files -------------------
    
    singleFile = false;
    
    %% Warn about MRI-files
    resp=questdlg({['WARNING! Exporting data to MRI-format is not advisable. ' ...
                    'The data will be scaled and rounded to 12-bit before ' ...
                    'writing into the MRI-file, and scaling errors of ' ...
                    'various severities will occur!'],'',...
                   'Are you sure you want to continue?'},'WARNING!',...
                  'Yes','No','No');
    if isempty(resp) || strcmpi(resp,'No')
      return
    end
    
    fileType = 'MRI';
    fext = '.mri';
  end
  
  % Determine file prefix
  [fp,fn,fe]=fileparts([fpath,fname]);
  if strcmpi(fe,fext)
    fprefix = fn;
  else
    fprefix = [fn,fe];
  end
  
  done=false;
  
  %% Write separate files
  if ~singleFile
    
    % Create file names
    filenames={};
    for ii=1:Dat.ImageDim(3)
      filenames{ii}=sprintf('%s%03d%s',[fprefix,'_'],ii,fext);
    end
    nFiles = length(filenames);
    
    % Check if files exist
    ind=aedes_check_file_exist(filenames,fpath);
    if any(ind)
      % Warn if some files are about to be overwritten
      overwrited_files = {filenames{ind}};
      
      % Limit the file list to 20 files
      if length(overwrited_files)>20
        overwrited_files = {overwrited_files{1:20}};
        overwrited_files{end+1}='...';
      end
      
      resp=questdlg({['The following files already exist in the output directory' ...
                      ' "',fpath,'"'],'','(NOTE: max. 20 files shown here)','',...
                     overwrited_files{:},'',...
                     'Do you want to overwrite these files?'},...
                    'Overwrite Existing Files?','Overwrite','Abort','Abort');
      if isempty(resp) || strcmpi(resp,'Abort')
        return
      end
    end
    
    if strcmpi(fileType,'MRI')
      done=true;
      h=aedes_wbar(0,sprintf('Saving slice 1/%d in MRI format...',nFiles));
      for ii=1:nFiles
        aedes_wbar(ii/nFiles,h,sprintf('Saving slice %d/%d in MRI format...',ii, ...
                                 nFiles))
        try
		  aedes_smiswrite(DATA{Dat.DataInd}.FTDATA(:,:,ii),...
			[fpath,filenames{ii}])
        catch
          done=false;
          break
        end
      end
      
    elseif strcmpi(fileType,'NIfTI')
      h=aedes_wbar(0,sprintf('Saving slice 1/%d in NIfTI format...',nFiles));
      for ii=1:nFiles
        aedes_wbar(ii/nFiles,h,sprintf('Saving slice %d/%d in NIfTI format...',ii, ...
                                 nFiles))
        [done,msg]=aedes_write_nifti(DATA{Dat.DataInd}.FTDATA(:,:,ii),...
                               [fpath,filenames{ii}],'DataType','single',...
                               'FileType',2,'Clim',Dat.SliceClim(Dat.DataInd,:));
      end
      
    elseif  strcmpi(fileType,'Analyze 7.5')
      h=aedes_wbar(0,sprintf('Saving slice 1/%d in Analyze 7.5 format...',nFiles));
      for ii=1:nFiles
        aedes_wbar(ii/nFiles,h,sprintf('Saving slice %d/%d in Analyze 7.5 format...',ii, ...
                                 nFiles))
        [done,msg]=aedes_write_nifti(DATA{Dat.DataInd}.FTDATA(:,:,ii),...
                               [fpath,filenames{ii}],'DataType','single',...
                               'FileType',0,'Clim',Dat.SliceClim(Dat.DataInd,:));
      end
    end
    
  else
    if strcmpi(fileType,'NIfTI')
      [h,txh]=aedes_calc_wait('Saving image data in NIfTI format...');
      [done,msg]=aedes_write_nifti(DATA{Dat.DataInd}.FTDATA,...
                             [fpath,fname],'filetype',2,...
                             'Clim',Dat.SliceClim(Dat.DataInd,:));
    elseif  strcmpi(fileType,'Analyze 7.5')
      [h,txh]=aedes_calc_wait('Saving image data in Analyze 7.5 format...');
      [done,msg]=aedes_write_nifti(DATA{Dat.DataInd}.FTDATA,...
                             [fpath,fname],'filetype',0,...
                             'Clim',Dat.SliceClim(Dat.DataInd,:));
    end
  end
  delete(h)
  
  %% Show warning if something went wrong
  if ~done
    h=errordlg(msg,'Error while saving image data','modal')
    uiwait(h);
    return
  end
  

  %% Save Slice Data
else
  
  %% Get default save file name
  if strcmpi(DATA{Dat.DataInd}.DataFormat,'vnmr')
    [fp,fn,fe]=fileparts([DATA{Dat.DataInd}.HDR.fpath,DATA{Dat.DataInd}.HDR.fname]);
    [fp,fn,fe]=fileparts(fp);
  else
    [fp,fn,fe]=fileparts([DATA{Dat.DataInd}.HDR.fpath,DATA{Dat.DataInd}.HDR.fname]);
  end
  
  %% Get saved file name and type
  [fname, fpath, findex] = uiputfile( ...
      {'*.mat:*.MAT','Matlab MAT-File (*.mat)';...
       '*.nii;*.NII','NIfTI-files (*.nii)'; ...
       '*.hdr;*.HDR','Analyze 7.5 files (*.hdr)'}, ...
    'Save Image Data as', [filepath,fn]);
  if ( all(fname==0) | all(fpath==0) ) % Cancel is pressed
    return
  end
  
  %% Set fpath to preferences
  setpref('Aedes','PutDataFileDir',fpath);
  
  if findex==1 %% Save data in MAT-file
    
    % Determine file prefix
    [fp,fn,fe]=fileparts([fpath,fname]);
    if strcmpi(fe,'.mat')
      fprefix = fn;
    else
      fprefix = [fn,fe];
    end
    
    %% Save data in Matlab MAT-file
    try
      DataRotation = Dat.DataRotation;
      DataFlip = Dat.DataFlip;
      SliceClim = Dat.SliceClim;
      save([fpath,fprefix,'.mat'],'DATA','DataRotation','DataFlip',...
           'SliceClim','-mat')
    catch
      hh=errodlg({'Error occurred while saving DATA to file',...
                 [fpath,fprefix,'.mat']},'Error.','modal');
      return
    end
    
  elseif any(findex==[2 3])
    if findex==2
      fext = '.nii';
      fileType = 2;
    else
      fext = '.hdr';
      fileType = 0;
    end
    if findex==2 && Dat.DataL>1
      %% Ask if NIfTI files are saved in separate files or the whole
      %% volume into a single file
      resp=questdlg('Save slices in separate NIfTI-files or in a single file?',...
                    'Separate files or single file?',...
                    'Separate files','Single file','Cancel',...
                    'Separate files');
      if isempty(resp) || strcmpi(resp,'Cancel')
        return
      elseif strcmpi(resp,'Separate files')
        singleFile = false;
      else
        singleFile = true;
      end
    else
      % Try to write Analyze75 files allways into a single volume file
      singleFile = true;
    end
    
    
    % Determine file prefix
    [fp,fn,fe]=fileparts([fpath,fname]);
    if strcmpi(fe,fext)
      fprefix = fn;
    else
      fprefix = [fn,fe];
    end
    
    if singleFile
      %% Check if all slices are the same size. Otherwise it is not
      %% possible to save all the slices in one file
      if all(Dat.ImageDim(:,1)/Dat.ImageDim(1,1)) && ...
          all(Dat.ImageDim(:,2)/Dat.ImageDim(1,2))
        
        [calc_h,txh]=aedes_calc_wait({'Writing image data to file',
                            [fpath,fprefix,fext]});
        
        % Generate a new data matrix. This can in extreme cases inflict
        % "out of memory" error.
        data_mtx = zeros(Dat.ImageDim(1,1),Dat.ImageDim(1,2),Dat.DataL,...
                         'single');
        for ii=1:Dat.DataL
          data_mtx(:,:,ii)=DATA{ii}.FTDATA;
        end
        
        % Write data to NIfTI/Analyze75-file
        [done,msg]=aedes_write_nifti(data_mtx,[fpath,fprefix,fext],...
                               'FileType',fileType,'DataType', ...
                               'single');
        clear('data_mtx')
        if ~done
          delete(calc_h)
          hh=errordlg({'An error occurred while writing file',...
                       [fpath,fprefix,fext]},'Error Writing File',...
                      'modal');
          return
        end
        delete(calc_h)
      else
        hh=errordlg(['Cannot write slices into a single file, because all ' ...
                     'the slices are not of the same size.'],...
                    'Cannot Write Single File','modal');
        return
      end
    else
      %% Write all slices into separate NIfTI-files
      % Create file names
      filenames={};
      for ii=1:Dat.DataL
        filenames{ii}=sprintf('%s%03d%s',[fprefix,'_'],ii,fext);
      end
      nFiles = length(filenames);
      
      % Check if files exist
      ind=aedes_check_file_exist(filenames,fpath);
      if any(ind)
        % Warn if some files are about to be overwritten
        overwrited_files = {filenames{ind}};
        
        % Limit the file list to 20 files
        if length(overwrited_files)>20
          overwrited_files = {overwrited_files{1:20}};
          overwrited_files{end+1}='...';
        end
        
        resp=questdlg({['The following files already exist in the output directory' ...
                        ' "',fpath,'"'],'','(NOTE: max. 20 files shown here)','',...
                       overwrited_files{:},'',...
                       'Do you want to overwrite these files?'},...
                      'Overwrite Existing Files?','Overwrite','Abort','Abort');
        if isempty(resp) || strcmpi(resp,'Abort')
          return
        end
      end
      
      h=aedes_wbar(0,sprintf('Saving slice 1/%d in NIfTI format...',nFiles));
      for ii=1:nFiles
        aedes_wbar(ii/nFiles,h,sprintf('Saving slice %d/%d in NIfTI format...',ii, ...
                                 nFiles))
        [done,msg]=aedes_write_nifti(DATA{ii}.FTDATA,...
                               [fpath,filenames{ii}],'DataType','single',...
                               'FileType',2,...
                               'Clim',Dat.SliceClim(ii,:));
        if ~done
          delete(h)
          hh=errordlg({'An error occurred while writing file',...
                      [fpath,filenames{ii}]},'Error Writing File', ...
                      'modal');
          return
        end
      end
      delete(h)
    end
    
  end
end

catch
  aedes_errordump(lasterror);
end
end % function l_SaveImageData(h,


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check if ROI(s) are saved
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cancel=l_CheckRoiSaved()

% Warn if ROIs are not saved
if ~isempty(ROI) && ~Dat.RoiSaved
  resp=questdlg({'Modified ROI(s) exist.',...
                 'Do you want save ROI(s) before closing file?'},...
                'Save ROI(s) before closing file?','Yes','No','Cancel','Yes');
  if strcmpi(resp,'Yes')
    l_RoiSave([],[])
    cancel = false;
    return
  elseif isempty(resp) || strcmpi(resp,'Cancel')
    cancel = true;
    return
  else
    cancel = false;
  end
else
  cancel = false;
end
figure(H.FIG)

end % function cancel=l_CheckRoiSaved()


%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Close File(s)
%%
%%%%%%%%%%%%%%%%%%%%%%%%
function l_CloseFile(h,evd,opt)
try

  if nargin>2 && strcmpi(opt,'PreserveData')
    PreserveData = true;
  else
    PreserveData = false;
  end

  % Warn if ROIs are not saved
  if ~isempty(h) && ishandle(h) && ~PreserveData
    cancel=l_CheckRoiSaved;
    if cancel
      return
    end
  end

% Reset uicontrol focus
%uicontrol(H.SL_SLIDER);

% Remember HG2 status
hg2_status = Dat.HG2graphics;

% Delete Voxel Time Series Figure if open
try
  close(H.TSFIG)
end
set(H.UIVIEW_TIMESERIES,'checked','off')

%% Remove callback for Mouse scroll wheel
if Dat.MatlabVersion>=7.05
  set(H.FIG,'WindowScrollWheelFcn','')
end

% Delete all ROIs and related images and axes
if ~PreserveData
  clear ROI DATA Dat
  ROI = [];
else
  clear Dat
end

try
  delete(H.IMROI)
end
try
  delete(H.ROIAX)
end
H.ROIAX=[];
H.IMROI=[];
if isfield(H,'IM')
  set(H.IM,'CData',[])
end

% Set mouse for changing window/level
set(get(H.UITOGGLE_ARROW,'userdata'),'state','off')
set(H.UITOGGLE_ARROW,'state','on')

% Clear ROI undo buffer and disable ROI uicontrols
l_RoiUndoBuffer('clear')
set(H.ROIVIEW_LBOX,'value',[],'string','')
set([H.ROIBTN_DEL,H.ROIBTN_DELALL,...
     H.ROIBTN_UNDO,H.UIPUSH_UNDOROI,H.ROIBTN_COPY],'enable','off')
set(H.ROI_EDIT,'value',1,'string',{''},'enable','off')

% Delete crossbars
try
	delete(H.CROSSBAR_LN)
end
H.CROSSBAR_LN=[];
%set(H.SHOW_CROSSBARS,'value',0)

% Reset image sliders
set(H.IMSLIDER,'visible','off','value',0)
set(H.IMSLIDER_FRAME,'visible','off')

% Reset uitoolbar buttons
set([H.UIPUSH_ZOOMIN,H.UIPUSH_ZOOMOUT,H.UITOGGLE_ZOOMNORM],'enable','off')
set([H.UITOGGLE_VIEW3D,H.UITOGGLE_VIEWY,H.UITOGGLE_VIEWZ,H.UITOGGLE_CROSSHAIRS],'state','off')
set([H.UITOGGLE_ARROW,H.UITOGGLE_PAN,...
     H.UITOGGLE_DRAWROI,H.UITOGGLE_ERASEROI,H.UITOGGLE_MOVEROI,...
     H.UITOGGLE_MOVEROI3D,H.UITOGGLE_GRID,H.UITOGGLE_CROSSHAIRS],...
    'enable','off')
set([H.UIPUSH_NEWROI,H.UIPUSH_UNDOROI],...
    'enable','off')
set(H.UICH_ROIENABLED,'enable','off')
set(H.UIMENU_ROTATEFLIP,'enable','off')
set(H.UIEDIT_IMSTACK,'enable','off')

% Delete overlay control window
if isfield(H,'OVERLAY_CONTROL_FIG') && ~isempty(H.OVERLAY_CONTROL_FIG) ...
    && ishandle(H.OVERLAY_CONTROL_FIG)
  delete(H.OVERLAY_CONTROL_FIG)
  H.OVERLAY_CONTROL_FIG = [];
end

% Disable overlay uimenus
set([H.UIOVERLAY_CONTROLS,...
	H.UIOVERLAY_SAVE,...
	H.UIOVERLAY_DELETE],'enable','off')

% Reset volume button and editbox
set([H.V_BTN,H.V_EDIT],'enable','off')
set(H.V_EDIT,'string','1')

% Delete overlays if loaded
if isfield(H,'IMOVERLAY') && ~isempty(H.IMOVERLAY)
  try
    delete(H.IMOVERLAY)
  catch
  end
  H.IMOVERLAY=[];
end

% Delete colorbar image and intensity line
try
  delete(H.COLORBAR_IM)
end
try
  delete(H.COLORBAR_LN)
end
try
  delete(H.COLORBAR_TX)
end
H.COLORBAR_IM=[];
H.COLORBAR_LN=[];
H.COLORBAR_TX=[];

% Delete data images
try
  delete(H.IM)
end

% Clear data from memory
if ~PreserveData
  DATA = [];
end

% disable uicontrols
set(H.UICH_ENABLED,'enable','off')

% Reset figure window title
set(H.FIG,'Name','Aedes 1.0')

% Update Voxel Value text
set(H.VOXEL_VALUE,'string','-')

%% Place blank frame in front of gui
set(H.BLANK_FRAME,'visible','on')
%drawnow

% Clear internal variables
Dat = [];
Dat.HG2graphics = hg2_status;

catch
  aedes_errordump(lasterror);
end
end % function l_CloseFile(h,


%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Initialize GUI
%%
%%%%%%%%%%%%%%%%%%%%%%%
function l_Initialize(h,evd)
try
  
%% Initialize "Mixed" type (i.e. slice per file) data -------------------
if length(DATA)>1 || ( length(DATA)==1 & any(ndims(DATA{1}.FTDATA)==[1 2]))
  Dat.isDataMixed = true;
  Dat.DataL = length(DATA);
  Dat.RealAxSize = cell(1,Dat.DataL);
  for ii=1:Dat.DataL
    Dat.ImageDim(ii,:) = [size(DATA{ii}.FTDATA) 1 1];
    Dat.RealAxSize{ii} = [Dat.ImageDim(ii,2),Dat.ImageDim(ii,1);...
                        Dat.ImageDim(ii,3),Dat.ImageDim(ii,1);...
                        Dat.ImageDim(ii,3),Dat.ImageDim(ii,2)];
  end
  Dat.DataInd = 1;
  if ~isfield(Dat,'DataRotation') || ~isfield(Dat,'DataFlip')
    Dat.DataRotation = zeros(1,Dat.DataL);
    Dat.DataFlip = zeros(1,Dat.DataL);
  end
  Dat.RotateFlip3d = {};
else %% Initialize normal 3D- or 4D data --------------
  Dat.isDataMixed = false;
  Dat.ImageDim = size(DATA{1}.FTDATA);
  if length(Dat.ImageDim)<4
    le=length(Dat.ImageDim);
    Dat.ImageDim = [Dat.ImageDim ones(1,4-le)];
  end
  Dat.DataL = 1;
  Dat.DataInd = 1;
  if ~isfield(Dat,'DataRotation') || ~isfield(Dat,'DataFlip')
    Dat.DataRotation = [];
    Dat.DataFlip = [];
  end
  if ~isfield(Dat,'RotateFlip3d')
	Dat.RotateFlip3d = {};
  end
  Dat.RealAxSize{Dat.DataInd} = [Dat.ImageDim(Dat.DataInd,2),Dat.ImageDim(Dat.DataInd,1);...
                      Dat.ImageDim(Dat.DataInd,3),Dat.ImageDim(Dat.DataInd,1);...
                      Dat.ImageDim(Dat.DataInd,3),Dat.ImageDim(Dat.DataInd,2)];
end

% % Detect if HG2 graphics are used
% if isa(H.FIG,'matlab.ui.Figure')
% 	Dat.HG2graphics = true;
% else
% 	Dat.HG2graphics = false;
% 	set(H.FIG,'DoubleBuffer','on')
% end




% Save data file headers -------------------------------
for ii=1:Dat.DataL
  Dat.HDR{ii} = DATA{ii}.HDR;
end

% Initialize volumes and slices ------------------------
Dat.Vols = Dat.ImageDim(Dat.DataInd,4);
Dat.Rcvrs = size(DATA{1}.FTDATA,5);
Dat.CurrentVol = 1;
% Dat.Slices = [1 1 1];
Dat.Slices = [ceil(Dat.ImageDim(1)/2), ceil(Dat.ImageDim(2)/2), ceil(Dat.ImageDim(3)/2)];

% Initialize Clim --------------------------------------
Dat.ScaleMin=[];
Dat.ScaleMax=[];
Dat.OrigClim=[];
for ii=1:Dat.DataL
  Dat.ScaleMin(ii) = double(min(min(min(DATA{ii}.FTDATA(:,:,:,Dat.CurrentVol)))));
  Dat.ScaleMax(ii) = double(max(max(max(DATA{ii}.FTDATA(:,:,:,Dat.CurrentVol)))));
  Dat.OrigClim(ii,:) = [Dat.ScaleMin(ii) Dat.ScaleMax(ii)];
  if diff(Dat.OrigClim(ii,:))==0
    Dat.OrigClim(ii,:) = [0 1];
  end
end

% Detect Matlab version
[Dat.MatlabVersion,Dat.isImageProc] = aedes_getmatlabversion;

% If any of the currently opened files is Analyze or NIfTI type, the
% cal_max and cal_min fields are checked and adjusted accordingly
Dat.Clim = Dat.OrigClim(1,:);
if ~isfield(Dat,'SliceClim')
  Dat.SliceClim = Dat.OrigClim;
  
  for ii=1:Dat.DataL
    if isfield(DATA{ii},'DataFormat') && ...
        any(strcmpi(DATA{ii}.DataFormat,{'analyze75','nifti(1)', ...
                          'nifti(2)'}))
      try
        cal_min=DATA{ii}.HDR.FileHeader.dime.cal_min;
        cal_max=DATA{ii}.HDR.FileHeader.dime.cal_max;
        if cal_min~=0
          Dat.SliceClim(ii,1)=cal_min;
        end
        if cal_max~=0
          Dat.SliceClim(ii,2)=cal_max;
        end
      catch
        continue
      end
    end
  end
end

% Get default Clim behavior from preferences -------------
% try
%   Dat.LockClim = getpref('Aedes','LockClim');
%   if Dat.LockClim==0
%     val=3;
%   elseif Dat.LockClim==1
%     val=1;
%   elseif Dat.LockClim==2
%     val=2;
%   else
%     val=4;
%   end
% catch
  Dat.LockClim = 0; % 0 = clim unlocked
                    % 1 = clim locked global
                    % 2 = clim locked slice
                    % 3 = clim unlocked and auto-balanced
  val = 3;
% end

% Set Clim values ----------------------------------
set(H.CLIMMAXEDIT,'string',num2str(Dat.Clim(2)),...
                  'userdata',Dat.Clim(2))
set(H.CLIMMINEDIT,'string',num2str(Dat.Clim(1)),...
                  'userdata',Dat.Clim(1))
set(H.CLIMRANGE_POPUP,'value',val)

% Default gap between axes --------------------
Dat.AxGap=4;

% Default view direction --------------------
Dat.AxView=0;

% Initialize colorbar -----------------------
Dat.ColMap = [];
Dat.ColMapInverted=false;
try
  Dat.ShowColorbar = getpref('Aedes','ShowColorbar');
  if Dat.ShowColorbar
    set(H.UITOGGLE_COLORBAR,'state','on')
  end
catch
  Dat.ShowColorbar = false;
  set(H.UITOGGLE_COLORBAR,'state','off')
end

% Try to determine field of view (FOV) ----------------------------
if isfield(DATA{Dat.DataInd},'PROCPAR') && isfield(DATA{Dat.DataInd}.PROCPAR,'lro') && ...
      isfield(DATA{Dat.DataInd}.PROCPAR,'lpe')
  Dat.FOV = [DATA{Dat.DataInd}.PROCPAR.lpe DATA{Dat.DataInd}.PROCPAR.lro];
else
  Dat.FOV = [];
end

% Get default Zoom levels from preferences -----------------------------
Dat.OldZoom = [];
Dat.ZoomStep = 0.2;
try
  Dat.ZoomLevel = getpref('Aedes','ZoomLevel');
  if Dat.ZoomLevel==0
    Dat.AxSize = 0;
    set(H.UITOGGLE_ZOOMNORM,'enable','on',...
                      'state','on')
  else
    %Dat.AxSize = Dat.RealAxSize*zoom_level;
  end
catch
  % Default axes sizes
  Dat.ZoomLevel=1;
  setpref('Aedes','ZoomLevel',1)
end

% Initialize axes ---------------------------------------------
tmp=get(H.SIDEBAR_FRAME,'position');
l_AxesPositions([],[],Dat.AxView,Dat.ZoomLevel);
set(H.IMAX1,'xlim',[0.5 Dat.ImageDim(Dat.DataInd,2)+0.5],...
                  'ylim',[0.5 Dat.ImageDim(Dat.DataInd,1)+0.49],...
                  'visible','off',...
                  'clim',Dat.Clim)
set(H.IMOVERLAYAX(1),'xlim',[0.5 Dat.ImageDim(Dat.DataInd,2)+0.5],...
                  'ylim',[0.5 Dat.ImageDim(Dat.DataInd,1)+0.49]);

set(H.IMAX2,'xlim',[0.5 Dat.ImageDim(Dat.DataInd,3)+0.5],...
            'ylim',[0.5 Dat.ImageDim(Dat.DataInd,1)+0.49],...
            'visible','off',...
            'clim',Dat.Clim)
set(H.IMOVERLAYAX(2),'xlim',[0.5 Dat.ImageDim(Dat.DataInd,3)+0.5],...
            'ylim',[0.5 Dat.ImageDim(Dat.DataInd,1)+0.49]);

set(H.IMAX3,'xlim',[0.5 Dat.ImageDim(Dat.DataInd,3)+0.5],...
            'ylim',[0.5 Dat.ImageDim(Dat.DataInd,2)+0.49],...
            'visible','off',...
            'clim',Dat.Clim)
set(H.IMOVERLAYAX(3),'xlim',[0.5 Dat.ImageDim(Dat.DataInd,3)+0.5],...
                  'ylim',[0.5 Dat.ImageDim(Dat.DataInd,2)+0.49])

% Get default View direction from preferences ------------------------------
try
  Dat.AxView = getpref('Aedes','AxView');
  if Dat.AxView==0
    set(H.UITOGGLE_VIEW3D,'state','on')
  elseif Dat.AxView==1
    set(H.UITOGGLE_VIEWX,'state','on')
  elseif Dat.AxView==2
    set(H.UITOGGLE_VIEWY,'state','on')
  else
    set(H.UITOGGLE_VIEWZ,'state','on')
  end
catch
  Dat.AxView = 0;
  setpref('Aedes','AxView',Dat.AxView)
  set(H.UITOGGLE_VIEW3D,'state','on')
end
if Dat.ImageDim(Dat.DataInd,3)==1
  Dat.AxView=3;
  set(H.UITOGGLE_VIEWZ,'state','on')
end

% Get grid info from preferences ---------------------------
try
  Dat.ShowGrid = getpref('Aedes','ShowGrid');
catch
  Dat.ShowGrid = false;
end
if Dat.ShowGrid
  set(H.UIVIEW_GRID,'checked','off')
  set(H.UITOGGLE_GRID,'state','off',...
                    'enable','on')
  l_ViewAxesUnits([],[],'pixel')
else
  set(H.UIVIEW_GRID,'checked','on')
  set(H.UITOGGLE_GRID,'state','on',...
                    'enable','on')
  l_ViewAxesUnits([],[],'pixel')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize toolbar buttons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize Zoom buttons ----------------------------
if Dat.ZoomLevel==0.5
  set(H.UIPUSH_ZOOMIN,'enable','on')
elseif Dat.ZoomLevel==10
  set(H.UIPUSH_ZOOMOUT,'enable','on')
else
  set([H.UIPUSH_ZOOMIN,H.UIPUSH_ZOOMOUT],'enable','on')
end
set(H.UITOGGLE_ZOOMNORM,'enable','on')
%Dat.ZoomLevel = zoom_level;

% Initialize mouse function buttons --------------------------
set([H.UITOGGLE_ARROW,H.UITOGGLE_PAN],'enable','on')

% Initialize ROI toolbar button -----------------------------
set(H.UIPUSH_NEWROI,'enable','on')

% Set colormap ---------------------------------------
try
  ColMap_str=getpref('Aedes','ColorMap');
catch
  ColMap_str='gray';
  Dat.ColMap=gray(256);
  setpref('Aedes','ColorMap',ColMap_str);
end
if strcmpi(ColMap_str,'Custom')
  try
	Dat.ColMap = getpref('Aedes','CustomColorMap');
	if ~isequal(size(Dat.ColMap),[256 3]) || ...
		~isnumeric(Dat.ColMap)
	  ColMap_str='gray';
	  Dat.ColMap=gray(256);
	  setpref('Aedes','ColorMap',ColMap_str);
	  h=errordlg(['Error initializing custom colormap!',...
		' Reverting to GRAY.'],...
		'Error!','modal');
	end
  catch
	ColMap_str='gray';
	Dat.ColMap=gray(256);
	setpref('Aedes','ColorMap',ColMap_str);
	h=errordlg(['Error initializing custom colormap!',...
	  ' Reverting to GRAY.'],...
	  'Error!','modal');
  end
  
end

% Set colorbar image -----------------------------
if Dat.ShowColorbar
  vis = 'on';
else
  vis = 'off';
end
H.COLORBAR_IM = image('parent',H.COLORBAR_AX,...
                      'cdatamapping','scaled',...
                      'cdata',gray(256),...
                      'visible',vis);
set(H.COLORBAR_AX,'xlim',[0.5 3.5],...
                  'ylim',[0.5 256+0.5],...
                  'ydir','normal',...
                  'layer','top',...
                  'xtick',[])
l_ChangeColorbarLimits([],[])

% Set default Window, Level, and Gamma.
% Store current valid value in the userdata...
Dat.Window = 100;
Dat.Level = 50;
set(H.WINDOW_EDIT,'string',num2str(Dat.Window),...
  'userdata',Dat.Window)
set(H.LEVEL_EDIT,'string',num2str(Dat.Level),...
  'userdata',Dat.Level)


% Set initial slices ---------------------------------
if ~Dat.isDataMixed && isfield(DATA{1},'DataFormat') && ...
	strncmpi(DATA{1}.DataFormat,'NIfTI',5)
  if DATA{1}.HDR.FileHeader.hist.qform_code>0
	q_b = DATA{1}.HDR.FileHeader.hist.quatern_b;
	q_c = DATA{1}.HDR.FileHeader.hist.quatern_c;
	q_d = DATA{1}.HDR.FileHeader.hist.quatern_d;
	q_a = sqrt(1.0-(q_b*q_b+q_c*q_c+q_d*q_d));
	R = [q_a*q_a+q_b*q_b-q_c*q_c-q_d*q_d,2*q_b*q_c-2*q_a*q_d,2*q_b*q_d+2*q_a*q_c;...
	  2*q_b*q_c+2*q_a*q_d,q_a*q_a+q_c*q_c-q_b*q_b-q_d*q_d,2*q_c*q_d-2*q_a*q_b; ...
	  2*q_b*q_d-2*q_a*q_c,2*q_c*q_d+2*q_a*q_b,q_a*q_a+q_d*q_d-q_c*q_c-q_b*q_b];
	qfac = DATA{1}.HDR.FileHeader.dime.pixdim(1);
	qoffset = [DATA{1}.HDR.FileHeader.hist.qoffset_x;...
	  DATA{1}.HDR.FileHeader.hist.qoffset_y;...
	  DATA{1}.HDR.FileHeader.hist.qoffset_z];
  end
else
  set(H.X_EDIT,'string',num2str(Dat.Slices(1)))
  set(H.Y_EDIT,'string',num2str(Dat.Slices(2)))
  set(H.Z_EDIT,'string',num2str(Dat.Slices(3)))
end

% Initialize volumes
if Dat.Vols>1
  set([H.V_BTN,H.V_EDIT],'enable','on')
  set(H.V_EDIT,'string','1');
  set(H.FLIP3D_V,'enable','on');
else
  set([H.V_BTN,H.V_EDIT],'enable','off')
  set(H.V_EDIT,'string','1');
  set(H.V_BTN,'fontweig','normal');
  set(H.FLIP3D_V,'enable','off');
end

% Set defaults for ROIs
Dat.RoiColors = [255 0 0;...
                0 255 0;...
                0 0 255;...
                255 255 0;...
                255 0 255;...
                0 255 255;...
                255 128 0;...
                128 255 0;...
                0 128 255;...
                128 0 255;...
                255 0 128;...
                0 255 128;...
                64 128 255;...
                128 64 255;...
                255 64 128;...
                255 128 64;...
                128 255 64;...
                64 255 128;...
                255 128 128;...
                128 255 128;...
                128 128 255;...
                255 255 128;...
                255 128 255;...
                128 255 255;...
                255 64 64;...
                64 255 64;...
                64 64 255;...
                255 255 64;...
                255 64 255;...
                64 255 255;...
                255 64 0;...
                64 255 0;...
                64 0 255;...
                0 64 255;...
                0 255 64;...
                255 0 64];
Dat.RoiView = [];
Dat.RoiUndo=[];
Dat.RoiSaved=1;
Dat.ShowRoiEdges = false;
Dat.CopiedROI = [];

% Set ROI transparency --------------------------
try
  Dat.RoiTransp=getpref('Aedes','RoiTransp');
catch
  Dat.RoiTransp=0.5;
end
set(H.ROI_TRANSP_SLIDER,'value',Dat.RoiTransp)
set(H.ROI_TRANSP_EDIT,'string',Dat.RoiTransp,...
                  'userdata',Dat.RoiTransp)

% Initialize Info Texts
try
  ShowInfoTexts = getpref('Aedes','ShowInfoTexts');
catch
  ShowInfoTexts = true;
end
if ShowInfoTexts
  set(H.UITOGGLE_INFO,'state','on')
  set([H.INFOTEXT,H.ROIINFOTEXT],...
      'visible','on')
else
  set(H.UITOGGLE_INFO,'state','off')
  set([H.INFOTEXT,H.ROIINFOTEXT],...
      'visible','off')
end

% Initialize Time Series view
Dat.ShowTimeSeries = false;

% Initialize Info Text Background
try
  ShowInfoBackGround = getpref('Aedes','ShowInfoTextBackGround');
catch
  ShowInfoBackGround = false;
end
if ShowInfoBackGround
  set([H.INFOTEXT,H.ROIINFOTEXT],...
      'backgroundcolor','w')
		if ~Dat.HG2graphics
			set([H.INFOTEXT,H.ROIINFOTEXT],...
				'erasemode','normal')
		end
  set(H.UIVIEW_INFOTXTBACKGROUND,'checked','on')
else
  set([H.INFOTEXT,H.ROIINFOTEXT],...
      'backgroundcolor','none')
		if ~Dat.HG2graphics
			set([H.INFOTEXT,H.ROIINFOTEXT],...
				'erasemode','normal')
		end
  set(H.UIVIEW_INFOTXTBACKGROUND,'checked','off')
end



% Initialize variables for imagemask
Dat.ImageOverlay = [];
Dat.ImageOverlayOrig = [];
Dat.isImageOverlayLoaded = false;

% Set figure resize function
set(H.FIG,'resizefcn',@l_ResizeFcn)

% Display data and initialize
l_DisplayData([],[],[1 2 3],1)

% Draw crassbars
l_ShowHideCrossbars;

% Set default zoom level
l_Zoom([],[],Dat.ZoomLevel)

% Set default colormap
l_ChangeColormap([],[],ColMap_str)

% Set default view direction
l_ChangeView([],[],Dat.AxView)

% Refresh Clim values
l_RefreshClim([],[]);

% Enable uicontrols
set(H.UICH_ENABLED,'enable','on')

% Disable "View VNMR Procpar" for non-VNMR data
%if ~isfield(DATA{1},'PROCPAR') || isempty(DATA{1}.PROCPAR)
%  set(H.viewprocpar_h,'enable','off')
%end

% If data is mixed type, disable rotate and image stack buttons
if ~Dat.isDataMixed
  set([H.UIPUSH_IMSTACK,H.UIPUSH_ROTATE],'enable','off')
end


if Dat.isDataMixed
  set(H.SL_SLIDER,'value',1,'min',0,'max',1,...
                  'sliderstep',[0.1 1],...
                  'enable','off')
  set([H.X_BTN,H.X_EDIT,H.Y_BTN,H.Y_EDIT,...
       H.Z_BTN,H.Z_EDIT,H.V_EDIT,H.V_BTN],'enable','off')
  set([H.X_BTN,H.Y_BTN,H.Z_BTN,H.V_BTN],'value',0,'fontweight','normal')
  set(H.Z_BTN,'value',1,'fontweight','bold')
  set([H.UIVIEW_3D,H.UIVIEW_Y,H.UIVIEW_X],'enable','off',...
                    'checked','off')
  set(H.UIVIEW_Z,'enable','off','checked','on')
  set([H.UITOGGLE_VIEW3D,H.UITOGGLE_VIEWY,H.UITOGGLE_VIEWX],'enable', ...
                    'off')
  set(H.UIMENU_ROTATEFLIP,'enable','on')
  set(H.UIMENU_ROTATEFLIP3D,'enable','off')
  set(H.UIEDIT_IMSTACK,'enable','on')
  Dat.AxView=3;
  
  % Setup slice slider for Mixed data
  if Dat.DataL>1
    set([H.Z_BTN,H.Z_EDIT],'enable','on')
    set(H.SL_SLIDER,'value',1,'min',1,'max',Dat.DataL,...
                    'sliderstep',[1/(Dat.DataL-1) ...
                        ceil(Dat.DataL/10)/(Dat.DataL-1)],...
                    'enable','on')
  end

	set(H.CLIMRANGE_POPUP,'enable','on','value',3)
else
  if Dat.ImageDim(Dat.DataInd,3)==1
    set(H.SL_SLIDER,'value',1,'min',1,...
      'max',Dat.ImageDim(Dat.DataInd,4),...
      'sliderstep',[1/max(1,Dat.ImageDim(Dat.DataInd,4)-1) ...
      ceil(Dat.ImageDim(Dat.DataInd,4)/10)/max(1,Dat.ImageDim(Dat.DataInd,4)-1)],...
      'enable','on')
    set([H.X_BTN,H.X_EDIT,H.Y_BTN,H.Y_EDIT,...
      H.Z_BTN,H.Z_EDIT],'enable','on')
    set([H.Z_BTN,H.Z_EDIT],'enable','off')
    set([H.X_BTN,H.Y_BTN,H.Z_BTN],'value',0,'fontweight','normal')
    set(H.V_BTN,'value',1,'fontweight','bold')
  else
    set(H.SL_SLIDER,'value',1,'min',1,...
      'max',Dat.ImageDim(Dat.DataInd,3),...
      'sliderstep',[1/(Dat.ImageDim(Dat.DataInd,3)-1) ...
      ceil(Dat.ImageDim(Dat.DataInd,3)/10)/(Dat.ImageDim(Dat.DataInd,3)-1)])
    set([H.X_BTN,H.X_EDIT,H.Y_BTN,H.Y_EDIT,...
      H.Z_BTN,H.Z_EDIT],'enable','on')
    set([H.X_BTN,H.Y_BTN,H.Z_BTN],'value',0,'fontweight','normal')
    set(H.Z_BTN,'value',1,'fontweight','bold')
  end
  set([H.UIVIEW_3D,H.UIVIEW_X,H.UIVIEW_Y,H.UIVIEW_Z],'enable','on')
  set([H.UITOGGLE_VIEW3D,H.UITOGGLE_VIEWX,H.UITOGGLE_VIEWY,H.UITOGGLE_VIEWZ],'enable','on')
  set(H.UIMENU_ROTATEFLIP3D,'enable','on')
  set(H.UIMENU_ROTATEFLIP,'enable','off')
	set(H.CLIMRANGE_POPUP,'value',3,'enable','off')
end


% Reset window/level sliders
set(H.WINDOW_SLIDER,'value',100)
set(H.LEVEL_SLIDER,'value',50)

% Set file name to figure window title
set(H.FIG,'Name',['Aedes 1.0 - (' DATA{Dat.DataInd}.HDR.fpath,DATA{Dat.DataInd}.HDR.fname ...
                  ')'])

%% Remove blank frame
set(H.BLANK_FRAME,'visible','off')

%% Run resize
l_ResizeFcn([],[])

% A quick hack option for disabling mouse wheel zooming, which apparently
% has issues in OSX with magic mouse...
try
if ispref('Aedes','DisableMouseWheel')
  Dat.DisableMouseWheel = getpref('Aedes','DisableMouseWheel');
else
        Dat.DisableMouseWheel = false;
        setpref('Aedes','DisableMouseWheel',false);
    end
catch
  Dat.DisableMouseWheel = false;
end

%% Set callback for Mouse scroll wheel
if Dat.MatlabVersion>=7.05 && ~Dat.DisableMouseWheel
  set(H.FIG,'WindowScrollWheelFcn',@l_MouseWheelFcn)
end



catch
  aedes_errordump(lasterror);
end
end % function l_Initialize(h,



%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Display data
%%
%%%%%%%%%%%%%%%%%%%%%%%%
function l_DisplayData(h,evd,update_ind,init)
try
if nargin<3
  update_ind=[1 2 3];
  init=0;
elseif nargin==3
  init=0;
end

if init
  update_ind=[1 2 3];
elseif Dat.isDataMixed
  update_ind = 3;
end

% Extract slice data
updateX = any(update_ind==1);
updateY = any(update_ind==2);
updateZ = any(update_ind==3);

% Get current slice for Z direction -------------
if updateZ
  slice1 = DATA{Dat.DataInd}.FTDATA(:,:,Dat.Slices(3),Dat.CurrentVol);
end

% Get current slice for Y direction -------------
if updateY
  slice2 = reshape(DATA{Dat.DataInd}.FTDATA(:,Dat.Slices(2),:,Dat.CurrentVol),...
                   Dat.ImageDim(Dat.DataInd,1),Dat.ImageDim(Dat.DataInd,3));
end

% Get current slice for Z direction -------------
if updateX
  slice3 = reshape(DATA{Dat.DataInd}.FTDATA(Dat.Slices(1),:,:,Dat.CurrentVol),...
                   Dat.ImageDim(Dat.DataInd,2),Dat.ImageDim(Dat.DataInd,3));
end


% Update slice editboxes
if Dat.isDataMixed
  set(H.X_EDIT,'string',sprintf('%d',Dat.Slices(1)))
  set(H.Y_EDIT,'string',sprintf('%d',Dat.Slices(2)))
  set(H.Z_EDIT,'string',sprintf('%d',Dat.DataInd))
else
  set(H.X_EDIT,'string',sprintf('%d',Dat.Slices(1)))
  set(H.Y_EDIT,'string',sprintf('%d',Dat.Slices(2)))
  set(H.Z_EDIT,'string',sprintf('%d',Dat.Slices(3)))
end



% Display slices in different orientations
%if ~isfield(H,'IM')
if init
  H.IM(1) = image('parent',H.IMAX1,...
                  'CData',slice1,...
                  'cdatamapping','scaled',...
                  'hittest','on',...
                  'buttondownfcn',@l_SetMouseGestures);
  
  H.IM(2) = image('parent',H.IMAX2,...
                  'CData',slice2,...
                  'cdatamapping','scaled',...
                  'hittest','on',...
                  'buttondownfcn',@l_SetMouseGestures);

  
  H.IM(3) = image('parent',H.IMAX3,...
                  'cdata',slice3,...
                  'cdatamapping','scaled',...
                  'hittest','on',...
                  'buttondownfcn',@l_SetMouseGestures);

else
  if updateZ
    set(H.IM(1),'CData',slice1)
  end
  if updateY
    set(H.IM(2),'CData',slice2)
  end
  if updateX
    set(H.IM(3),'CData',slice3)
  end
  % Show ROI
  l_RefreshRoi(update_ind)
end

% Change figure window title if necessary
if Dat.isDataMixed
  set(H.FIG,'Name',['Aedes 1.0 - (' ...
                    DATA{Dat.DataInd}.HDR.fpath,...
                    DATA{Dat.DataInd}.HDR.fname ')'])
end

% Update Info Text
l_UpdateInfoText([],[])

% Update ROI Info Text
if Dat.isDataMixed
  l_UpdateRoiInfoText([],[])
  l_AxesPositions([],[],3,Dat.ZoomLevel)
end

% Draw ROI edges
if Dat.ShowRoiEdges
  l_ShowRoiEdges([],[],'all')
end

% Update Image Mask
if Dat.isImageOverlayLoaded
  l_UpdateImageOverlay([],[],update_ind)
end


catch
  aedes_errordump(lasterror);
end
end % function l_DisplayData(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Edit Image Stack
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_EditImageStack(h,evd)
try

% Call stack editor GUI
[filelist,rem_ind,add_ind,sort_ind] = aedes_editstack(DATA);
if isempty(filelist)
  return
end

% Look for a quick return
if ~any(add_ind) && ~any(rem_ind) && ...
      length(filelist)==length(DATA) && all(sort_ind==(1:length(sort_ind)))
  % Nothing to do...
  return
end

% Construct a new data structure
DATA_tmp = cell(1,length(filelist));
old_sort = sort_ind;
old_sort(old_sort==0)=[];
DataRotation = zeros(1,length(filelist));
DataFlip = zeros(1,length(filelist));
SliceClim = zeros(length(filelist),2);

% Assign data from the old structure to the new one in the right order
DATA_tmp(sort_ind~=0) = {DATA{old_sort}};
DataRotation(sort_ind~=0) = Dat.DataRotation(old_sort);
DataFlip(sort_ind~=0) = Dat.DataFlip(old_sort);
SliceClim(sort_ind~=0,:) = Dat.SliceClim(old_sort,:);

% Read new data
skippedFiles={};
if any(add_ind)
  [hCalcWait,txt_h] = aedes_calc_wait('Reading new data files...'); 
  ind=find(add_ind);
  for ii=ind
    data=aedes_data_read(filelist{ii});
    if isempty(data)
      fprintf(1,'Could not read data file "%s". Skipping file...\n',filelist{ii})
      skippedFiles{end+1}=filelist{ii};
    elseif length(size(data.FTDATA))>2
      fprintf(1,'Data file "%s" contains multiple slices. Skipping file...\n',filelist{ii})
      skippedFiles{end+1}=filelist{ii};
    else
      DATA_tmp{ii} = data;
    end
  end
  
  % Remove possible empty fields
  empty_ind = cellfun('isempty',DATA_tmp);
  DATA_tmp(empty_ind)=[];
  DataRotation(empty_ind)=[];
  DataFlip(empty_ind)=[];
  SliceClim(empty_ind,:)=[];
  
  % Remove aedes_calc_wait
  delete(hCalcWait)
end

% Check that the new DATA structure is not empty
if all(cellfun('isempty',DATA_tmp))
  h=warndlg({'The stack editing has failed!','',...
             'Could not read new files to the stack!'},...
            'Stack editing failed','modal');
  return
end


% Update ROIs if they exist
ROI_tmp=[];
if ~isempty(ROI)
  
  % Get file and path names
  for ii=1:length(DATA_tmp)
    new_fnames{ii} = DATA_tmp{ii}.HDR.fname;
    new_fpaths{ii} = DATA_tmp{ii}.HDR.fpath;
    new_fullfiles{ii} = [DATA_tmp{ii}.HDR.fpath,...
                        DATA_tmp{ii}.HDR.fname];
  end
  
  % Construct a new ROI structure
  ROI_tmp = ROI;
  
  for ii=1:length(ROI)
    ROI_tmp(ii).voxels = cell(1,length(filelist));
    ROI_tmp(ii).fpath = cell(1,length(DATA_tmp));
    ROI_tmp(ii).fname = cell(1,length(DATA_tmp));
    
    ROI_tmp(ii).voxels(sort_ind~=0) = {ROI(ii).voxels{old_sort}};
    ROI_tmp(ii).fpath = new_fpaths;
    ROI_tmp(ii).fname = new_fpaths;
    
    % Check added files
    for kk=find(add_ind)
      ROI_tmp(ii).voxels{kk} = false(size(DATA_tmp{kk}.FTDATA));
    end
    
    % Delete for empty cell values (files that could not be readed)
    empty_ind=find(cellfun('isempty',ROI_tmp(ii).voxels));
    ROI_tmp(ii).voxels(empty_ind)=[];
  end
end

Clim = Dat.Clim;

% Close file
l_CloseFile([],[])

% Update structures
DATA = DATA_tmp;
ROI = ROI_tmp;


% Re-initialize GUI
l_Initialize([],[])

% Update rotation/flip
Dat.DataRotation = DataRotation;
Dat.DataFlip = DataFlip;

% Replace old slice Clim values
clim_ind = ~all(SliceClim==0,2);
Dat.SliceClim(clim_ind,:)=SliceClim(clim_ind,:);

% Refresh data
l_DisplayData([],[])

% Refresh window/level
l_SetWindowLevel(0,[],Clim,0)

% Load ROI(s)
if ~isempty(ROI)
  l_RoiLoad([],[],'')
end

% Warn about skipped files
if ~isempty(skippedFiles)
  h=warndlg({'The following file(s) could not be read and were skipped:','',skippedFiles{:}},...
            'Some files could not be read','modal');
end

catch
  aedes_errordump(lasterror);
end
end % function l_EditImageStack(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Update Image Mask
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_UpdateImageOverlay(h,evd,update_ind)
try
  
% Return immediately if mask is not currently loaded
if ~Dat.isImageOverlayLoaded || Dat.isImageOverlayHidden
  return
end

% If the data underlying the 3D overlay data is 4D, use only the first
% volume of the overlay...
if size(Dat.ImageOverlay,4)==1
  volInd = 1;
else
  volInd = Dat.CurrentVol;
end

% Extract slice data
updateX = any(update_ind==1);
updateY = any(update_ind==2);
updateZ = any(update_ind==3);

if ~Dat.isOverlayRGB
  cmap = Dat.ImageOverlayCmap;
  clim = round(((Dat.ImageOverlayClim-Dat.ImOverlayMin)*256)./...
	(Dat.ImOverlayMax-Dat.ImOverlayMin));
  if clim(1)==clim(2)
	if clim(1)==0
	  clim(2)=1;
	elseif clim(1)==256
	  clim(1)=255;
	end
  end
  thold = round(((Dat.ImageOverlayThold-Dat.ImOverlayMin)*256)/...
	(Dat.ImOverlayMax-Dat.ImOverlayMin));
  
end
alpha_val = Dat.ImageOverlayAlpha;

% Get current slice for X direction -------------
if updateZ
  if ~Dat.isOverlayRGB
    % Convert indexed image to RGB image
    slice1_ind=Dat.ImageOverlay(:,:,Dat.Slices(3),volInd);
    slice1_ind=double(slice1_ind);
	
	% Get thresholded alpha indices
	if Dat.ImageOverlayTholdDirPos==1
	  slice1_alpha_th = slice1_ind<thold;
	else
	  slice1_alpha_th = slice1_ind>thold;
	end
	
	% Get clim alpha indices
	%slice1_alpha_clim = ( slice1_ind>=clim(1) & slice1_ind<=clim(2) );
	
    slice1_ind(slice1_ind<clim(1))=clim(1);
    slice1_ind(slice1_ind>clim(2))=clim(2);
	
    slice1_ind=ceil((slice1_ind-clim(1))./diff(clim)*255+1);
    %disp(['1: ',num2str([max(slice1_ind(:)) min(slice1_ind(:))])])
	
    sz = size(slice1_ind);
    slice1 = zeros([sz(1) sz(2) 3],'single');
    slice1(:,:,1)= reshape(cmap(slice1_ind,1),sz);
    slice1(:,:,2)= reshape(cmap(slice1_ind,2),sz);
    slice1(:,:,3)= reshape(cmap(slice1_ind,3),sz);
  else
    %slice1 = zeros([sz(1) sz(2) 3],'single');
    slice1 = squeeze(Dat.ImageOverlay(:,:,Dat.Slices(3),:));
  end
end

% Get current slice for Y direction -------------
if updateY
  if ~Dat.isOverlayRGB
    slice2_ind=reshape(Dat.ImageOverlay(:,Dat.Slices(2),:,volInd),...
	  Dat.ImageDim(Dat.DataInd,1), ...
	  Dat.ImageDim(Dat.DataInd,3));
	slice2_ind=double(slice2_ind);
	
	% Get thresholded alpha indices
	if Dat.ImageOverlayTholdDirPos==1
	  slice2_alpha_th = slice2_ind<thold;
	else
	  slice2_alpha_th = slice2_ind>thold;
	end
	
	% Get clim alpha indices
	%slice2_alpha_clim = ( slice2_ind>=clim(1) & slice2_ind<=clim(2) );
	
	
    slice2_ind(slice2_ind<clim(1))=clim(1);
    slice2_ind(slice2_ind>clim(2))=clim(2);
   
    slice2_ind=ceil((slice2_ind-clim(1))./diff(clim)*255+1);
    %disp(['2: ',num2str([max(slice2_ind(:)) min(slice2_ind(:))])])
    sz = size(slice2_ind);
    slice2 = zeros([sz(1) sz(2) 3],'single');
    slice2(:,:,1)= reshape(cmap(slice2_ind,1),sz);
    slice2(:,:,2)= reshape(cmap(slice2_ind,2),sz);
    slice2(:,:,3)= reshape(cmap(slice2_ind,3),sz);
  else
    %slice2 = zeros([sz(1) sz(2) 3],'single');
    slice2 = squeeze(Dat.ImageOverlay(:,Dat.Slices(2),:,:));
  end
end

% Get current slice for Z direction -------------
if updateX
  if ~Dat.isOverlayRGB
    slice3_ind=reshape(Dat.ImageOverlay(Dat.Slices(1),:,:,volInd),...
	  Dat.ImageDim(Dat.DataInd,2),Dat.ImageDim(Dat.DataInd,3));
    slice3_ind=double(slice3_ind);
	
	% Get thresholded alpha indices
	if Dat.ImageOverlayTholdDirPos==1
	  slice3_alpha_th = slice3_ind<thold;
	else
	  slice3_alpha_th = slice3_ind>thold;
	end
	
	% Get clim alpha indices
	%slice3_alpha_clim = ( slice3_ind>=clim(1) & slice3_ind<=clim(2) );
	
    slice3_ind(slice3_ind<clim(1))=clim(1);
    slice3_ind(slice3_ind>clim(2))=clim(2);
    
    slice3_ind=ceil((slice3_ind-clim(1))./diff(clim)*255+1);
    %disp(['3: ',num2str([max(slice3_ind(:)) min(slice3_ind(:))])])
    sz = size(slice3_ind);
    slice3 = zeros([sz(1) sz(2) 3],'single');
    slice3(:,:,1)= reshape(cmap(slice3_ind,1),sz);
    slice3(:,:,2)= reshape(cmap(slice3_ind,2),sz);
    slice3(:,:,3)= reshape(cmap(slice3_ind,3),sz);
  else
    %slice3 = zeros([sz(1) sz(2) 3],'single');
    slice3 = squeeze(Dat.ImageOverlay(Dat.Slices(1),:,:,:));
  end
end


% Update image mask
%if isempty(H.IMOVERLAY)
  
  % Create alpha data
  if ~Dat.isOverlayRGB
	if updateZ
 	  slice1_alpha = ones(size(slice1_ind))*alpha_val;
      %slice1_alpha(slice1_alpha_clim) = alpha_val;
	  slice1_alpha(slice1_alpha_th) = 0;
	end
    
    if updateY
	  slice2_alpha = ones(size(slice2_ind))*alpha_val;
      %slice2_alpha(slice2_alpha_clim) = alpha_val;
	  slice2_alpha(slice2_alpha_th) = 0;
    end
    
    if updateX
	  slice3_alpha = ones(size(slice3_ind))*alpha_val;
      %slice3_alpha(slice3_alpha_clim) = alpha_val;
	  slice3_alpha(slice3_alpha_th) = 0;
    end
  else
	if updateZ
	  slice1_alpha = zeros([size(slice1,1),size(slice1,2)]);
	  slice1_alpha(:,:) = alpha_val;
	end
    
    if updateY
      slice2_alpha = zeros([size(slice2,1),size(slice2,2)]);
      slice2_alpha(:,:) = alpha_val;
    end
    
    if updateX
      slice3_alpha = zeros([size(slice3,1),size(slice3,2)]);
      slice3_alpha(:,:) = alpha_val;
    end
  end

  % Update image mask
  if isempty(H.IMOVERLAY)
	if updateZ
	  H.IMOVERLAY(1) = image('parent',H.IMOVERLAYAX(1),...
		'CData',slice1,...
		'cdatamapping','scaled',...
		'hittest','off',...
		'AlphaDataMapping','none',...
		'AlphaData',slice1_alpha,...
		'visible','on');
	end
	if updateY
	  H.IMOVERLAY(2) = image('parent',H.IMOVERLAYAX(2),...
		'CData',slice2,...
		'cdatamapping','scaled',...
		'hittest','off',...
		'AlphaDataMapping','none',...
		'AlphaData',slice2_alpha,...
		'visible','on');
	end
	
	if updateX
	  H.IMOVERLAY(3) = image('parent',H.IMOVERLAYAX(3),...
		'cdata',slice3,...
		'cdatamapping','scaled',...
		'hittest','off',...
		'AlphaDataMapping','none',...
		'AlphaData',slice3_alpha,...
		'visible','on');
	end
	
	if Dat.HG2graphics
    ind = find([updateZ,updateY,updateX]);
		for ii=ind
			set(H.IMOVERLAY(ii),'hittest','on',...
				'buttondownfcn',@l_SetMouseGestures)
		end
		try
			delete(H.CROSSBAR_LN)
			H.CROSSBAR_LN=[];
		end
		l_ShowHideCrossbars([],[]);
	end
	
  else
	if updateZ
	  set(H.IMOVERLAY(1),'cdata',slice1,...
		'AlphaDataMapping','none',...
		'AlphaData',slice1_alpha,...
		'visible','on')
	end
	
	if updateY
	  set(H.IMOVERLAY(2),'cdata',slice2,...
		'AlphaDataMapping','none',...
		'AlphaData',slice2_alpha,...
		'visible','on')
	end
	
	if updateX
	  set(H.IMOVERLAY(3),'cdata',slice3,...
		'AlphaDataMapping','none',...
		'AlphaData',slice3_alpha,...
		'visible','on')
	end
	
  end

catch
  aedes_errordump(lasterror);
end
end % function l_UpdateImageOverlay(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Load Image Overlay Data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_LoadImageOverlay(h,evd,opt)
%try
if Dat.isDataMixed && length(DATA)>1
  h=warndlg('Image overlay is not supported image stack data',...
            'Image overlay not supported','modal');
  return
end


% Load file or variable
if strcmpi(opt,'file')
  
  % Get default dir
  try
    default_dir = getpref('Aedes','GetOverlayFileDir');
  catch
    default_dir = [pwd,filesep];
  end
  
  % Read data
  try
    tmp_data=aedes_data_read([],'default_dir',default_dir);
  catch
    h=errordlg({'Error loading image overlay!','',lasterr},...
               'Error loading image overlay','modal');
    return
  end
  if isempty(tmp_data)
	% Action cancelled
	return
  end
  
  % Save the overlay directory to preferences
  if isfield(tmp_data,'HDR')
    if iscell(tmp_data.HDR.fpath)
      setpref('Aedes','GetOverlayFileDir',tmp_data.HDR.fpath{1})
    else
      setpref('Aedes','GetOverlayFileDir',tmp_data.HDR.fpath)
    end
  end
  
elseif strcmpi(opt,'var') % Load from variable
  
    % Get workspace variables
  str = evalin('base','who');
  if isempty(str)
    h=errordlg({'Error loading image overlay!','',...
               'Matlab workspace doesn''t contain any variables.'},...
               'Error loading image overlay','modal');
    return
  end
  
  % Show load dialog
  [ind,isok]=listdlg('PromptString',...
                     {'Select a workspace variable to','load as image overlay'},...
                     'SelectionMode','single',...
                     'ListString',str);
  if isok==0
    % Action cancelled
    return
  end
  
  % Load variable
  tmp=evalin('base',str{ind});
  if isnumeric(tmp) || islogical(tmp)
    tmp_data.FTDATA = tmp;
    clear tmp;
  elseif isstruct(tmp) && isfield(tmp,'FTDATA')
    tmp_data = tmp;
    clear tmp;
  else
    clear tmp;
    h=errordlg({'Error loading image overlay!','',...
                'Loaded variable invalid.'},...
               'Error loading image overlay','modal');
    return
  end
  
end

% If overlay was given as input argument to Aedes
if strcmpi(opt,'inputarg') && ~isempty(inputOverlay)
  if isstruct(inputOverlay) && isfield(inputOverlay,'FTDATA')
    tmp_data.FTDATA = inputOverlay.FTDATA;
  else
    tmp_data.FTDATA = inputOverlay;
  end
end

% Check image size
sz1 = size(tmp_data.FTDATA);
tmp=ones(1,4);
tmp(:,1:length(sz1))=sz1;
sz1=tmp;
sz2 = size(DATA{1}.FTDATA);
tmp=ones(1,4);
tmp(:,1:length(sz2))=sz2;
sz2=tmp;
dim1 = ndims(tmp_data.FTDATA);
dim2 = ndims(DATA{1}.FTDATA);



% Check if overlay data is RGB data
if dim1==4 && sz1(4)==3
  Dat.isOverlayRGB = true;
else
  Dat.isOverlayRGB = false;
end

% If only one slice is loaded and there is size mismatch, ask if the user
% wants resize the overlay...
if sz1(1)==sz1(2) && sz2(1)==sz2(2) && ...
	~ismember(sz1,sz2,'rows') && sz2(3)==1 && sz2(4)==1
  resp=questdlg({'Data/overlay size mismatch!',...
	'',['Data size: ',num2str(sz2(1)),'x',num2str(sz2(2)),'x',num2str(sz2(3))],...
	['Overlay size: ',num2str(sz1(1)),'x',num2str(sz1(2)),'x',num2str(sz1(3))],...
	'','','Do you want to resize the overlay to match the data?'},...
	'Resize overlay to data size?','Resize Overlay','Abort','Abort');
  if isempty(resp) || strcmpi(resp,'Abort')
	clear tmp_data;
	return
  end
  tmp_data.FTDATA=imresize(tmp_data.FTDATA,[sz2(1) sz2(2)]);
else
  if ~Dat.isOverlayRGB && ~ismember(sz1,sz2,'rows')
    % Allow Overlay loading if only 4th dimension differs
    if ~ismember(sz1(1:3),sz2(1:3),'rows')
      clear tmp_data;
      h=errordlg({'Data/overlay size mismatch!',...
        '',['Data size: ',num2str(sz2(1)),'x',num2str(sz2(2)),'x',num2str(sz2(3))],...
        ['Overlay size: ',num2str(sz1(1)),'x',num2str(sz1(2)),'x',num2str(sz1(3))]},...
        'Error loading overlay!','modal');
      return
    end
  elseif Dat.isOverlayRGB && ~ismember(sz1(1:3),sz2(1:3),'rows')
    clear tmp_data;
    h=errordlg({'Data/overlay size mismatch!',...
      '',['Data size: ',num2str(sz2(1)),'x',num2str(sz2(2)),'x',num2str(sz2(3))],...
      ['Overlay size: ',num2str(sz1(1)),'x',num2str(sz1(2)),'x',num2str(sz1(3))]},...
      'Error loading overlay!','modal');
    return
  end
end

Dat.ImageOverlay=tmp_data.FTDATA;
Dat.ImageOverlayOrig = tmp_data.FTDATA;
clear tmp_data;

Dat.ImOverlayMax = double(max(max(max(max(Dat.ImageOverlay)))));
Dat.ImOverlayMin = double(min(min(min(min(Dat.ImageOverlay)))));
if ~Dat.isOverlayRGB
  if ~any(strcmpi(class(Dat.ImageOverlay),{'single','double'}))
	Dat.ImageOverlay = single(Dat.ImageOverlay);
  end
  
  if any(strcmpi(class(Dat.ImageOverlay),{'single','double'}))
    % Scale data to range 0..255
    Dat.ImageOverlay = Dat.ImageOverlay-min(min(min(min(Dat.ImageOverlay))));
    Dat.ImageOverlay = Dat.ImageOverlay./ ...
      (max(max(max(max(Dat.ImageOverlay)))));
    Dat.ImageOverlay = Dat.ImageOverlay*255;
    Dat.ImageOverlay = uint8(round(Dat.ImageOverlay));
  else
    % Make sure that the minimum of data is 1
    Dat.ImageOverlay = max(1,Dat.ImageOverlay);
    
    % Generate indices in uint8 range
    ind=round(linspace(1,256,Dat.ImOverlayMax));
    ind = uint8(ind);
    
    % Scale overlay data to the uint8 range (0:255)
    Dat.ImageOverlay(:)=ind(Dat.ImageOverlay);
    Dat.ImageOverlay = uint8(Dat.ImageOverlay);
	end
end

Dat.isImageOverlayLoaded = true;
Dat.isImageOverlayHidden = false;

% Set defaults for overlay from preferences (if available)
if ~Dat.isOverlayRGB
  Dat.ImageOverlayThold=Dat.ImOverlayMin;
  Dat.ImageOverlayClim=[Dat.ImOverlayMin Dat.ImOverlayMax];
else
  Dat.ImageOverlayThold = 0;
  Dat.ImageOverlayClim = [0 1];
end

if ispref('Aedes','ImageOverlayAlpha')
  Dat.ImageOverlayAlpha=getpref('Aedes','ImageOverlayAlpha');
else
  Dat.ImageOverlayAlpha=0.6;
end

if ispref('Aedes','ImageOverlayTholdDirPos')
  Dat.ImageOverlayTholdDirPos = getpref('Aedes','ImageOverlayTholdDirPos');
else
  Dat.ImageOverlayTholdDirPos = 1;
end

if ispref('Aedes','ImageOverlayCmap')
  Dat.ImageOverlayCmapStr = getpref('Aedes','ImageOverlayCmap');
  Dat.ImageOverlayCmap = feval(Dat.ImageOverlayCmapStr,256);
else
  Dat.ImageOverlayCmap = jet(256);
  Dat.ImageOverlayCmapStr = 'jet';
end

% If options to the overlay were given as input argument...
if strcmpi(opt,'inputarg') && ~isempty(inputOverlay)
  if isstruct(inputOverlay)
    if isfield(inputOverlay,'ImageOverlayCmapStr')
      Dat.ImageOverlayCmapStr = inputOverlay.ImageOverlayCmapStr;
      Dat.ImageOverlayCmap = feval(Dat.ImageOverlayCmapStr,256);
    end
    if isfield(inputOverlay,'ImageOverlayTholdDirPos')
      Dat.ImageOverlayTholdDirPos = inputOverlay.ImageOverlayTholdDirPos;
    end
    if isfield(inputOverlay,'ImageOverlayAlpha')
      Dat.ImageOverlayAlpha = inputOverlay.ImageOverlayAlpha;
    end
    if isfield(inputOverlay,'ImageOverlayThold')
      Dat.ImageOverlayThold = inputOverlay.ImageOverlayThold;
    end
    if isfield(inputOverlay,'ImOverlayMin')
      Dat.ImageOverlayClim(1) = inputOverlay.ImOverlayMin;
    end
    if isfield(inputOverlay,'ImOverlayMax')
      Dat.ImageOverlayClim(2) = inputOverlay.ImOverlayMax;
    end
    if isfield(inputOverlay,'fMRIonsets') && ...
        isfield(inputOverlay,'fMRIdurats')
      Dat.fMRIonsets = inputOverlay.fMRIonsets;
      Dat.fMRIdurats = inputOverlay.fMRIdurats;
    end
  end
  inputOverlay = [];
end

% Enable uimenus
set([H.UIOVERLAY_CONTROLS,...
	H.UIOVERLAY_SAVE,...
	H.UIOVERLAY_DELETE],'enable','on')


% Refresh screen
if Dat.isDataMixed
  l_DisplayData([],[],1)
else
  l_DisplayData([],[])
end

% Draw crossbars
l_ShowHideCrossbars([],[])

% Show Overlay controls
l_OverlayControls([],[],'show');

%catch
%  aedes_errordump(lasterror);
%end
end % function l_LoadImageOverlay(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Image Overlay Controls
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_OverlayControls(h,evd,opt)
%try
if strcmpi(opt,'show')
  
  %% If overlay control window is already opened...
  if isfield(H,'OVERLAY_CONTROL_FIG') && ~isempty(H.OVERLAY_CONTROL_FIG) ...
      && ishandle(H.OVERLAY_CONTROL_FIG)
    figure(H.OVERLAY_CONTROL_FIG)
    return
  end
  
  %% Load default font and colors
  FigColor=get(0,'DefaultUicontrolBackgroundcolor');
  GD=aedes_gui_defaults;
  
  % Overlay controls figure
  fig_h = 250;
  fig_w = 285;
	fig_location = aedes_dialoglocation([fig_w,fig_h]);
	fig_left = fig_location(1);
	fig_bottom = fig_location(2);
  
  % Use previous position for the overlay controls
  if ispref('Aedes','OverlayControlPos')
		tmp_pos=getpref('Aedes','OverlayControlPos');
		
		% Check that the figure is on the screen
		scrsz = get(0,'Screensize');
		if scrsz(3)>tmp_pos(1) && scrsz(4)>tmp_pos(2)
			fig_left = tmp_pos(1);
			fig_bottom = tmp_pos(2);
		end
  end
  
  H.OVERLAY_CONTROL_FIG = figure('units','pixel',...
    'position',...
    [fig_left fig_bottom ...
    fig_w fig_h],...
    'Name','Image Overlay Controls',...
    'numbertitle','off',...
    'tag','aedes_overlay_controls_fig',...
    'Toolbar','none',...
    'Color',FigColor,...
    'Menubar','none',...
    'DockControls','off',...
    'renderer','painters',...
    'resize','off',...
    'CloseRequestFcn',{@l_OverlayControlCB,'close_window'},...
    'Handlevisibility','off');
  if ~GD.HG2graphics
		set(H.OVERLAY_CONTROL_FIG, 'DoubleBuffer','on')
	end
  % Suppress warning from get(fh,'javaFrame') generated in Matlab R2008a->
  if Dat.MatlabVersion>=7.06
    warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
  end
  
  % Try to set the figure as floating i.e. "always on top"
  jf = get(H.OVERLAY_CONTROL_FIG,'javaFrame');
	drawnow
	drawnow
  pause(0.2) % Make Matlab more stable...
	drawnow
	drawnow
	if ~GD.HG2graphics
		if Dat.MatlabVersion>=7.13
			wh = jf.fHG1Client.getWindow;
		else
			wh = jf.fFigureClient.getWindow;
		end
	else
		wh = handle(jf.fHG2Client.getWindow);
	end
	drawnow
	drawnow
  pause(0.1)
  wh.setAlwaysOnTop(true);
	drawnow
	drawnow
  pause(0.1)
							   
  uipanel_h = uipanel('parent',H.OVERLAY_CONTROL_FIG,...
                      'units','pixel',...
                      'position',[10 10 fig_w-20 fig_h-20]);
  
  % Threshold controls
  th_tx = uicontrol('parent',uipanel_h,...
                    'units','pixel',...
                    'position',[10 200 150 15],...
                    'style','text',...
                    'string','Threshold',...
                    'horizontalalign','left');
  tmp = get(th_tx,'position');
  if Dat.isOverlayRGB
	sl_value = 0;
  else
	sl_value = round(((Dat.ImageOverlayThold-Dat.ImOverlayMin)*256)/...
	  (Dat.ImOverlayMax-Dat.ImOverlayMin));
  end
  th_slider = uicontrol('parent',uipanel_h,...
                        'units','pixel',...
                        'position',[tmp(1) tmp(2)-20 190 20],...
                        'style','slider',...
                        'Min',0,...
                        'Max',255,...
                        'SliderStep',[1/255 10/255],...
                        'value',sl_value,...
                        'backgroundcolor',GD.col.slider);
	if ~GD.HG2graphics
		thSliderListener = handle.listener(th_slider,...
			'ActionEvent',...
			{@l_OverlayControlCB,'th_slider'});
		setappdata(th_slider,'ThresholdSliderListener',...
			thSliderListener)
	else
		thSliderListener = addlistener(th_slider,...
			'ContinuousValueChange',...
			@l_OverlayControlCB);
		setappdata(th_slider,'ThresholdSliderListener',...
			thSliderListener)
	end
  
  tmp = get(th_slider,'position');
  th_edit = uicontrol('parent',uipanel_h,...
                      'units','pixel',...
                      'position',[tmp(1)+tmp(3)+5 tmp(2) 55 20],...
                      'style','edit',...
                      'string',num2str(Dat.ImageOverlayThold),...
                      'backgroundcolor','w',...
                      'callback',{@l_OverlayControlCB,'th_edit'});
  
  % Threshold direction checkbox					
  tmp = get(th_slider,'position');
  th_direction = uicontrol('parent',uipanel_h,...
	'units','pixel',...
	'position',[tmp(1) tmp(2)-tmp(4)-5 240 20],...
	'style','checkbox',...
	'string','Show values over the threshold',...
	'value',Dat.ImageOverlayTholdDirPos,...
	'callback',{@l_OverlayControlCB,'th_direction'});
					
  % Clim edit boxes
  tmp = get(th_direction,'position');
  clim_tx = uicontrol('parent',uipanel_h,...
                      'units','pixel',...
                      'position',[10 tmp(2)-tmp(4)-10 120 15],...
                      'style','text',...
                      'string','Range (min-max):',...
                      'horizontalalign','left');
  tmp=get(clim_tx,'position');
  climmin_h = uicontrol('parent',uipanel_h,...
                        'units','pixel',...
                        'position',[tmp(1)+tmp(3) tmp(2)-2 55 20],...
                        'style','edit',...
                        'string',num2str(Dat.ImageOverlayClim(1)),...
                        'backgroundcolor','w',...
                        'userdata',Dat.ImageOverlayClim(1),...
                        'callback',{@l_OverlayControlCB,'climmin_edit'});
  tmp=get(climmin_h,'position');
  climmax_h = uicontrol('parent',uipanel_h,...
                        'units','pixel',...
                        'position',[tmp(1)+tmp(3)+5 tmp(2) 55 20],...
                        'style','edit',...
                        'string',num2str(Dat.ImageOverlayClim(2)),...
                        'backgroundcolor','w',...
                        'userdata',Dat.ImageOverlayClim(2),...
                        'callback',{@l_OverlayControlCB,'climmax_edit'});
  
  
  % Transparency controls
  tmp = get(th_slider,'position');
  tmp2 = get(clim_tx,'position');
  alpha_tx = uicontrol('parent',uipanel_h,...
                       'units','pixel',...
                       'position',[tmp(1) tmp2(2)-30 150 15],...
                       'style','text',...
                       'string','Transparency',...
                       'horizontalalign','left');
  tmp = get(alpha_tx,'position');
  alpha_slider = uicontrol('parent',uipanel_h,...
                           'units','pixel',...
                           'position',[tmp(1) tmp(2)-20 200 20],...
                           'style','slider',...
                           'Min',0,...
                           'Max',1,...
                           'SliderStep',[0.1 0.3],...
                           'value',Dat.ImageOverlayAlpha,...
                           'backgroundcolor',GD.col.slider);
  if ~GD.HG2graphics
		alphaSliderListener = handle.listener(alpha_slider,...
			'ActionEvent',...
			{@l_OverlayControlCB,'alpha_slider'});
		setappdata(alpha_slider,'alphaSliderListener',...
			alphaSliderListener)
	else
		alphaSliderListener = addlistener(alpha_slider,...
			'ContinuousValueChange',...
			@l_OverlayControlCB);
		setappdata(alpha_slider,'alphaSliderListener',...
			alphaSliderListener)
	end
  
  tmp = get(alpha_slider,'position');
  alpha_edit = uicontrol('parent',uipanel_h,...
                         'units','pixel',...
                         'position',[tmp(1)+tmp(3)+5 tmp(2) 35 20],...
                         'style','edit',...
                         'string',num2str(Dat.ImageOverlayAlpha),...
                         'backgroundcolor','w',...
						 'userdata',Dat.ImageOverlayAlpha,...
                         'callback',{@l_OverlayControlCB,'alpha_edit'});
  
  % Colormap controls
  colmap_tx = uicontrol('parent',uipanel_h,...
                        'units','pixel',...
                        'position',[tmp(1) tmp(2)-20 150 15],...
                        'style','text',...
                        'string','Colormap',...
                        'horizontalalign','left');
  tmp = get(colmap_tx,'position');
  colmap_popup = uicontrol('parent',uipanel_h,...
                           'units','pixel',...
                           'position',[tmp(1) tmp(2)-20 100 20],...
                           'style','popup',...
                           'string',{'Gray','Jet','HSV','Cool','Hot','Spring',...
                      'Summer','Autumn','Winter','Bone','Copper','Pink'},...
                           'backgroundcolor','w',...
                           'callback',{@l_OverlayControlCB,'colmap'});
  set(colmap_popup,'value',...
                   find(strcmpi(Dat.ImageOverlayCmapStr,get(colmap_popup,'string'))));
  
  % Hide image overlay
  tmp = get(colmap_popup,'position');
  hide_chbox = uicontrol('parent',uipanel_h,...
                         'units','pixel',...
                         'position',[tmp(1) tmp(2)-25 150 20],...
                         'style','checkbox',...
                         'string','Hide Image Overlay',...
                         'value',0,...
                         'callback',{@l_OverlayControlCB,'hide_overlay'});
  if Dat.isImageOverlayHidden
    set(hide_chbox,'value',1)
  else
    set(hide_chbox,'value',0)
  end
  
  if Dat.isOverlayRGB
	set([th_slider,th_edit,th_direction,climmin_h,climmax_h,colmap_popup],...
	  'enable','off')
  end
  
elseif strcmpi(opt,'delete')
  
  % Ask confirmation
  resp = questdlg('Are you sure you want to delete image overlay?',...
                  'Delete image overlay?','Yes','No','No');
  if isempty(resp) || strcmpi(resp,'No')
    return
  end
  
  % Delete image overlay
  Dat.ImageOverlay = [];
	Dat.ImageOverlayOrig = [];
  Dat.isImageOverlayLoaded = false;
  
  % Delete axes and image objects
  delete(H.IMOVERLAY)
  H.IMOVERLAY = [];
  
  % Delete overlay control window
  if isfield(H,'OVERLAY_CONTROL_FIG') && ~isempty(H.OVERLAY_CONTROL_FIG) && ishandle(H.OVERLAY_CONTROL_FIG)
    delete(H.OVERLAY_CONTROL_FIG)
  end
  
  % Disable overlay uimenus
  set([H.UIOVERLAY_CONTROLS,...
		H.UIOVERLAY_SAVE,...
		H.UIOVERLAY_DELETE],'enable','off')
  
  % Refresh crossbars
  l_ShowHideCrossbars([],[])
  
elseif strcmpi(opt,'save')
	% Save overlay to file ------------------------------
	
	if isempty(Dat.ImageOverlayOrig)
		errordlg('Cannot save overlay. Overlay is not loaded.',...
			'Cannot save overlay.','modal');
		return
	end
	
	% Get default save directory
	if ispref('Aedes','PutOverlayFileDir')
		default_dir = getpref('Aedes','PutOverlayFileDir');
	else
		default_dir = [pwd,filesep];
	end
	default_fname = 'overlay';
	
	% Prompt for a file
	[fname,fpath,findex] = uiputfile({'*.nii','NIfTI files (*.nii)';...
		'*.mat','Matlab MAT files (*.mat)'},...
		'Save overlay as...',[default_dir,default_fname]);
	if isequal(fname,0)
		% Canceled
		return
	end
	
	% Check if NIfTI or MAT file should be saved
	if findex == 1
		saveformat = 'nifti';
	else
		saveformat = 'mat';
	end
	
	switch saveformat
		case 'nifti' % -----------------------
			
			% Check filename
			if length(fname)<4 || ~strcmpi(fname(end-3:end),'.nii')
				fname = [fname,'.nii'];
			end
			
			try
				aedes_write_nifti(Dat.ImageOverlayOrig,[fpath,fname]);
			catch
				errordlg({'Error writing NIfTI file.',...
					[fpath,fname],...
					'',...
					lasterr},...
					'Could not write NIfTI file','modal');
				return
			end
			
		case 'mat' % -----------------------
			
			% Check filename
			if length(fname)<4 || ~strcmpi(fname(end-3:end),'.mat')
				fname = [fname,'.mat'];
			end
			
			try
				data.DataFormat = 'mat';
				data.HDR.fname = '';
				data.HDR.fpath = '';
				data.HDR.FileHeader.ImageOverlayThold = Dat.ImageOverlayThold;
				data.HDR.FileHeader.ImageOverlayTholdDirPos = Dat.ImageOverlayTholdDirPos;
				data.HDR.FileHeader.ImageOverlayClim = Dat.ImageOverlayClim;
				data.FTDATA = Dat.ImageOverlayOrig;
				data.KSPACE = [];
				data.PROCPAR = [];
				data.PHASETABLE = [];
				save([fpath,fname],'data','-mat');
				clear data
			catch
				errordlg({'Error saving overlay to MAT file.',...
					[fpath,fname],...
					'',...
					lasterr},...
					'Could not write MAT file','modal');
				return
			end
			
		otherwise
			error('Unknown save format for overlay "%s".',saveformat);
	end
	
	% Put save directory to preferences
	try
		setpref('Aedes','PutOverlayFileDir',fpath);
	catch
		warning('Could not save overlay directory preferences.');
	end
	
end
  
%catch
%  aedes_errordump(lasterror);
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overlay controls callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function l_OverlayControlCB(h,evd,opt)
  try
		
	if ~isempty(h)
		if h==th_slider
			opt = 'th_slider';
		elseif h== alpha_slider
			opt = 'alpha_slider';
		end
	end
	
  % Threshold slider callback
  if strcmpi(opt,'th_slider')
    val=get(th_slider,'value');
    val=round(val)*((Dat.ImOverlayMax-Dat.ImOverlayMin)/256)+Dat.ImOverlayMin;
	val=fix(val*100)/100;
    set(th_edit,'string',num2str(val))
    Dat.ImageOverlayThold = val;
    
    % Alpha slider callback
  elseif strcmpi(opt,'alpha_slider')
    val=get(alpha_slider,'value');
    val=fix(val*100)/100;
    set(alpha_edit,'string',num2str(val))
    Dat.ImageOverlayAlpha = val;
	setpref('Aedes','ImageOverlayAlpha',val)
    
	% Threshold direction callback
  elseif strcmpi(opt,'th_direction')
	val=get(th_direction,'value');
	Dat.ImageOverlayTholdDirPos = val;
	setpref('Aedes','ImageOverlayTholdDirPos',val)
	
	% Threshold editbox callback
  elseif strcmpi(opt,'th_edit')
    
    str = get(th_edit,'string');
    val=str2num(str);
    if isempty(val) || ~isreal(val) || val<Dat.ImOverlayMin || val>Dat.ImOverlayMax
      if val<Dat.ImOverlayMin
        val = Dat.ImOverlayMin;
        value=round(((val-Dat.ImOverlayMin)*256)/(Dat.ImOverlayMax-Dat.ImOverlayMin));
        set(th_slider,'value',value)
        Dat.ImageOverlayThold = val;
      elseif val>Dat.ImOverlayMax
        val = Dat.ImOverlayMax;
        value=round(((val-Dat.ImOverlayMin)*256)/(Dat.ImOverlayMax-Dat.ImOverlayMin));
        set(th_slider,'value',value)
        Dat.ImageOverlayThold = val;
      else
        set(th_edit,'string',num2str(Dat.ImageOverlayThold))
        set(th_slider,'value',...
          round(((Dat.ImageOverlayThold-Dat.ImOverlayMin)*256)/(Dat.ImOverlayMax-Dat.ImOverlayMin)))
      end
    else
      value=round(((val-Dat.ImOverlayMin)*256)/(Dat.ImOverlayMax-Dat.ImOverlayMin));
      set(th_slider,'value',value)
      Dat.ImageOverlayThold = val;
    end
    
    % Alpha editbox callback
  elseif strcmpi(opt,'alpha_edit')
    
    str = get(alpha_edit,'string');
    val=str2num(str);
    if isempty(val) || ~isreal(val) || val<0 || val>1
	  if val<0
		Dat.ImageOverlayAlpha = 0;
		set(alpha_edit,'string',num2str(Dat.ImageOverlayAlpha))
		set(alpha_slider,'value',Dat.ImageOverlayAlpha)
	  elseif val>1
		Dat.ImageOverlayAlpha = 1;
		set(alpha_edit,'string',num2str(Dat.ImageOverlayAlpha))
		set(alpha_slider,'value',Dat.ImageOverlayAlpha)
	  else
		set(alpha_edit,'string',num2str(Dat.ImageOverlayAlpha))
		set(alpha_slider,'value',Dat.ImageOverlayAlpha)
	  end
    else
      set(alpha_slider,'value',val)
      Dat.ImageOverlayAlpha = val;
	  setpref('Aedes','ImageOverlayAlpha',val)
    end
    
  elseif strcmpi(opt,'hide_overlay')

    if get(hide_chbox,'value')
      Dat.isImageOverlayHidden = true;
      set(H.IMOVERLAY,'visible','off')
    else
      Dat.isImageOverlayHidden = false;
    end
    
    
  elseif strcmpi(opt,'colmap')
    
    str=lower(popupstr(colmap_popup));
    Dat.ImageOverlayCmapStr = str;
    eval(['Dat.ImageOverlayCmap=',str,'(256);']);
    setpref('Aedes','ImageOverlayCmap',str)
	
  elseif strcmpi(opt,'close_window')
	tmp_pos = get(H.OVERLAY_CONTROL_FIG,'position');
    setpref('Aedes','OverlayControlPos',[tmp_pos(1) tmp_pos(2)]);
	delete(H.OVERLAY_CONTROL_FIG)
  H.OVERLAY_CONTROL_FIG = [];
	
  elseif strcmpi(opt,'climmin_edit')
    val = get(climmin_h,'string');
    if strcmpi(val,'min')
      val = Dat.ImOverlayMin;
      val=fix(val*100)/100;
      set(climmin_h,'string',num2str(val))
    else
      val = str2num(val);
      val=fix(val*100)/100;
    end
    
    val_2 = get(climmax_h,'userdata');
    if isempty(val) || ~isreal(val) || ...
        val>=val_2
      set(climmin_h,'string',num2str(get(climmin_h,'userdata')))
%     elseif val<Dat.ImOverlayMin
%       set(climmin_h,'userdata',Dat.ImOverlayMin)
%       Dat.ImageOverlayClim(1)=Dat.ImOverlayMin;
%       set(climmin_h,'string',num2str(Dat.ImOverlayMin))
    else
      set(climmin_h,'userdata',val)
      Dat.ImageOverlayClim(1)=val;
    end
  elseif strcmpi(opt,'climmax_edit')
    val = get(climmax_h,'string');
    if strcmpi(val,'max')
      val = Dat.ImOverlayMax;
      val=fix(val*100)/100;
      set(climmax_h,'string',num2str(val))
    else
      val = str2num(val);
      val=fix(val*100)/100;
    end
    val_2 = get(climmin_h,'userdata');
    if isempty(val) || ~isreal(val) || ...
        val<=val_2
      set(climmax_h,'string',num2str(get(climmax_h,'userdata')))
%     elseif val>Dat.ImOverlayMax
%       set(climmax_h,'userdata',Dat.ImOverlayMax)
%       Dat.ImageOverlayClim(2)=Dat.ImOverlayMax;
%       set(climmax_h,'string',num2str(Dat.ImOverlayMax))
    else
      set(climmax_h,'userdata',val)
      Dat.ImageOverlayClim(2)=val;
    end
  end
  
  % Refresh screen
  if Dat.AxView==0
    l_DisplayData([],[])
  else
    l_DisplayData([],[],Dat.AxView)
  end
  
  catch
	aedes_errordump(lasterror);
  end
  end % function l_OverlayControlCB(h,


end % function l_OverlayControls(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Update info text
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_UpdateInfoText(h,evd)
try
% Called from uibutton or uimenu
if ~isempty(h)
  if strcmpi(get(H.UITOGGLE_INFO,'state'),'on')
    set(H.INFOTEXT,'visible','on')
    l_UpdateRoiInfoText([],[])
    set(H.ROIINFOTEXT,'visible','on')
    setpref('Aedes','ShowInfoTexts',true);
  else
    set([H.INFOTEXT,H.ROIINFOTEXT],'visible','off')
    setpref('Aedes','ShowInfoTexts',false);
  end
  %return
elseif strcmpi(get(H.UITOGGLE_INFO,'state'),'off')
  return
end

if Dat.ZoomLevel==0
  zoom_txt = 'normalized';
else
  zoom_txt = [num2str(Dat.ZoomLevel),'x'];
end

% file_str = fliplr(sprintf('%12s',fliplr('File:')));
% slice_str = fliplr(sprintf('%12s',fliplr('Slice:')));
% imagesize_str = fliplr(sprintf('%12s',fliplr('Image Size:')));
% zoom_str = fliplr(sprintf('%12s',fliplr('Zoom:')));
% rotation_str = fliplr(sprintf('%12s',fliplr('Rotation:')));
% flip_str = fliplr(sprintf('%12s',fliplr('Flip:')));

if Dat.isDataMixed
  % Data rotation text
  if Dat.DataRotation(Dat.DataInd)==0
    rot_txt = '-';
  elseif Dat.DataRotation(Dat.DataInd)==1
    rot_txt = ['90',char(176)];
  elseif Dat.DataRotation(Dat.DataInd)==2
    rot_txt = ['180',char(176)];
  else
    rot_txt = ['270',char(176)];
  end
  
  % Data flipping text
  if Dat.DataFlip(Dat.DataInd)==0
    flip_txt = '-';
  elseif Dat.DataFlip(Dat.DataInd)==1
    flip_txt = 'UD';
  elseif Dat.DataFlip(Dat.DataInd)==2
    flip_txt = 'LR';
  end
  
  
%   str=sprintf('%s%s%s\n%s%d/%d\%s%dx%d\n%s%s\n%s%s\n%s%s',...
% 	file_str,DATA{Dat.DataInd}.HDR.fpath,DATA{Dat.DataInd}.HDR.fname,...
% 	slice_str,Dat.DataInd,Dat.DataL,...
% 	imagesize_str,Dat.ImageDim(Dat.DataInd,1),Dat.ImageDim(Dat.DataInd,2),...
% 	zoom_str,zoom_txt,...
% 	rotation_str,rot_txt,...
% 	flip_str,flip_txt);
%   
%   set(H.INFOTEXT,'string',str);
  set(H.INFOTEXT,'string',...
                 {['File:       ',DATA{Dat.DataInd}.HDR.fpath,...
                   DATA{Dat.DataInd}.HDR.fname],...
                  ['Slice:      ',num2str(Dat.DataInd),'/',num2str(Dat.DataL)],...
                  ['Image Size: ',num2str(Dat.ImageDim(Dat.DataInd,1)),'x',...
                   num2str(Dat.ImageDim(Dat.DataInd,2))],...
                  ['Zoom:       ',zoom_txt],...
                  ['Rotation:   ',rot_txt],...
                  ['Flip:       ',flip_txt]});
else
  if Dat.ImageDim(Dat.DataInd,4)>1
	set(H.INFOTEXT,'string',...
	  {['File:       ',DATA{Dat.DataInd}.HDR.fpath,...
	  DATA{Dat.DataInd}.HDR.fname],...
	  ['Volume:     ',...
	  num2str(Dat.CurrentVol),'/',num2str(Dat.ImageDim(Dat.DataInd,4))],...
	  ['Image Size: ',num2str(Dat.ImageDim(Dat.DataInd,1)),'x',...
	  num2str(Dat.ImageDim(Dat.DataInd,2)),'x',...
	  num2str(Dat.ImageDim(Dat.DataInd,3)),'x',...
	  num2str(Dat.ImageDim(Dat.DataInd,4))],...
	  ['Zoom:       ',zoom_txt]});
  else
	set(H.INFOTEXT,'string',...
	  {['File:       ',DATA{Dat.DataInd}.HDR.fpath,...
	  DATA{Dat.DataInd}.HDR.fname],...
	  ['Slice:      ',...
	  num2str(Dat.Slices(3)),'/',num2str(Dat.ImageDim(Dat.DataInd,3))],...
	  ['Image Size: ',num2str(Dat.ImageDim(Dat.DataInd,1)),'x',...
	  num2str(Dat.ImageDim(Dat.DataInd,2)),'x',...
	  num2str(Dat.ImageDim(Dat.DataInd,3))],...
	  ['Zoom:       ',zoom_txt]});
  end
end


catch
  aedes_errordump(lasterror);
end
end % function l_UpdateInfoText(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Update ROI Info Text
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_UpdateRoiInfoText(h,evd)
try
if strcmpi(get(H.UITOGGLE_INFO,'state'),'off')
  return
end

roi_ind = get(H.ROI_EDIT,'value');

if isempty(ROI)
  set(H.ROIINFOTEXT,'string',...
                    {'Current ROI:',...
                     'Mean:',...
                     'STD:',...
                     'Sum:',...
                     'Min:',...
                     'Max:',...
                     'Pixel Count:'})
else
  vals = DATA{Dat.DataInd}.FTDATA(ROI(roi_ind).voxels{Dat.DataInd});
  vals = double(vals);
  if isempty(vals)
    set(H.ROIINFOTEXT,'string',...
                      {['Current ROI: ',ROI(roi_ind).label],...
                       ['Mean:        ',...
                        '-'],...
                       ['STD:         ',...
                        '-'],...
                       ['Sum:         ',...
                        '-'],...
                       ['Min:         ',...
                        '-'],...
                       ['Max:         ',...
                        '-'],...
                       ['Pixel Count: ',...
                        '0']})
  else
    set(H.ROIINFOTEXT,'string',...
                      {['Current ROI: ',ROI(roi_ind).label],...
                       ['Mean:        ',...
                        sprintf('%.2f',mean(vals))],...
                       ['STD:         ',...
                        sprintf('%.2f',std(vals))],...
                       ['Sum:         ',...
                        sprintf('%.2f',sum(vals))],...
                       ['Min:         ',...
                        sprintf('%.2f',min(vals))],...
                       ['Max:         ',...
                        sprintf('%.2f',max(vals))],...
                       ['Pixel Count: ',...
                        sprintf('%.0f',numel(vals))]})
  end
end

catch
  aedes_errordump(lasterror);
end
end % function l_UpdateRoiInfoText(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Show/Hide Info Text Background
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_InfoTextBackground(h,evd)
try

if strcmpi(get(H.UIVIEW_INFOTXTBACKGROUND,'checked'),'off')
  set([H.INFOTEXT,H.ROIINFOTEXT],...
                    'backgroundcolor','w','color','k')
	if ~Dat.HG2graphics
		set([H.INFOTEXT,H.ROIINFOTEXT],'erasemode','normal')
	end
  set(H.UIVIEW_INFOTXTBACKGROUND,'checked','on')
  setpref('Aedes','ShowInfoTextBackGround',true)
else
  set([H.INFOTEXT,H.ROIINFOTEXT],...
                    'backgroundcolor','none','color','w')
	if ~Dat.HG2graphics
		set([H.INFOTEXT,H.ROIINFOTEXT],'erasemode','normal')
	end
  set(H.UIVIEW_INFOTXTBACKGROUND,'checked','off')
  setpref('Aedes','ShowInfoTextBackGround',false)
end

catch
  aedes_errordump(lasterror);
end
end % function l_InfoTextBackground(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Show/hide colorbar
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_ShowColorbar(h,evd)
try
if strcmpi(get(h,'state'),'on')
  Dat.ShowColorbar = true;
  setpref('Aedes','ShowColorbar',true)
  vis= 'on';
else
  Dat.ShowColorbar = false;
  setpref('Aedes','ShowColorbar',false)
  vis = 'off';
end

%% Refresh intensity value
l_UpdateIntensityValue([],[])

l_AxesPositions([],[],Dat.AxView,Dat.ZoomLevel)
set([H.COLORBAR_AX,H.COLORBAR_IM,],'visible',vis)
if isfield(H,'COLORBAR_LN') && isfield(H,'COLORBAR_TX') && ...
    ~isempty(H.COLORBAR_LN) && ~isempty(H.COLORBAR_TX)
  set([H.COLORBAR_LN,H.COLORBAR_TX],'visible',vis)
end

catch
  aedes_errordump(lasterror);
end
end % function l_ShowColorbar(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Calculate axis positions
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_AxesPositions(h,evd,ax_ind,zoom_level)
try
fig_pos = get(H.FIG,'position');
%sidebar_pos = 232;
largegap=Dat.AxGap;
smallgap = 4;

if Dat.ShowColorbar
  colorbar_gap = 80;
else
  colorbar_gap = 0;
end
ax_size=round(Dat.RealAxSize{Dat.DataInd}*Dat.ZoomLevel);

% If called from image sliders, use current values for input arguments
% ax_ind and ax_size
if ~isempty(h) && ishandle(h) && any(h==H.IMSLIDER)
  ImSliderCB = 1;
else
  ImSliderCB = 0;
end
if ImSliderCB
  ax_ind = Dat.AxView;
  ax_size = round(Dat.RealAxSize{Dat.DataInd}*Dat.ZoomLevel);
end

if ax_ind~=0
	tmp_ind = [3 2 1];
	ax_ind = tmp_ind(ax_ind);
end

% Axes positions are normalized to the figure size
if all(ax_size==0)
  
  % Get axes positions without sliders
  if ax_ind==0
    ax_h1 = round((fig_pos(4)-smallgap-2*largegap)*...
                    (Dat.RealAxSize{Dat.DataInd}(1,2)/...
                     (Dat.RealAxSize{Dat.DataInd}(1,2)+Dat.RealAxSize{Dat.DataInd}(3,2))));
    prop_sz = ax_h1/Dat.RealAxSize{Dat.DataInd}(1,2);
  else
    prop_sz = (fig_pos(4)-largegap-smallgap)/Dat.RealAxSize{Dat.DataInd}(ax_ind,2);
  end
  axpos=l_axpos(ax_ind,ax_size,0,0,prop_sz);
  
  % Determine the need of image sliders
  if ~ImSliderCB
    isYSlider = 0;
    if ax_ind==0
      x_over = fig_pos(3)-232-2*largegap-smallgap-axpos(1,3)-axpos(2,3)-colorbar_gap;
    else
      x_over = fig_pos(3)-232-largegap-smallgap-axpos(ax_ind,3)-colorbar_gap;
    end
    y_over = 1;
    if x_over < 0
      isXSlider = 1;
    else
      isXSlider = 0;
    end
  else
    isXSlider = 1;
    y_over = 1;
    x_over = get(H.IMSLIDER(1),'Max');
  end
  
  % Calculate proportional sizes to axes
  if ax_ind==0
    if isXSlider
      ax_h1 = round((fig_pos(4)-2*largegap-smallgap-17)*...
                    (Dat.RealAxSize{Dat.DataInd}(1,2)/...
                     (Dat.RealAxSize{Dat.DataInd}(1,2)+Dat.RealAxSize{Dat.DataInd}(3,2))));
    else
      ax_h1 = round((fig_pos(4)-2*largegap-smallgap)*...
                    (Dat.RealAxSize{Dat.DataInd}(1,2)/...
                     (Dat.RealAxSize{Dat.DataInd}(1,2)+Dat.RealAxSize{Dat.DataInd}(3,2))));
    end
    prop_sz = ax_h1/Dat.RealAxSize{Dat.DataInd}(1,2);
  else
    if isXSlider
      prop_sz = (fig_pos(4)-largegap-smallgap-17)/Dat.RealAxSize{Dat.DataInd}(ax_ind,2);
    else
      prop_sz = (fig_pos(4)-largegap-smallgap)/Dat.RealAxSize{Dat.DataInd}(ax_ind,2);
    end
  end
  
  % Set image sliders if necessary
  if ~ImSliderCB
    if isXSlider
      l_SetImageSliders([],[],isXSlider,isYSlider,x_over,y_over)
    else
      set(H.IMSLIDER,'visible','off')
      set(H.IMSLIDER_FRAME,'visible','off')
    end
  end
  
  % Get image slider value
  if isXSlider
    x=get(H.IMSLIDER(1),'value');
  else
    x=0;
  end
  y=0;
  
  % Get axes positions
  axpos=l_axpos(ax_ind,ax_size,x,y,prop_sz);
    
else
  
  % Get axes positions without sliders
  axpos=l_axpos(ax_ind,ax_size,0,0);
  
  % Determine the need of image sliders
  if ax_ind==0
    x_over = fig_pos(3)-232-2*largegap-smallgap-axpos(1,3)-axpos(2,3)-colorbar_gap;
    y_over = fig_pos(4)-axpos(1,4)-axpos(3,4)-2*largegap-smallgap;
  else
    x_over = fig_pos(3)-232-largegap-smallgap-axpos(ax_ind,3)-colorbar_gap;
    y_over = fig_pos(4)-axpos(ax_ind,4)-largegap-smallgap;
  end
  isXSlider = x_over < 0;
  isYSlider = y_over < 0;
  if isXSlider && ~isYSlider
    if y_over-17-largegap < 0
      y_over = y_over-17-largegap;
      x_over = x_over-17-smallgap;
      isYSlider = 1;
    end
  elseif isYSlider && ~isXSlider
    if x_over-17-smallgap < 0
      x_over = x_over-17-smallgap;
      y_over = y_over-17-largegap;
      isXSlider = 1;
    end
  elseif isYSlider && isXSlider
    y_over = y_over-17-largegap;
    x_over = x_over-17-smallgap;
  end
  
  
  % Set image sliders if necessary
  if ~ImSliderCB
    if isXSlider || isYSlider
      l_SetImageSliders([],[],isXSlider,isYSlider,x_over,y_over)
    else
      set(H.IMSLIDER,'visible','off')
      set(H.IMSLIDER_FRAME,'visible','off')
    end
  end
  
  % Get image slider value
  if isXSlider
    x=get(H.IMSLIDER(1),'value');
  else
    x=0;
  end
  if isYSlider
    y=get(H.IMSLIDER(2),'value');
  else
    y=0;
  end
  
  % Get axes positions
  axpos=l_axpos(ax_ind,ax_size,x,y);
  
end
%axpos
% Set new positions to the axes
set(H.IMAX1,'position',axpos(1,:))
set(H.IMAX2,'position',axpos(2,:))
set(H.IMAX3,'position',axpos(3,:))
set(H.IMOVERLAYAX(1),'position',axpos(1,:))
set(H.IMOVERLAYAX(2),'position',axpos(2,:))
set(H.IMOVERLAYAX(3),'position',axpos(3,:))
if Dat.ShowColorbar
  set(H.COLORBAR_AX,'position',axpos(4,:))
end
if ~isempty(H.ROIAX)
  set(H.ROIAX(1,:),'position',axpos(1,:))
  set(H.ROIAX(2,:),'position',axpos(2,:))
  set(H.ROIAX(3,:),'position',axpos(3,:))
end
%drawnow

catch
  aedes_errordump(lasterror);
end
end % function l_AxesPositions(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Helper function for l_AxesPositions
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function axpos=l_axpos(ax_ind,ax_size,x,y,prop_sz)
try
axpos = zeros(3,4);

fig_pos = get(H.FIG,'position');
largegap=Dat.AxGap;
colbar_val = 20;
smallgap = 4;
if Dat.ShowColorbar
  tightinset = get(H.COLORBAR_AX,'TightInset');
  ttop = tightinset(4);
  if ttop<10
    ttop=0;
  end
end


%% Calculate axes positions
if all(ax_size==0) % Axes sizes are normalized
  
  % Positions for 3D-view ---------------------
  if ax_ind==0
    % Calculate axes positions
    % X-axis [left,bottom,width,height]
    axpos(1,:) = [largegap+232,...
                  fig_pos(4)-prop_sz*Dat.RealAxSize{Dat.DataInd}(1,2)-smallgap,...
                  prop_sz*Dat.RealAxSize{Dat.DataInd}(1,1),...
                  prop_sz*Dat.RealAxSize{Dat.DataInd}(1,2)];
    % Y-axis [width,height]
    axpos(2,:) = [2*largegap+prop_sz*Dat.RealAxSize{Dat.DataInd}(1,1)+232,...
                  fig_pos(4)-prop_sz*Dat.RealAxSize{Dat.DataInd}(2,2)-smallgap,...
                  prop_sz*Dat.RealAxSize{Dat.DataInd}(2,1),...
                  prop_sz*Dat.RealAxSize{Dat.DataInd}(2,2)];
    % Z-axis [width,height]
    axpos(3,:) = [largegap+232,...
                  fig_pos(4)-largegap-smallgap-prop_sz*Dat.RealAxSize{Dat.DataInd}(1,2)-prop_sz*Dat.RealAxSize{Dat.DataInd}(3,2),...
                  prop_sz*Dat.RealAxSize{Dat.DataInd}(3,1),...
                  prop_sz*Dat.RealAxSize{Dat.DataInd}(3,2)];
    
    % Colorbar position
    if Dat.ShowColorbar
      axpos(4,:) = [2*largegap+prop_sz*Dat.RealAxSize{Dat.DataInd}(1,1)+232+...
                    prop_sz*Dat.RealAxSize{Dat.DataInd}(2,1)+5,...
                    smallgap+colbar_val,...
                    35*fig_pos(4)/661,...
                    fig_pos(4)-2*smallgap-colbar_val-ttop];
    end
  else % Positions for X, Y, or Z views ------------------------
    
    % Calculate axes positions.
    % Set the position of the viewed axes
    axpos(ax_ind,:) = [largegap+232,... 
                       fig_pos(4)-prop_sz*Dat.RealAxSize{Dat.DataInd}(ax_ind,2)-smallgap,...
                       prop_sz*Dat.RealAxSize{Dat.DataInd}(ax_ind,1),...
                       prop_sz*Dat.RealAxSize{Dat.DataInd}(ax_ind,2)];
    
    % Set all other axis to ones
    tmp = [1 2 3];
    axpos(tmp(tmp~=ax_ind),:) = ones(2,4);
    
    % Colorbar position
    if Dat.ShowColorbar
      axpos(4,:) = [largegap+232+prop_sz*Dat.RealAxSize{Dat.DataInd}(ax_ind,1)+5,...
                    smallgap+colbar_val,...
                    35*fig_pos(4)/661,...
                    fig_pos(4)-2*smallgap-colbar_val-ttop];
    end
  end
  
else % Axes sizes are determined with pixels
  
  %% View only X, Y, or Z
  if any(ax_ind==[1 2 3])
    
    % Calculate axes positions.
    % Set the position of the viewed axes
    axpos(ax_ind,:) = [largegap+232,...
                       fig_pos(4)-smallgap-ax_size(ax_ind,2),...
                       ax_size(ax_ind,1),...
                       ax_size(ax_ind,2)];
    
    % Set all other axis to ones
    tmp = [1 2 3];
    axpos(tmp(tmp~=ax_ind),:) = ones(2,4);
    
    % Colorbar position
    if Dat.ShowColorbar
      axpos(4,:) = [largegap+232+ax_size(ax_ind,1)+5,...
                    fig_pos(4)-smallgap-ax_size(ax_ind,2)+colbar_val,...
                    35*ax_size(ax_ind,2)/661,...
                    ax_size(ax_ind,2)-colbar_val-ttop];
      if axpos(4,4)<1
        axpos(4,2)=axpos(ax_ind,2);
        axpos(4,4)=axpos(ax_ind,4);
      end
    end
    
  else % View 3D
    
    % Calculate axes positions
    % X-axis [left,bottom,width,height]
    axpos(1,:) = [largegap+232,...
                  fig_pos(4)-ax_size(1,2)-smallgap,...
                  ax_size(1,1),...
                  ax_size(1,2)];
    % Y-axis [width,height]
    axpos(2,:) = [2*largegap+ax_size(1,1)+232,...
                  fig_pos(4)-ax_size(2,2)-smallgap,...
                  ax_size(2,1),...
                  ax_size(2,2)];
    % Z-axis [width,height]
    axpos(3,:) = [largegap+232,...
                  fig_pos(4)-largegap-smallgap-ax_size(1,2)-ax_size(3,2),...
                   ax_size(3,1),...
                  ax_size(3,2)];
    
    % Colorbar position
    if Dat.ShowColorbar
      axpos(4,:) = [2*largegap+ax_size(1,1)+232+ax_size(2,1)+5,...
                    fig_pos(4)-largegap-smallgap-ax_size(1,2)-ax_size(3,2)+colbar_val,...
                    35*(ax_size(1,2)+ax_size(3,2)+smallgap)/661,...
                    ax_size(1,2)+ax_size(3,2)+smallgap-colbar_val-ttop];
      
    end
  end
  
end

% Move axes according to image sliders
axpos(:,1)=axpos(:,1)-x;
axpos(:,2)=axpos(:,2)-y;
return

catch
  aedes_errordump(lasterror);
end
end % function axpos=l_axpos(ax_ind,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Change the displayed slice
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_ChangeSlice(h,evd,opt,add_input)
try
% Get indices to updated images
update_ind = [1 2 3];

% Chewing gum and ironwire for HG2
if ~isempty(h)
	if h == H.SL_SLIDER
		opt = 'slider';
	end
end

% Check if function is called from image buttondownfcn
switch opt
 case 'mouse' % ------------------------------------------
  
  % Get indice to image
  ind=[3 2 1];
  update_ind=update_ind(update_ind~=add_input);
  
  % Get handles to the image axes
  ax_h = get(H.IM(ind(add_input)),'parent');
  cp = get(ax_h,'currentpoint');
	
  % Set slice indices
  old_slices = Dat.Slices;
  if add_input==3
    yy=min(max(round(cp(1)),1),Dat.ImageDim(Dat.DataInd,2));
    xx=min(max(round(cp(3)),1),Dat.ImageDim(Dat.DataInd,1));
    Dat.Slices = [xx yy old_slices(add_input)];
  elseif add_input==2
    yy=min(max(round(cp(1)),1),Dat.ImageDim(Dat.DataInd,3));
    xx=min(max(round(cp(3)),1),Dat.ImageDim(Dat.DataInd,1));
    Dat.Slices = [xx old_slices(add_input) yy];
  else
    yy=min(max(round(cp(1)),1),Dat.ImageDim(Dat.DataInd,3));
    xx=min(max(round(cp(3)),1),Dat.ImageDim(Dat.DataInd,2));
    Dat.Slices = [old_slices(add_input) xx yy];
	end
	
 case 'wheel' % -----------------------------------------
	 
	 if ~Dat.isDataMixed
		 dim = evd;
		 if isempty(H.ROIAX)
			 roiax = [0;0;0];
		 else
			 roiax = H.ROIAX;
		 end
		 if Dat.AxView==0
			 if dim == 3
				 Dat.Slices(3) = min(max(Dat.Slices(3)+add_input,1),Dat.ImageDim(Dat.DataInd,3));
			 elseif dim == 2
				 Dat.Slices(2) = min(max(Dat.Slices(2)+add_input,1),Dat.ImageDim(Dat.DataInd,2));
			 else
				 Dat.Slices(1) = min(max(Dat.Slices(1)+add_input,1),Dat.ImageDim(Dat.DataInd,1));
			 end
		 else
			 ind = [1 2 3];
			 Dat.Slices(Dat.AxView) = min(max(Dat.Slices(Dat.AxView)+add_input,1),Dat.ImageDim(Dat.DataInd,ind(Dat.AxView)));
		 end
	 else
		 Dat.DataInd = min(max(Dat.DataInd+add_input,1),Dat.DataL);
		 
		 % Check that the slice coordinates don't exceed the data
		 Dat.Slices(2) = min(max(Dat.Slices(2),1),Dat.ImageDim(Dat.DataInd,2));
		 Dat.Slices(1) = min(max(Dat.Slices(1),1),Dat.ImageDim(Dat.DataInd,1));
	 end
	 
 case 'slider' % -----------------------------------------
  
   % Change volume
   if get(H.V_BTN,'value')==1
     l_ChangeVolume(H.SL_SLIDER,[],'slider');
     return
   else % Change slice
     
     % Get current slider value
     sl_val = round(get(h,'value'));
     tmp=get([H.X_BTN,H.Y_BTN,H.Z_BTN],'value');
     update_ind=update_ind(logical([tmp{:}]));
     
     %% Determine which orientation to change
     if ~Dat.isDataMixed
       if get(H.X_BTN,'value')
         Dat.Slices(1) = sl_val;
       elseif get(H.Y_BTN,'value')
         Dat.Slices(2) = sl_val;
       else
         Dat.Slices(3) = sl_val;
       end
     else
       dataind = Dat.DataInd;
       Dat.DataInd = sl_val;
       
       Dat.Slices(1)=round(Dat.Slices(1)/Dat.ImageDim(dataind,1)*...
                           Dat.ImageDim(Dat.DataInd,1));
       Dat.Slices(2)=round(Dat.Slices(2)/Dat.ImageDim(dataind,2)*...
                           Dat.ImageDim(Dat.DataInd,2));
       
       
       %% Check that the slice coordinates don't exceed
       %% the data
       if Dat.Slices(3)<1
         Dat.Slices(3)=1;
       end
       if Dat.Slices(2)<1
         Dat.Slices(2)=1;
       end
       if Dat.ImageDim(Dat.DataInd,1)<Dat.Slices(1) 
         Dat.Slices(1) = Dat.ImageDim(Dat.DataInd,1);
       end
       if Dat.ImageDim(Dat.DataInd,2)<Dat.Slices(2)
         Dat.Slices(2) = Dat.ImageDim(Dat.DataInd,2);
       end
     end
   end
  
 case 'editbox' % -----------------------------------------
  
   if get(H.V_BTN,'value')==1
     update_ind=[1 2 3];
   else
	 update_ind=[1 2 3];
	 %tmp=get([H.X_BTN,H.Y_BTN,H.Z_BTN],'value');
     %update_ind=update_ind(logical([tmp{:}]));
   end
  
  %% Check that all values are within slice limits
  x_val = str2num(get(H.X_EDIT,'string'));
  y_val = str2num(get(H.Y_EDIT,'string'));
  z_val = str2num(get(H.Z_EDIT,'string'));
  
  % Check Z
  if Dat.isDataMixed
    z_valnew=min(max(z_val,1),Dat.DataL);
    Dat.DataInd = z_valnew;
    set(H.Z_EDIT,'string',num2str(z_valnew))
  else
    z_valnew=min(max(z_val,1),Dat.ImageDim(Dat.DataInd,3));
    set(H.Z_EDIT,'string',num2str(z_valnew))
  end
  
  % Check Y
  y_valnew=min(max(y_val,1),Dat.ImageDim(Dat.DataInd,2));
  set(H.Y_EDIT,'string',num2str(y_valnew))
  
  % Check X
  x_valnew=min(max(x_val,1),Dat.ImageDim(Dat.DataInd,1));
  set(H.X_EDIT,'string',num2str(x_valnew))
  
  % Set slice coordinates
  if Dat.isDataMixed
    Dat.Slices = [1,1,1];
  else
    Dat.Slices = [x_valnew,y_valnew,z_valnew];
  end
end

% Update slider value
if get(H.V_BTN,'value')==0
  tmp=get([H.X_BTN,H.Y_BTN,H.Z_BTN],'value');
  if Dat.isDataMixed
    set(H.SL_SLIDER,'value',Dat.DataInd)
  else
    set(H.SL_SLIDER,'value',Dat.Slices(logical([tmp{:}])))
  end
end

% Display Data
l_DisplayData([],[],update_ind)

% Change image axes x- and y-limits if data is in mixed slice form
if Dat.isDataMixed
  set(H.IMAX1,'xlim',[0.5 Dat.ImageDim(Dat.DataInd,2)+0.5],...
              'ylim',[0.5 Dat.ImageDim(Dat.DataInd,1)+0.49])
  set(H.IMOVERLAYAX(1),'xlim',[0.5 Dat.ImageDim(Dat.DataInd,2)+0.5],...
                    'ylim',[0.5 Dat.ImageDim(Dat.DataInd,1)+0.49])
  try
    set(H.ROIAX(1,:),'xlim',[0.5 Dat.ImageDim(Dat.DataInd,2)+0.5],...
                     'ylim',[0.5 Dat.ImageDim(Dat.DataInd,1)+0.49])
  catch
  end
end

% Draw crossbars
if strcmpi(get(H.UITOGGLE_CROSSHAIRS,'state'),'on')
  l_ShowHideCrossbars;
end

% Refresh intensity value
l_UpdateIntensityValue([],[])

% Refresh time series plot
if Dat.ShowTimeSeries
  l_ShowTimeSeries([],[],'');
end

% Refresh Clim if unlocked or slice locked
if any(Dat.LockClim==[0 1 2 3]) && Dat.isDataMixed
  l_RefreshClim([],[]);
end

catch
  aedes_errordump(lasterror);
end
end % function l_ChangeSlice(h,


function l_WindowButtonUpFcn(h,evd)
try
% Stop click-and-drag
set(H.FIG,'WindowButtonMotionFcn','',...
          'WindowButtonUpFcn','')

% Set pointer style to normal (arrow)
set(H.FIG,'pointer','arrow')

l_ShowRoiEdges([],[],'')

% Set crossbar lines erasemode back to normal
%set(H.CROSSBAR_LN,'erasemode','normal')

catch
  aedes_errordump(lasterror);
end
end % function l_WindowButtonUpFcn(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change Viewed orientation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_ChangeSliderOrient(h,evd,ind)
try  
% Function changes slice slider to to work in 
% X, Y, Z, or V directions. The input argument "ind"
% determines the orientation as follows:
%
% ind=1 --> X
% ind=2 --> Y
% ind=3 --> Z
% ind=4 --> V

% If there is only one slice read to Aedes there is nothing to do here
% -> return
if Dat.ImageDim(Dat.DataInd,3)==1 && Dat.ImageDim(Dat.DataInd,4)==1
  return
end



% Button handles
btn_h=[H.X_BTN,H.Y_BTN,H.Z_BTN,H.V_BTN];
if Dat.ImageDim(Dat.DataInd,4)>1 && ind==1 && Dat.ImageDim(Dat.DataInd,3)==1 
  ind = 4;
end
set(btn_h(ind),'value',1,'fontweight','bold')
set(btn_h(btn_h~=btn_h(ind)),'value',0,...
  'fontweight','normal')

%% Determine which button was pressed
if ind==1  % X-button pressed
  % Set slider callback
	set(H.SL_SLIDER,'min',1,...
	  'max',Dat.ImageDim(Dat.DataInd,1),...
	  'sliderstep',[1/(Dat.ImageDim(Dat.DataInd,1)-1) ...
	  ceil(Dat.ImageDim(Dat.DataInd,1)/10)/(Dat.ImageDim(Dat.DataInd,1)-1)],...
	  'value',str2num(get(H.X_EDIT,'string')),...
	  'enable','on')
  
elseif ind==2  % Y-button pressed
  % Set slider callback
  set(H.SL_SLIDER,'min',1,...
                  'max',Dat.ImageDim(Dat.DataInd,2),...
                  'sliderstep',[1/(Dat.ImageDim(Dat.DataInd,2)-1) ...
                      ceil(Dat.ImageDim(Dat.DataInd,2)/10)/(Dat.ImageDim(Dat.DataInd,2)-1)],...
                  'value',str2num(get(H.Y_EDIT,'string')),...
				  'enable','on')
  
elseif ind==3  % Z-button pressed
  % Set slider callback
	if Dat.ImageDim(Dat.DataInd,3)~=1
		set(H.SL_SLIDER,'min',1,...
			'max',Dat.ImageDim(Dat.DataInd,3),...
			'sliderstep',[1/(Dat.ImageDim(Dat.DataInd,3)-1) ...
			ceil(Dat.ImageDim(Dat.DataInd,3)/10)/(Dat.ImageDim(Dat.DataInd,3)-1)],...
			'value',str2num(get(H.Z_EDIT,'string')),...
			'enable','on')
	end
elseif ind==4
  % Set slider callback
  set(H.SL_SLIDER,'min',1,...
                  'max',Dat.ImageDim(Dat.DataInd,4),...
                  'sliderstep',[1/(Dat.ImageDim(Dat.DataInd,4)-1) ...
                      ceil(Dat.ImageDim(Dat.DataInd,4)/10)/(Dat.ImageDim(Dat.DataInd,4)-1)],...
                  'value',str2num(get(H.V_EDIT,'string')),...
				  'enable','on')
end

catch
  aedes_errordump(lasterror);
end
end % function l_ChangeSliderOrient(h,

%%%%%%%%%%%%%%%%%%%%%%%
% Change Volume
%%%%%%%%%%%%%%%%%%%%%%%
function l_ChangeVolume(h,evd,opt)
try
  
if strcmpi(opt,'slider')
  volnum = round(get(h,'value'));
  Dat.CurrentVol = volnum;
  set(H.V_EDIT,'string',num2str(volnum),...
               'userdata',num2str(volnum))
elseif strcmpi(opt,'editbox')
  str=get(h,'string');
  volnum = str2num(str);
  if isempty(volnum) || ~isreal(volnum) 
    set(h,'string',get(h,'userdata'))
    return
  end
  volnum=round(volnum);
  if volnum<1 | volnum>Dat.Vols
    set(h,'string',get(h,'userdata'))
    return
  else
    Dat.CurrentVol = volnum;
    set(h,'userdata',num2str(volnum))
  end
  if get(H.V_BTN,'value')==1
    set(H.SL_SLIDER,'value',volnum)
  end
end

% Display Data
l_DisplayData([],[])

% Refresh intensity value
l_UpdateIntensityValue([],[])

catch
  aedes_errordump(lasterror);
end
end % function l_ChangeVolume(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Change viewed axis
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_ChangeView(h,evd,vtype)
try
% Reset checked status in uimenus and toggletools
set([H.UIVIEW_3D,H.UIVIEW_X,H.UIVIEW_Y,H.UIVIEW_Z],...
    'checked','off')
set([H.UITOGGLE_VIEW3D,H.UITOGGLE_VIEWX,H.UITOGGLE_VIEWY,H.UITOGGLE_VIEWZ],...
    'state','off')
if vtype==1
  set(H.UIVIEW_X,'checked','on')
  set(H.UITOGGLE_VIEWX,'state','on')
elseif vtype==2
  set(H.UIVIEW_Y,'checked','on')
  set(H.UITOGGLE_VIEWY,'state','on')
elseif vtype==3
  set(H.UIVIEW_Z,'checked','on')
  set(H.UITOGGLE_VIEWZ,'state','on')
else
  set(H.UIVIEW_3D,'checked','on')
  set(H.UITOGGLE_VIEW3D,'state','on')
end

ind=[3 2 1];
if any(vtype==ind) 
  % View only X, Y, or Z
  Dat.AxView=vtype;
  set(H.IM(ind(vtype)),'visible','on')
  set(H.IM(ind(ind~=ind(vtype))),'visible','off')
  
  % Hide crossbars from other axis if they are selected
  if strcmpi(get(H.UITOGGLE_CROSSHAIRS,'state'),'on')
    set(H.CROSSBAR_LN(:,:),'visible','off')
    set(H.CROSSBAR_LN(ind(vtype),:),'visible','on')
  end

  % Set the slice slider to control the selected view
  %btn_h=[H.X_BTN,H.Y_BTN,H.Z_BTN];
  %set(btn_h,'value',0,'fontweight','normal')
  %set(btn_h(vtype),'value',1,'fontweight','bold')
  l_ChangeSliderOrient([],[],vtype)
  
else % View XYZ (3D -view, default)
  Dat.AxView=0;
  set(H.IM,'visible','on')
  if  strcmpi(get(H.UITOGGLE_CROSSHAIRS,'state'),'on')
    set(H.CROSSBAR_LN(:,:),'visible','on')
  end
end
setpref('Aedes','AxView',vtype)

% Save Dat structure
%setappdata(H.FIG,'Dat',Dat)

% Set axis positions
%l_ResizeFcn([],[])
l_AxesPositions([],[],Dat.AxView,Dat.ZoomLevel);

% Check if sliders are needed
%l_SetImageSliders(h,evd)

% Refresh Images and ROI(s)
if vtype==0
  l_DisplayData([],[],[1 2 3])
else
  l_DisplayData([],[],vtype)
end
%l_RefreshRoi(H,Dat)

catch
  aedes_errordump(lasterror);
end
end % function l_ChangeView(h,


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%
  %% Display axes units
  %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function l_ViewAxesUnits(h,evd,opt)
  try
  
  %tmp=get([H.UIVIEW_GRID,H.UIVIEW_FOV_UNITS],...
  %        'checked');
  tmp=get(H.UIVIEW_GRID,...
          'checked');
  
  %% Show pixel units
  if strcmpi(opt,'pixel')
    if strcmpi(tmp,'off')
      Dat.AxGap=30;
      set([H.IMAX1,H.IMAX2,H.IMAX3],...
          'visible','on',...
          'xgrid','on',...
          'ygrid','on',...
          'box','on',...
          'xtickmode','auto',...
          'ytickmode','auto',...
          'xticklabelmode','auto',...
          'yticklabelmode','auto',...
          'layer','top')
      set(H.UIVIEW_GRID,...
          'checked','on')
      set(H.UITOGGLE_GRID,...
          'state','on')
      Dat.ShowGrid = true;
      setpref('Aedes','ShowGrid',Dat.ShowGrid)
    else
      Dat.AxGap=4;
      set([H.IMAX1,H.IMAX2,H.IMAX3],...
          'visible','off')
      set(H.UIVIEW_GRID,...
          'checked','off')
      set(H.UITOGGLE_GRID,...
          'state','off')
      Dat.ShowGrid = false;
      setpref('Aedes','ShowGrid',Dat.ShowGrid)
    end
  
  %% Show FOV units
  elseif strcmpi(opt,'FOV')
    if isempty(Dat.FOV)
      return
    end
    
    if strcmpi(tmp{2},'off')
      xticklabel = (get(H.IMAX1,'xtick')/Dat.ImageDim(Dat.DataInd,2)*Dat.FOV(1)).';
      xticklabel = fix(xticklabel*100)/100;
      yticklabel = (get(H.IMAX1,'ytick')/Dat.ImageDim(Dat.DataInd,1)*Dat.FOV(2)).';
      yticklabel = fix(yticklabel*100)/100;
      
      Dat.AxGap=30;
      set(H.IMAX1,...
          'visible','on',...
          'xgrid','on',...
          'ygrid','on',...
          'box','on',...
          'xticklabel',xticklabel,...
          'yticklabel',yticklabel,...
          'layer','top')
      set(H.UIVIEW_GRID,...
          'checked','off')
      set(H.UIVIEW_FOV_UNITS,...
          'checked','on')
    else
      Dat.AxGap=4;
      set([H.IMAX1,H.IMAX2,H.IMAX3],...
          'visible','off')
      set(H.UIVIEW_GRID,...
          'checked','off')
      set(H.UIVIEW_FOV_UNITS,...
          'checked','off')
    end
    
  end
  
  l_AxesPositions([],[],Dat.AxView,Dat.ZoomLevel)
  
  catch
    aedes_errordump(lasterror);
  end
  end % function l_ViewAxesUnits(h,

%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Show fMRI time series
%%
%%%%%%%%%%%%%%%%%%%%%%%%%
function l_ShowTimeSeries(h,evd,opt)
try

% if data is mixed type, return
if Dat.isDataMixed
  
  isEPI = false;
  
  % Check that all images are the same size
  if all(Dat.ImageDim(:,1)./Dat.ImageDim(1,1)) && ...
    all(Dat.ImageDim(:,2)./Dat.ImageDim(1,2))
    UseDim = 5;
  else
    return
  end
else
  
  % Check if EPI or RASER data are used
  if ( ~Dat.isDataMixed && isfield(DATA{1},'PROCPAR') && ...
      isfield(DATA{1}.PROCPAR,'readres') ) || ...
      ( ~Dat.isDataMixed && isfield(DATA{1},'PROCPAR') && isfield(DATA{1}.PROCPAR,'teType') )
    isEPI = true;
  else
    isEPI = false;
  end
  
  if isfield(DATA{1},'PROCPAR') && isfield(DATA{1}.PROCPAR,'apptype') && strcmp(DATA{1}.PROCPAR.apptype,'imEPI')
    isEPI=true;
  end

  if isfield(Dat,'fMRIonsets') && ~isempty(Dat.fMRIonsets)
    isEPI=true;
  end
  
  % Determine time series dimension
  if Dat.ImageDim(4)>1
    % Assume the time series of a 4D Matrix to be in the 4th dimension
    UseDim = 4;
  elseif Dat.ImageDim(3)>1
    % Assume the time series of a 3D Matrix to be in the 3th dimension
    UseDim = 3;
  else
    return
  end
end

if ~isempty(opt) && strcmpi(opt,'toggle')
  ischecked=get(H.UIVIEW_TIMESERIES,'checked');
  if strcmpi(ischecked,'on')
    try
      close(H.TSFIG);
    end
    set(H.UIVIEW_TIMESERIES,'checked','off')
    Dat.ShowTimeSeries = false;
    return
  else
    set(H.UIVIEW_TIMESERIES,'checked','on')
    Dat.ShowTimeSeries = true;
  end
end


    

% Get voxel time series data
if UseDim == 4
  tsdata = squeeze(DATA{1}.FTDATA(Dat.Slices(1),Dat.Slices(2),Dat.Slices(3),:));
elseif UseDim == 3
  tsdata = squeeze(DATA{1}.FTDATA(Dat.Slices(1),Dat.Slices(2),:,1));
elseif UseDim == 5
  % For mixed data
  tsdata = [];
  for ii=1:length(DATA)
    tsdata(ii) = DATA{ii}.FTDATA(Dat.Slices(1),Dat.Slices(2));
  end
end

tsdata = double(tsdata);
tsdata=tsdata(:).';

% Dont show data from reference image and calculate persential change
if isEPI
  baseline = mean(tsdata(1:min(20,length(tsdata))));
  tsdata = (tsdata./baseline-1)*100;
  tt=aedes_trendest(tsdata,10);
  tt=tt(:).';
end

% normalize
%tsdata=tsdata./max(tsdata);

% Open time series figure, if not already open
if isfield(H,'TSFIG') && ~isempty(H.TSFIG) && ishandle(H.TSFIG)
  % Refresh line data
  set(H.TSLINE,'xdata',1:length(tsdata),...
    'ydata',tsdata);
  data_mx = max(tsdata);
  data_mn = min(tsdata);
  tmp = abs(diff([data_mn data_mx]))*0.05;
  ylim = [data_mn-tmp data_mx+tmp];
  if ylim(1)==ylim(2)
    ylim(1)=ylim(1)-1;
    ylim(2)=ylim(2)+1;
  end
  set(H.TSAX,'ylim',ylim)
  if isEPI
    set(H.TSTREND,'xdata',1:length(tsdata),...
      'ydata',tt);
    
    if isfield(Dat,'fMRIonsets') && isfield(Dat,'fMRIdurats') && ...
        ~isempty(Dat.fMRIonsets) && ~isempty(Dat.fMRIdurats)
      
      ydata = [ylim(1) ylim(1) ylim(2) ylim(2)];
      set(H.TSPATCH,'ydata',ydata)
    end
  end
else
  % Get default figure position
  try
    default_pos = getpref('Aedes','TSfig_position');
		scrsz = get(0,'Screensize');
		if default_pos(1)>scrsz(3) || default_pos(2)>scrsz(4)
			default_pos = get(0,'defaultfigureposition');
		end
  catch
    default_pos = get(0,'defaultfigureposition');
  end
  
  H.TSFIG = figure('NumberTitle','off',...
    'HandleVisibility','off',...
    'Name','Aedes: Voxel Time Series Plot',...
    'position',default_pos,...
    'closereq',@l_SaveTSfigPosition);
  
  % Suppress warning from get(fh,'javaFrame') generated in Matlab R2008a->
  if Dat.MatlabVersion>=7.06
    warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
  end
  
  % Try to set the figure as floating i.e. "always on top"
  jf = get(H.TSFIG,'javaFrame');
  pause(0.05) % Make Matlab more stable...
	if ~Dat.HG2graphics
		if Dat.MatlabVersion>=7.13
			wh = jf.fHG1Client.getWindow;
		else
			wh = jf.fFigureClient.getWindow;
		end
	else
		wh = handle(jf.fHG2Client.getWindow);
	end
  pause(0.05)
  wh.setAlwaysOnTop(true);
  pause(0.05)
  
  H.TSAX = axes('parent',H.TSFIG,...
                'xlim',[1 length(tsdata)]);%,'ylim',[0.7 1]);
  xlabel(H.TSAX,'Slice Number');
  if isEPI
    ylabel(H.TSAX,'BOLD (%)')
  else
    ylabel(H.TSAX,'Intensity (a.u.)')
  end
  
  title(H.TSAX,'')
  H.TSLINE = line('parent',H.TSAX,'xdata',1:length(tsdata),'ydata',tsdata);
  
  % Show trend and paradigm over EPI data
  if isEPI
    H.TSTREND = line('parent',H.TSAX,'xdata',1:length(tsdata),...
      'ydata',tt,'linewidth',2,'color','r');
    
    
    if isfield(Dat,'fMRIonsets') && isfield(Dat,'fMRIdurats') && ...
        ~isempty(Dat.fMRIonsets) && ~isempty(Dat.fMRIdurats)
      H.TSPATCH = [];
      for ii=1:length(Dat.fMRIonsets)
        xdata = [Dat.fMRIonsets(ii) Dat.fMRIonsets(ii)+Dat.fMRIdurats(ii) ...
          Dat.fMRIonsets(ii)+Dat.fMRIdurats(ii) Dat.fMRIonsets(ii)];
        tmp = get(H.TSAX,'ylim');
        ydata = [tmp(1) tmp(1) tmp(2) tmp(2)];
        H.TSPATCH(ii) = patch('parent',H.TSAX,...
          'xdata',xdata,'ydata',ydata,'FaceColor','r',...
          'FaceAlpha',0.3,'LineStyle','none');
      end
    end
  end
end



catch
  aedes_errordump(lasterror);
end

end % function l_ShowTimeSeries(h,

function l_SaveTSfigPosition(h,evd)

% Set current position to preferences
setpref('Aedes','TSfig_position',get(h,'position'))

% Delete figure
delete(h)

% Set checked status in uimenu
set(H.UIVIEW_TIMESERIES,'checked','off')
Dat.ShowTimeSeries = false;
H.TSFIG = [];

end % function l_SaveTSfigPostion(h,  

%%%%%%%%%%%%%%%%%%%%
%%
%% Quit Aedes
%%
%%%%%%%%%%%%%%%%%%%%
function l_quit(h,evd)
try
fh=H.FIG;

% Check if ROI(s) should be saved
cancel=l_CheckRoiSaved;
if cancel
  return
end

% Try to close file
if ~isempty(DATA)
  l_CloseFile([],[]);
end

% Clear public variables
clear ROI DATA Dat H

% Close GUI
delete(fh)


catch
  aedes_errordump(lasterror);
end
end % function l_quit(h,

%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Show/Hide Crossbars
%%
%%%%%%%%%%%%%%%%%%%%%%%
function l_ShowHideCrossbars(h,evd)
try
if  strcmpi(get(H.UITOGGLE_CROSSHAIRS,'state'),'on') || ...
		~isfield(H,'CROSSBAR_LN') || isempty(H.CROSSBAR_LN)
  % Draw crossbars
  ax_h = [H.IMAX1,H.IMAX2,H.IMAX3];
  if ~isempty(H.ROIAX)
    ax_h = H.ROIAX(:,end);
  elseif ~isempty(H.IMOVERLAY)
    ax_h = H.IMOVERLAYAX;
  end
  tmp_sl = [2 1;
           3 1;
           3 2];
  
  % Default visibility for lines is always off
    if strcmpi(get(H.UITOGGLE_CROSSHAIRS,'state'),'on')
      vis = 'on';
      %vis = 'off';
    else
      vis = 'off';
    end
  
  if ~isfield(H,'CROSSBAR_LN') || ~any(any(ishandle(H.CROSSBAR_LN)))
    for ii=1:3
      H.CROSSBAR_LN(ii,1) = line('parent',ax_h(ii),...
                                 'xdata',get(ax_h(ii),'xlim'),...
                                 'ydata',[Dat.Slices(tmp_sl(ii,2)) ...
                          Dat.Slices(tmp_sl(ii,2))],...
                                 'color','r',...
                                 'hittest','off',...
                                 'visible',vis);
      
      % Y-line
      H.CROSSBAR_LN(ii,2) = line('parent',ax_h(ii),...
                                 'xdata',[Dat.Slices(tmp_sl(ii,1)) ...
                          Dat.Slices(tmp_sl(ii,1))],...
                                 'ydata',get(ax_h(ii),'ylim'),...
                                 'color','r',...
                                 'hittest','off',...
                                 'visible',vis);
			if ~Dat.HG2graphics
				set([H.CROSSBAR_LN(ii,1),H.CROSSBAR_LN(ii,2)],...
					'erasemode','xor')
			else
				set([H.CROSSBAR_LN(ii,1),H.CROSSBAR_LN(ii,2)],'buttondownfcn',@l_SetMouseGestures,...
					'hittest','on')
				set([H.CROSSBAR_LN(ii,1),H.CROSSBAR_LN(ii,2)],'userdata',H.IM(ii))
				if exist('ROI','var') && length(ROI)>0
                    % keyboard;
                    % ROIAX is missing _if ROI loaded from command line
                    % because it is processed later. Let's add a check here
                    % and see what happens.. mjn 01/2021
                    if ~isempty(H.ROIAX)
					            set(H.CROSSBAR_LN(ii,:),'parent',H.ROIAX(ii,length(ROI)))
				            else
					            set(H.CROSSBAR_LN(ii,:),'parent',H.IMOVERLAYAX(ii))
				            end
        else
					set(H.CROSSBAR_LN(ii,:),'parent',H.IMOVERLAYAX(ii))
				end
			end
      
    end
  else
     for ii=1:3
      set(H.CROSSBAR_LN(ii,1),'xdata',get(ax_h(ii),'xlim'),...
                                 'ydata',[Dat.Slices(tmp_sl(ii,2)) ...
                          Dat.Slices(tmp_sl(ii,2))],...
                        'visible',vis,...
                        'parent',ax_h(ii));
      
      % Y-line
      set(H.CROSSBAR_LN(ii,2),'xdata',[Dat.Slices(tmp_sl(ii,1)) ...
                          Dat.Slices(tmp_sl(ii,1))],...
                                 'ydata',get(ax_h(ii),'ylim'),...
                        'visible',vis,...
                        'parent',ax_h(ii));
    end
  end

  % Set visibility for lines
  if  strcmpi(get(H.UITOGGLE_CROSSHAIRS,'state'),'on')
    if any(Dat.AxView==[1 2 3])
			ind = [3 2 1];
      set(H.CROSSBAR_LN(ind(Dat.AxView),:),'visible','on')
    else
      set(H.CROSSBAR_LN,'visible','on')
    end
  end
else
  try
    % Hide crossbars
    set(H.CROSSBAR_LN,'visible','off')
  end
end

catch
  aedes_errordump(lasterror);
end
end % function l_ShowHideCrossbars(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Update intensity value
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_UpdateIntensityValue(h,evd)
try
	
	% Get the current intensity value
	IntVal=DATA{Dat.DataInd}.FTDATA(Dat.Slices(1),...
		Dat.Slices(2),...
		Dat.Slices(3),...
		Dat.CurrentVol);
	
	% Update Voxel Value text
	set(H.VOXEL_VALUE,'string',...
		sprintf('%.2f',IntVal))
	
	% Look for a quick return
	if ~Dat.ShowColorbar
		return
  end

% Get the current intensity value
IntVal=DATA{Dat.DataInd}.FTDATA(Dat.Slices(1),...
                               Dat.Slices(2),...
                               Dat.Slices(3),...
                               Dat.CurrentVol);
            


% Make sure that it is double-class
IntVal=double(IntVal);

% The value for the current intensity identifier in the colorbar has to
% be between Clim(1) and Clim(2)
LnIntVal = min(max(Dat.Clim(1),IntVal),Dat.Clim(2));

%ydata = ones(1,2)*(((LnIntVal-Dat.Clim(1))/(Dat.Clim(2)-Dat.Clim(1)))*256);

lineHeight = (1.5/256)*diff(Dat.Clim);

if isfield(H,'COLORBAR_LN') && ~isempty(H.COLORBAR_LN)
  
  % Update colorbar-lin position
  set(H.COLORBAR_LN,...
      'position',[0.5 LnIntVal-(lineHeight/2) 3 lineHeight])
  
  % Update intensity value text
  set(H.COLORBAR_TX,...
      'position',[0 -3],...
      'string',sprintf('%.2f',IntVal))
else
  if Dat.ShowColorbar
    vis = 'on';
  else
    vis = 'off';
  end
  H.COLORBAR_LN = rectangle('parent',H.COLORBAR_AX,...
                            'position',[0.5 LnIntVal-(lineHeight/2) 3 lineHeight],...
                            'edgecolor','k',...
                            'facecolor','w',...
                            'visible',vis);
	if ~Dat.HG2graphics
		set(H.COLORBAR_LN,'erasemode','normal')
	end
  H.COLORBAR_TX = text('parent',H.COLORBAR_AX,...
                       'units','pixel',...
                       'position',[0 -3],...
                       'string',sprintf('%.2f',IntVal),...
                       'fontsize',8,...
                       'color','k',...
                       'fontweig','bold',...
                       'horizontalalign','left',...
                       'verticalalign','top',...
                       'visible',vis);
end


% Refresh axes
l_AxesPositions([],[],Dat.AxView,Dat.ZoomLevel);


catch
  aedes_errordump(lasterror);
end
end % function l_UpdateIntensityValue(h,

%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure Resize Function
%%%%%%%%%%%%%%%%%%%%%%%%%
function l_ResizeFcn(h,evd)
try
% Set sidebar height to the current figure height
%sidebar_pos = get(H.SIDEBAR,'position');
sidebar_pos = get(H.SIDEBAR_FRAME,'position');
fig_pos=get(H.FIG,'position');

% Don't let position be smaller than 300 by 200 px
if fig_pos(3)<250 | fig_pos(4)<150
  if fig_pos(3)<250
    fig_pos(3)=250;
    
    % Check that figure doesn't go off the screen
    scrsz = get(0,'screensize');
    if ( fig_pos(1)+fig_pos(3) ) > scrsz(3)
      fig_pos(1)=scrsz(3)-fig_pos(3);
    end
  end
  if fig_pos(4)<150
    fig_pos(4)=150;
    
    % Check that figure doesn't go off the screen
    scrsz = get(0,'screensize');
    if ( fig_pos(2)+fig_pos(4) ) > scrsz(4)
      fig_pos(2)=scrsz(4)-fig_pos(4)-50;
    end
  end
  set(H.FIG,'position',fig_pos)
end
%set(H.SIDEBAR,'position',[sidebar_pos(1:3) fig_pos(4)+1])
set(H.SIDEBAR_FRAME,'position',[0 0 232 fig_pos(4)+2])
%H.IMSLIDER_FRAME


% Get uicontrol handles
uich = findobj(H.FIG,'type','uicontrol');

% Exclude a few uicontrols
uich(find(uich==H.IMSLIDER_FRAME))=[];
uich(find(uich==H.SIDEBAR_FRAME))=[];
uich(find(uich==H.BLANK_FRAME))=[];
uich(find(uich==H.IMSLIDER(1)))=[];
uich(find(uich==H.IMSLIDER(2)))=[];

% Set uicontrol positions
uich_pos = get(uich,'position');
tmp=reshape([uich_pos{:}],4,length(uich_pos))';
tmp(:,2)=ones(size(tmp,1),1)*fig_pos(4)-...
  ones(size(tmp,1),1)*H.ORIG_FIG_POS(4)+...
  H.UICH_ORIG_POS(:,2);
set(H.BLANK_FRAME,'position',[0 0 1 1])
set(uich,{'position'},mat2cell(tmp,ones(1,length(tmp)),4))


% Set axis positions
if ~isempty(Dat) && isfield(Dat,'AxView') && isfield(Dat,'ZoomLevel')
  l_AxesPositions([],[],Dat.AxView,Dat.ZoomLevel);
end

catch
  aedes_errordump(lasterror);
end
end % function l_ResizeFcn(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% View image sliders if needed
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_SetImageSliders(h,evd,isXSlider,isYSlider,x_over,y_over)
try
fig_pos = get(H.FIG,'position');
if ~isXSlider & ~isYSlider
  %% Return immediately if slider are not needed
  set(H.IMSLIDER,'visible','off')
  set(H.IMSLIDER_FRAME,'visible','off')
  return
elseif isXSlider & ~isYSlider
  %% Show only X-slider
  if ~strcmpi(get(H.IMSLIDER(1),'visible'),'on') % Slider is not currently
    sl1_val = 0;                                 % visible
    set(H.IMSLIDER(1),'visible','on',...
      'value',sl1_val)
    set(H.IMSLIDER(2),'visible','off')
    set(H.IMSLIDER_FRAME,'visible','off')
    % Check if we should now need Y-Slider too
    %l_SetImageSliders([],[])
    %return
  else % Slider is already visible
    val = get(H.IMSLIDER(1),'value');
    s_max = get(H.IMSLIDER(1),'max');
    sl1_val = min((val/s_max)*abs(x_over),abs(x_over));
    set(H.IMSLIDER(1),'visible','on')
    set(H.IMSLIDER(2),'visible','off')
    set(H.IMSLIDER_FRAME,'visible','off')
  end
  tmp=get(H.IMSLIDER(1),'position');
  tmp_scale = min(20,abs(x_over));
% $$$   set(H.IMSLIDER(1),'Min',0,...
% $$$                     'Max',abs(x_over),...
% $$$                     'sliderstep',...
% $$$                     [tmp_scale/10/abs(x_over) ...
% $$$ 		    tmp_scale/abs(x_over)],...
% $$$                     'position',[tmp(1:2) fig_pos(3)-tmp(1)+1 17],...
% $$$                     'value',sl1_val)
  set(H.IMSLIDER(1),'Min',0,...
                    'Max',abs(x_over),...
                    'sliderstep',...
                    [tmp_scale/10/abs(x_over) ...
                     tmp_scale/abs(x_over)*20],...
                    'position',[tmp(1:2) fig_pos(3)-tmp(1)+1 17],...
                    'value',sl1_val)
  
elseif ~isXSlider & isYSlider
  %% Show only Y-Slider
  if ~strcmpi(get(H.IMSLIDER(2),'visible'),'on')
    sl2_val = 0;
    set(H.IMSLIDER(2),'visible','on',...
      'value',sl2_val)
    set(H.IMSLIDER(1),'visible','off')
    set(H.IMSLIDER_FRAME,'visible','off')
    % Check if we should now need X-Slider too
    %l_SetImageSliders([],[])
    %return
  else
    val = get(H.IMSLIDER(2),'value');
    s_min = get(H.IMSLIDER(2),'min');
    sl2_val = max(val/s_min*y_over,y_over);
    set(H.IMSLIDER(2),'visible','on')
    set(H.IMSLIDER(1),'visible','off')
    set(H.IMSLIDER_FRAME,'visible','off')
  end
  tmp=get(H.IMSLIDER(2),'position');
  tmp_scale = min(20,abs(y_over));
  set(H.IMSLIDER(2),'Min',y_over,...
                    'Max',0,...
                    'sliderstep',[tmp_scale/10/abs(y_over) ...
		    tmp_scale/abs(y_over)*20],...
                    'position',[fig_pos(3)-17 0 18 fig_pos(4)+1],...
                    'value',sl2_val)
else
  %% Show both sliders
  if ~strcmpi(get(H.IMSLIDER(1),'visible'),'on')
    sl1_val=0;
  else
    val = get(H.IMSLIDER(1),'value');
    s_max = get(H.IMSLIDER(1),'max');
    sl1_val = min((val/s_max)*abs(x_over),abs(x_over));
  end
  if ~strcmpi(get(H.IMSLIDER(2),'visible'),'on')
    sl2_val=0;
  else
    val = get(H.IMSLIDER(2),'value');
    s_min = get(H.IMSLIDER(2),'min');
    sl2_val = max(val/s_min*y_over,y_over);
  end
  tmp=get(H.IMSLIDER(1),'position');
  tmp_scale = min(20,abs(x_over));
  set(H.IMSLIDER(1),'Min',0,...
                    'Max',abs(x_over),...
                    'sliderstep',...
                    [tmp_scale/10/abs(x_over) ...
		    tmp_scale/abs(x_over)*20],...
                    'position',[tmp(1:2) fig_pos(3)-tmp(1)-17 17],...
                    'value',sl1_val)
  tmp_scale = min(20,abs(y_over));
  tmp=get(H.IMSLIDER(2),'position');
  set(H.IMSLIDER(2),'Min',y_over,...
                    'Max',0,...
                    'sliderstep',[tmp_scale/10/abs(y_over) ...
		    tmp_scale/abs(y_over)*20],...
                    'position',[fig_pos(3)-17 17 18 fig_pos(4)-17+1],...
                    'value',sl2_val)
  tmp=get(H.IMSLIDER_FRAME,'position');
  set(H.IMSLIDER_FRAME,'position',...
    [tmp(1:2) fig_pos(3)-tmp(1)+2 17])
  set(H.IMSLIDER,'visible','on')
  set(H.IMSLIDER_FRAME,'visible','on')
end

catch
  aedes_errordump(lasterror);
end
end % function l_SetImageSliders(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Drag zoomed view with mouse
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_DragView(h,evd,startpoint,val1,val2)
try
% Return if Image sliders are not currently visible...
vis1 = get(H.IMSLIDER(1),'visible');
vis2 = get(H.IMSLIDER(2),'visible');
if strcmpi(vis1,'off') && strcmpi(vis2,'off')
  return
end

% Get current figure point
cp = get(H.FIG,'currentpoint');
df=startpoint-cp;

if strcmpi(vis1,'on')
  val_max=get(H.IMSLIDER(1),'max');
  %val = get(H.IMSLIDER(1),'value');
  new_val = min(max(0,val1+df(1)),val_max);
  set(H.IMSLIDER(1),'value',new_val)
end
if strcmpi(vis2,'on')
  val_min=get(H.IMSLIDER(2),'min');
  val = get(H.IMSLIDER(2),'value');
  new_val = min(max(val_min,val2+df(2)),0);
  set(H.IMSLIDER(2),'value',new_val)
end
l_AxesPositions(H.IMSLIDER(1),[])

catch
  aedes_errordump(lasterror);
end
end % function l_DragView(h,

function l_PanView()

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Zoom View
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_Zoom(h,evd,input_factor)
try
% Zoom in
if ischar(input_factor) && strcmpi(input_factor,'+')
  
  % Return immediately if axes zoom level is
  % normalized
  if Dat.ZoomLevel==0
		% Calculate current absolute (normalized) zoom level
		if Dat.AxView==0 || Dat.AxView==3
			pos = get(H.IMAX1,'position');
			pos_cur = Dat.RealAxSize{Dat.DataInd}(1,1);
		elseif Dat.AxView==1
			pos = get(H.IMAX3,'position');
			pos_cur = Dat.RealAxSize{Dat.DataInd}(3,1);
		else
			pos = get(H.IMAX2,'position');
			pos_cur = Dat.RealAxSize{Dat.DataInd}(2,1);
		end
		pos = pos(3);
		currentScale = pos/pos_cur;
		zoom_level = round(currentScale*10)/10;
  else
		% Get current zoom level
		zoom_level=Dat.ZoomLevel;
	end

	if zoom_level>=20
		return
	end
	set(H.UITOGGLE_ZOOMNORM,'state','off')
	
  % Set new zoom_level
  if zoom_level==0.5
    zoom_level=1;
  else
    zoom_level=zoom_level+Dat.ZoomStep;
    if zoom_level>=20
      set(H.UIPUSH_ZOOMIN,'enable','off')
    end
  end
  set(H.UIPUSH_ZOOMOUT,'enable','on')
  
  % Zoom out
elseif ischar(input_factor) && strcmpi(input_factor,'-')
  
  % Return immediately if axes zoom level is
  % normalized
  if Dat.ZoomLevel==0
    % Calculate current absolute (normalized) zoom level
		if Dat.AxView==0 || Dat.AxView==3
			pos = get(H.IMAX1,'position');
			pos_cur = Dat.RealAxSize{Dat.DataInd}(1,1);
		elseif Dat.AxView==1
			pos = get(H.IMAX3,'position');
			pos_cur = Dat.RealAxSize{Dat.DataInd}(3,1);
		else
			pos = get(H.IMAX2,'position');
			pos_cur = Dat.RealAxSize{Dat.DataInd}(2,1);
		end
		pos = pos(3);
		currentScale = pos/pos_cur;
		zoom_level = round(currentScale*10)/10;
	else
		 % Get current zoom level
		 zoom_level=Dat.ZoomLevel;
	end
 
  if zoom_level==0.5
    return
  end
	set(H.UITOGGLE_ZOOMNORM,'state','off')

  % Set new zoom_level
  if zoom_level-Dat.ZoomStep<0.5
    zoom_level=0.5;
    set(H.UIPUSH_ZOOMOUT,'enable','off')
  else
    zoom_level=zoom_level-Dat.ZoomStep;
  end
  set(H.UIPUSH_ZOOMIN,'enable','on')
  %Dat.AxSize=round(zoom_level*Dat.RealAxSize);
  
elseif ischar(input_factor) && strcmpi(input_factor,'normalize')
  
  if Dat.ZoomLevel==0
    if isempty(Dat.OldZoom)
      set([H.UIPUSH_ZOOMIN,H.UIPUSH_ZOOMOUT],'enable','on')
      l_Zoom([],[],1)
    else
      set([H.UIPUSH_ZOOMIN,H.UIPUSH_ZOOMOUT],'enable','on')
      l_Zoom([],[],Dat.OldZoom)
      if isfield(Dat,'OldImSliderVals') && ~isempty(Dat.OldImSliderVals)
        h1max = get(H.IMSLIDER(1),'max');
        h1min = get(H.IMSLIDER(1),'min');
        h2max = get(H.IMSLIDER(2),'max');
        h2min = get(H.IMSLIDER(2),'min');
        set(H.IMSLIDER(1),'value',max(min(Dat.OldImSliderVals(1),h1max),h1min))
        set(H.IMSLIDER(2),'value',max(min(Dat.OldImSliderVals(2),h2max),h2min))
        l_AxesPositions([],[],Dat.AxView,Dat.ZoomLevel);
      end
    end
  else
    l_Zoom([],[],0)
  end
  return
  
elseif input_factor==0
  
  % Get previous zoom level
  %ind=find(strcmpi(get(H.UIZOOM,'checked'),'on'));
  zoom_level=Dat.ZoomLevel;%get(H.UIZOOM(ind),'userdata');
  if zoom_level~=0
    Dat.OldZoom = zoom_level;
    % Save image slider values
    if strcmpi(get(H.IMSLIDER(1),'visible'),'on') || ...
          strcmpi(get(H.IMSLIDER(2),'visible'),'on')
      Dat.OldImSliderVals(1) = get(H.IMSLIDER(1),'value');
      Dat.OldImSliderVals(2) = get(H.IMSLIDER(2),'value');
    else
      Dat.OldImSliderVals=[];
    end
  end
  %Dat.AxSize=0;
  zoom_level=0;
  set([H.UIPUSH_ZOOMIN,H.UIPUSH_ZOOMOUT],'enable','on')
  set(H.UITOGGLE_ZOOMNORM,'state','on')
  
else
  zoom_level=input_factor;
  %Dat.AxSize=round(zoom_level*Dat.RealAxSize);
  if zoom_level==20
    set(H.UIPUSH_ZOOMOUT,'enable','on')
    set(H.UIPUSH_ZOOMIN,'enable','off')
  elseif zoom_level==0.5
    set(H.UIPUSH_ZOOMIN,'enable','on')
    set(H.UIPUSH_ZOOMOUT,'enable','off')
  else
    set(H.UIPUSH_ZOOMIN,'enable','on')
    set(H.UIPUSH_ZOOMOUT,'enable','on')
  end
  set(H.UITOGGLE_ZOOMNORM,'state','off')
end

% Set uimenu checked status
setpref('Aedes','ZoomLevel',zoom_level)
Dat.ZoomLevel=zoom_level;

% Set axis positions
l_AxesPositions([],[],Dat.AxView,Dat.ZoomLevel);

% Update info text
l_UpdateInfoText([],[])

catch
  aedes_errordump(lasterror);
end
end % function l_Zoom(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Rotate Images
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_RotateFlip(h,evd,input_factor)
try

% Don't allow rotation when an overlay is loaded
if Dat.isImageOverlayLoaded
  h = errordlg(['Rotating or flipping is not allowed when an',...
	' overlay is loaded. Sorry...'],'Rotate/flip not allowed!',...
	'modal');
  return
end

% Check if all images are square
squareImages = all((Dat.ImageDim(:,1)./Dat.ImageDim(:,2))==1);

% Reset rotation
if strcmpi(input_factor,'reset')
  
  % Look for a quick return
  if all(Dat.DataRotation==0) && all(Dat.DataFlip==0)
    % Nothing rotated or flipped -> reset not needed
    if ~isempty(h) && h==H.UIROTRESET
      hh=warndlg({'Images are already in their original orientations.',...
                  'Nothing to do...'},'Nothing to do...','modal');
    end
    return
  end
  
  % Ask for confirmation if called from uimenu
  if ~isempty(h) && h==H.UIROTRESET
    resp=questdlg({'This will reset ALL images to their original orientations?',...
                  'Are you sure?'},...
                  'Reset image orientations?','Yes','No','No');
    if isempty(resp) || strcmpi(resp,'No')
      return
    end
  end
  
  % Loop over images
  for ii=1:Dat.DataL
    if Dat.DataRotation(ii)~=0
      rot=abs(Dat.DataRotation(ii)-4);
      isRotEven = aedes_iseven(rot);
      
      if ~isRotEven && Dat.DataFlip(ii)~=0
        if Dat.DataFlip(ii)==1
          Dat.DataFlip(ii)=2;
        elseif Dat.DataFlip(ii)==2
          Dat.DataFlip(ii)=1;
        end
      end
      
      % Rotate image
      DATA{ii}.FTDATA = rot90(DATA{ii}.FTDATA,rot);
      % Rotate ROI(s)
      if ~isempty(ROI)
        for kk=1:length(ROI)
          ROI(kk).voxels{ii} = rot90(ROI(kk).voxels{ii},rot);
        end
      end
      Dat.DataRotation(ii)=0;
    end
    
    if Dat.DataFlip(ii)~=0
      % Flip image
      if Dat.DataFlip(ii)==1
        DATA{ii}.FTDATA = flipud(DATA{ii}.FTDATA);
      elseif Dat.DataFlip(ii)==2
        DATA{ii}.FTDATA = fliplr(DATA{ii}.FTDATA);
      end
      
      % Flip ROI(s)
      if ~isempty(ROI)
        for kk=1:length(ROI)
          if Dat.DataFlip(ii)==1
            ROI(kk).voxels{ii} = flipud(ROI(kk).voxels{ii});
          elseif Dat.DataFlip(ii)==2
            ROI(kk).voxels{ii} = fliplr(ROI(kk).voxels{ii});
          end
        end
      end
      Dat.DataFlip(ii)=0;
    end
  end
  
  % Update view if called from uimenu
  if ~isempty(h) && h==H.UIROTRESET
    l_DisplayData([],[])
  end
  
elseif strcmpi(input_factor,'custom')

  % Rotate/Flip Custom
  Out = aedes_rotateflip(DATA,squareImages);
  if isempty(Out)
    return
  end
  
  % Reset images before rotating/flipping
  %l_RotateFlip([],[],'reset')
  
  % Rotate and Flip images
  for ii=1:Dat.DataL
    % Rotate
    if Out.Rotate(ii)~=0
      DATA{ii}.FTDATA = rot90(DATA{ii}.FTDATA,Out.Rotate(ii));
      if ~isempty(ROI)
        % Rotate ROI(s)
        for kk=1:length(ROI)
          ROI(kk).voxels{ii} = rot90(ROI(kk).voxels{ii},Out.Rotate(ii));
        end
      end
      Dat.DataRotation(ii)=Dat.DataRotation(ii)+Out.Rotate(ii);
      if Dat.DataRotation(ii)>=4
        Dat.DataRotation(ii)=Dat.DataRotation(ii)-4;
      end
      
      if Out.Rotate~=2
        if Dat.DataFlip(ii)==1
          Dat.DataFlip(ii)=2;
        elseif Dat.DataFlip(ii)==2;
          Dat.DataFlip(ii)=1;
        end
      end
    end
    
    % Flip
    if Out.Flip(ii)==1
      DATA{ii}.FTDATA = flipud(DATA{ii}.FTDATA);
      if ~isempty(ROI)
        % Flip ROI(s)
        for kk=1:length(ROI)
          ROI(kk).voxels{ii} = flipud(ROI(kk).voxels{ii});
        end
      end
      if Dat.DataFlip(ii)==0
        Dat.DataFlip(ii)=1;
      elseif Dat.DataFlip(ii)==1
        Dat.DataFlip(ii)=0;
      elseif Dat.DataFlip(ii)==2
        Dat.DataFlip(ii)=0;
        Dat.DataRotation(ii)=Dat.DataRotation(ii)+2;
        if Dat.DataRotation(ii)>=4
          Dat.DataRotation(ii)=Dat.DataRotation(ii)-4;
        end
      end
      
    elseif Out.Flip(ii)==2
      DATA{ii}.FTDATA = fliplr(DATA{ii}.FTDATA);
      if ~isempty(ROI)
        % Flip ROI(s)
        for kk=1:length(ROI)
          ROI(kk).voxels{ii} = fliplr(ROI(kk).voxels{ii});
        end
      end
      if Dat.DataFlip(ii)==0
        Dat.DataFlip(ii)=2;
      elseif Dat.DataFlip(ii)==1
        Dat.DataFlip(ii)=0;
        Dat.DataRotation(ii)=Dat.DataRotation(ii)+2;
        if Dat.DataRotation(ii)>=4
          Dat.DataRotation(ii)=Dat.DataRotation(ii)-4;
        end
      elseif Dat.DataFlip(ii)==2
        Dat.DataFlip(ii)=0;
      end
      
    end
  end
  
  % Update view
  l_DisplayData([],[])
  
elseif any(strcmpi(input_factor,{'fliplr','flipud'}))
  
  if strcmpi(input_factor,'FlipLR')
    % FlipLR Current Slice
    DATA{Dat.DataInd}.FTDATA = fliplr(DATA{Dat.DataInd}.FTDATA);
    if ~isempty(ROI)
      % Flip ROI(s)
      for kk=1:length(ROI)
        ROI(kk).voxels{Dat.DataInd} = fliplr(ROI(kk).voxels{Dat.DataInd});
      end
    end
    if Dat.DataFlip(Dat.DataInd)==0
      Dat.DataFlip(Dat.DataInd)=2;
    elseif Dat.DataFlip(Dat.DataInd)==1
      Dat.DataFlip(Dat.DataInd)=0;
      Dat.DataRotation(Dat.DataInd)=Dat.DataRotation(Dat.DataInd)+2;
      if Dat.DataRotation(Dat.DataInd)>=4
        Dat.DataRotation(Dat.DataInd)=Dat.DataRotation(Dat.DataInd)-4;
      end
    elseif Dat.DataFlip(Dat.DataInd)==2
      Dat.DataFlip(Dat.DataInd)=0;
    end
    
        
  elseif strcmpi(input_factor,'FlipUD')
    % FlipUD Current Slice
    DATA{Dat.DataInd}.FTDATA = flipud(DATA{Dat.DataInd}.FTDATA);
    if ~isempty(ROI)
      % Flip ROI(s)
      for kk=1:length(ROI)
        ROI(kk).voxels{Dat.DataInd} = flipud(ROI(kk).voxels{Dat.DataInd});
      end
    end
    if Dat.DataFlip(Dat.DataInd)==0
      Dat.DataFlip(Dat.DataInd)=1;
    elseif Dat.DataFlip(Dat.DataInd)==1
      Dat.DataFlip(Dat.DataInd)=0;
    elseif Dat.DataFlip(Dat.DataInd)==2
      Dat.DataFlip(Dat.DataInd)=0;
      Dat.DataRotation(Dat.DataInd)=Dat.DataRotation(Dat.DataInd)+2;
      if Dat.DataRotation(Dat.DataInd)>=4
        Dat.DataRotation(Dat.DataInd)=Dat.DataRotation(Dat.DataInd)-4;
      end
    end
        
  end
  
  % Update view
  l_DisplayData([],[])
  
elseif strncmpi(input_factor,'3d',2)
  % Rotate volume data ------------------------------
  
  
  if strncmpi(input_factor,'3dx',3)
		dim = 1;
	elseif strncmpi(input_factor,'3dy',3)
		dim = 2;
	elseif strncmpi(input_factor,'3dz',3)
		dim = 3;
	end
	if strcmpi(input_factor(5:end),'90')
		k = 1;
	elseif strcmpi(input_factor(5:end),'180')
		k = 2;
	elseif strcmpi(input_factor(5:end),'270')
		k = 3;
	end
  
  % Rotate DATA
  DATA{Dat.DataInd}.FTDATA = ...
	aedes_rot3d(DATA{Dat.DataInd}.FTDATA,k,dim);
  
  % Rotate ROI
  if ~isempty(ROI)
		for ii=1:length(ROI)
			ROI(ii).voxels{Dat.DataInd} = ...
				aedes_rot3d(ROI(ii).voxels{Dat.DataInd},k,dim);
		end
	end
	
  clim = Dat.Clim;
	RotateFlip3d = Dat.RotateFlip3d;
  
  % Close the current file but preserve DATA and ROI structures
  l_CloseFile([],[],'PreserveData')
  
  % Re-initialize DATA
  l_Initialize([],[])
  
  % Load ROI(s)
  if ~isempty(ROI)
	l_RoiLoad([],[],'')
  end
  
  % Set rotation and flipping information
  Dat.RotateFlip3d{end+1} = {'rotate',k,dim};
  
  % Set Clim value
  Dat.Clim = clim;
  l_ChangeClim(h,evd)
  
elseif any(strcmpi(input_factor,{'flipx','flipy','flipz','flipv'}))
  % Flip volume data --------------------------------
  
  switch lower(input_factor)
		case 'flipx'
			dim = 1;
		case 'flipy'
			dim = 2;
		case 'flipz'
			dim = 3;
		case 'flipv'
			dim = 4;
	otherwise
	  % You shouldn't ever get here...
		return
  end
  
  % Flip DATA
  DATA{Dat.DataInd}.FTDATA = ...
		flipdim(DATA{Dat.DataInd}.FTDATA,dim);
  
  % Flip ROIs
  if ~isempty(ROI)
		for ii=1:length(ROI)
			ROI(ii).voxels{Dat.DataInd} = ...
				flipdim(ROI(ii).voxels{Dat.DataInd},dim);
		end
  end
  
  clim = Dat.Clim;
  RotateFlip3d = Dat.RotateFlip3d;
  
  % Close the current file but preserve DATA and ROI structures
  l_CloseFile([],[],'PreserveData')
  
  % Re-initialize DATA
  l_Initialize([],[])
  
  % Load ROI(s)
  if ~isempty(ROI)
	l_RoiLoad([],[],'')
  end
  
  % Set rotation and flipping information
  Dat.RotateFlip3d{end+1} = {'flip',dim};
  
  % Set Clim value
  Dat.Clim = clim;
  l_ChangeClim(h,evd)
  
else
  
  % Rotate Current slice
  if strcmpi(input_factor,'90')
    rotcount = 1;
  elseif strcmpi(input_factor,'180')
    rotcount = 2;
  else
    rotcount = 3;
  end
  
  if any(rotcount==[1 3]) && ~squareImages
    h=warndlg({'Rotation operation not allowed!','',...
               ['The 90 and 270 degree rotations are only permitted when',...
                ' all images in the image stack are square.']},...
              'Rotation not allowed!','modal');
    return
  end
  
  DATA{Dat.DataInd}.FTDATA = rot90(DATA{Dat.DataInd}.FTDATA,rotcount);
  if ~isempty(ROI)
    % Rotate ROI(s)
    for kk=1:length(ROI)
      ROI(kk).voxels{Dat.DataInd} = rot90(ROI(kk).voxels{Dat.DataInd},rotcount);
    end
  end
  Dat.DataRotation(Dat.DataInd)=Dat.DataRotation(Dat.DataInd)+rotcount;
  if Dat.DataRotation(Dat.DataInd)>=4
    Dat.DataRotation(Dat.DataInd)=Dat.DataRotation(Dat.DataInd)-4;
  end
  
  % Set flipping
  if rotcount~=2
    if Dat.DataFlip(Dat.DataInd)==1
      Dat.DataFlip(Dat.DataInd)=2;
    elseif Dat.DataFlip(Dat.DataInd)==2
      Dat.DataFlip(Dat.DataInd)=1;
    end
  end
  
  % Update view
  l_DisplayData([],[])
  
end


catch
  aedes_errordump(lasterror);
end
end % function l_Rotate(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Change Colorbar limits
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_ChangeColorbarLimits(h,evd)
try

set(H.COLORBAR_AX,'ylim',Dat.Clim)
set(H.COLORBAR_IM,'YData',Dat.Clim)

% Refresh the intensity value bar
l_UpdateIntensityValue([],[])

% $$$ if Dat.Clim(2)>1 || floor(diff(Dat.Clim))<1
% $$$   yticklabel = ceil(Dat.Clim(1)):floor(diff(Dat.Clim)*0.1):floor(Dat.Clim(2));
% $$$   %yticklabel = round(linspace(Dat.Clim(1),Dat.Clim(2),11));
% $$$ else
% $$$   yticklabel = linspace(Dat.Clim(1),Dat.Clim(2),11);
% $$$ end
% $$$ yticklabel = unique(yticklabel);
% $$$ ytick = ((yticklabel-Dat.Clim(1))./diff(Dat.Clim))*256;
% $$$ 
% $$$ %ytick = yticklabel-yticklabel(1);
% $$$ %ytick=(ytick./ytick(end))*256;
% $$$ %ytick
% $$$ set(H.COLORBAR_AX,'ytick',ytick)
% $$$ %ytick=get(H.COLORBAR_AX,'ytick');
% $$$ 
% $$$ set(H.COLORBAR_AX,'yticklabel',yticklabel,...
% $$$                   'yminortick','on')
% $$$ %set(H.COLORBAR_AX,'yticklabel',Dat.Clim(1)+(diff(Dat.Clim)./256).*ytick,...
% $$$ %                  'yminortick','on')

catch
  aedes_errordump(lasterror);
end
end % function l_ChangeColorbarLimits(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Modify Clim directly
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function l_ChangeClim(h,evd)
  try
  
  %% Called from Clim Min
  if h==H.CLIMMINEDIT
    val = str2num(get(H.CLIMMINEDIT,'string'));
		if isempty(val) && strcmpi(get(H.CLIMMINEDIT,'string'),'min')
			val = Dat.OrigClim(Dat.DataInd,1);
			set(H.CLIMMINEDIT,'string',num2str(val));
		end
    if isempty(val) || ~isreal(val) || val>=Dat.Clim(2)
      set(H.CLIMMINEDIT,'string',num2str(get(H.CLIMMINEDIT,'userdata')))
      return
    end
    Dat.Clim(1) = val;
    set(H.CLIMMINEDIT,'userdata',Dat.Clim(1))
  elseif h==H.CLIMMAXEDIT
    val = str2num(get(H.CLIMMAXEDIT,'string'));
		if isempty(val) && strcmpi(get(H.CLIMMAXEDIT,'string'),'max')
			val = Dat.OrigClim(Dat.DataInd,2);
			set(H.CLIMMAXEDIT,'string',num2str(val));
		end
    if isempty(val) || ~isreal(val) || val<=Dat.Clim(1)
      set(H.CLIMMAXEDIT,'string',num2str(get(H.CLIMMAXEDIT,'userdata')))
      return
    end
    Dat.Clim(2) = val;
    set(H.CLIMMAXEDIT,'userdata',Dat.Clim(2))
  else
		set(H.CLIMMINEDIT,'userdata',Dat.Clim(1))
		set(H.CLIMMAXEDIT,'userdata',Dat.Clim(2))
		set(H.CLIMMINEDIT,'string',num2str(get(H.CLIMMINEDIT,'userdata')))
		set(H.CLIMMAXEDIT,'string',num2str(get(H.CLIMMAXEDIT,'userdata')))
  end
  
  % Refresh axes with new Clim value
  set([H.IMAX1,H.IMAX2,H.IMAX3],'Clim',Dat.Clim)
  
%   % Save Clim values and change the range to slice locked position
% 	if Dat.isDataMixed
% 		if Dat.LockClim==2
% 			Dat.SliceClim(Dat.DataInd,:)=Dat.Clim;
% 		elseif Dat.LockClim==0 || Dat.LockClim==3
% 			Dat.LockClim = 1;
% 			set(H.CLIMRANGE_POPUP,'value',1)
% 		end
% 	end
	
	% Save Clim values and change the range to slice locked position
	if Dat.isDataMixed
		if any(Dat.LockClim==[0,2,3])
			Dat.LockClim = 2;
			Dat.SliceClim(Dat.DataInd,:)=Dat.Clim;
			set(H.CLIMRANGE_POPUP,'value',2)
		end
	end

  % Refresh colorbar limits
  l_ChangeColorbarLimits([],[])

	% Refresh window/level sliders
	l_SetWindowLevel([],0,Dat.Clim,0)
  
  catch
    aedes_errordump(lasterror);
  end
  end % function l_ChangeClim(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Change Clim Range
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_ChangeClimRange(h,evd)
try
  
val=get(H.CLIMRANGE_POPUP,'value');
if val==1
  Dat.LockClim = 1; % Global locked
elseif val==2
  Dat.LockClim = 2; % Slice locked
  l_RefreshClim([],[])
elseif val==3
  Dat.LockClim = 0; % Unlocked min-max
  l_RefreshClim([],[])
elseif val==4
  Dat.LockClim = 3; % Unlocked with auto-balance
  l_RefreshClim([],[])
end
%setpref('Aedes','LockClim',Dat.LockClim);

catch
  aedes_errordump(lasterror);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Refresh Clim values (locked or unlocked)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_RefreshClim(h,evd)
try

if Dat.LockClim==0
  %% Clim unlocked
  Dat.Clim = Dat.OrigClim(Dat.DataInd,:);
  l_SetWindowLevel(0,[],Dat.Clim,0)
  
elseif Dat.LockClim==1
  %% Clim locked global
	l_SetWindowLevel(0,[],Dat.Clim,0)
  

elseif Dat.LockClim==2
  %% Clim locked slice
  Dat.Clim = Dat.SliceClim(Dat.DataInd,:);
  %set([H.IMAX1,H.IMAX2,H.IMAX3],'Clim',Dat.Clim)
  %set(H.CLIMMINEDIT,'userdata',Dat.Clim(1),...
  %                  'string',num2str(Dat.Clim(1)))
  %set(H.CLIMMAXEDIT,'userdata',Dat.Clim(2),...
  %                  'string',num2str(Dat.Clim(2)))
  
  % Refresh colorbar limits
  %l_ChangeColorbarLimits([],[])

	l_SetWindowLevel(0,[],Dat.Clim,0)
  
elseif Dat.LockClim==3
  %% Clim unlocked with auto-balance
  l_SetAutoWindowLevel([],[])
end

catch
  aedes_errordump(lasterror);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Modify Window/Level
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_SetWindowLevel(h,evd,clim_in,wbm)
try
if nargin==4
  if ~isempty(clim_in) & wbm==0
		Dat.Clim = clim_in;

		% Set Window and Level according to range
		[Dat.Window,Dat.Level]=l_RangeToWindowLevel(clim_in);
    
    % Set slider values
    set(H.WINDOW_SLIDER,'value',Dat.Window)
    set(H.LEVEL_SLIDER,'value',Dat.Level)
  else % The function is called from WindowButtonMotion
    
    % Get figure's current point
    cp=get(H.FIG,'currentpoint');
    
    % Get temporary point and save current point
    cp_old = Dat.tmp_point;
    Dat.tmp_point = cp;
    
    % Compare current point to the initial point
    xy=cp_old-cp;
    
    % Use "tol" pixels to represent 1% shift in window
    % and level...
    tolx = 1; % Set 1% tolerance
    toly = 1;
    
    shift_window=xy(1)/tolx;
    shift_level=xy(2)/toly;
    Dat.Window = Dat.Window-shift_window;
    Dat.Level = Dat.Level+shift_level;
    Dat.Window = floor(Dat.Window);
    Dat.Level = floor(Dat.Level);
    
    % Make sure that values are within limits
    Dat.Window = max(min(Dat.Window,100),0);
    Dat.Level = max(min(Dat.Level,100),0);
		Dat.Clim = l_WindowLevelToRange([Dat.Window,Dat.Level]);

    % Set slider values
    set(H.WINDOW_SLIDER,'value',Dat.Window)
    set(H.LEVEL_SLIDER,'value',Dat.Level)

		if Dat.isDataMixed && any(Dat.LockClim==[0,2,3])
			Dat.LockClim = 2;
			Dat.SliceClim(Dat.DataInd,:)=Dat.Clim;
			set(H.CLIMRANGE_POPUP,'value',2)
		end
  end
else
  Clim = Dat.Clim;
end

% If called from Window Slider ---------------
if h==H.WINDOW_SLIDER
  val = get(h,'value');

	Dat.Window = floor(val);
  
% If called from Level Slider ----------------
elseif h==H.LEVEL_SLIDER
  val = get(h,'value');

  Dat.Level = floor(val);

% If called from Window Editbox -------------------
elseif h==H.WINDOW_EDIT
  
  val=str2num(get(h,'string'));
  % Check input values
  if ~isreal(val) || isempty(val) || val<0 || val>100
    set(h,'string',num2str(get(h,'userdata')))
  else
		Dat.Window = floor(val);
  end
  set(H.WINDOW_SLIDER,'value',Dat.Window)
  
% If called from Level Editbox ----------------
elseif h==H.LEVEL_EDIT
  
  val=str2num(get(h,'string'));
  % Check input values
  if ~isreal(val) || isempty(val) || val<-100 || val>100
    set(h,'string',num2str(get(h,'userdata')))
  else
		Dat.Level = floor(val);
  end
  set(H.LEVEL_SLIDER,'value',Dat.Level)
  
end

% Update Edit box values
set(H.WINDOW_EDIT,'string',sprintf('%.0f',Dat.Window))
set(H.LEVEL_EDIT,'string',sprintf('%.0f',Dat.Level))

% If called from Clim edit, do not update clim
if evd==0
	return
end

% Save Clim values and change the range to slice locked position
if any(h==[H.WINDOW_SLIDER,H.LEVEL_SLIDER,...
		H.WINDOW_EDIT,H.LEVEL_EDIT])
	Dat.Clim = l_WindowLevelToRange([Dat.Window,Dat.Level]);
	if Dat.isDataMixed
		if any(Dat.LockClim==[0,2,3])
			Dat.LockClim = 2;
			Dat.SliceClim(Dat.DataInd,:)=Dat.Clim;
			set(H.CLIMRANGE_POPUP,'value',2)
		end
	end
end

% Verify that Clim is valid
if diff(Dat.Clim)<=0
  Dat.Clim(2) = Dat.Clim(2)+abs(diff(Dat.Clim))+1;
end
  
% Refresh axes with new Clim value
set([H.IMAX1,H.IMAX2,H.IMAX3],'Clim',Dat.Clim)


% Update Clim editboxes
set(H.CLIMMINEDIT,'string',num2str(Dat.Clim(1)),...
                  'userdata',Dat.Clim(1))
set(H.CLIMMAXEDIT,'string',num2str(Dat.Clim(2)),...
                  'userdata',Dat.Clim(2))



% Refresh colorbar limits
l_ChangeColorbarLimits([],[])

catch
  aedes_errordump(lasterror);
end
end % function l_SetWindowLevel(h,

% Map Window/Level to Range ----------------------
function range = l_WindowLevelToRange(WL_in)

	% Make sure that Window is 0..100
	window = min(max(WL_in(1),0),100);
	window_w = diff(Dat.OrigClim(Dat.DataInd,:))*(window/100);

	% Make sure that level is 0..100
	level = min(max(WL_in(2),0),100);
	level_abs = (level)/100*diff(Dat.OrigClim(Dat.DataInd,:))+Dat.OrigClim(Dat.DataInd,1);

	% Calculate range
	range = [level_abs-window_w/2 level_abs+window_w/2];

end

% Map Range to Window/Level ----------------------
function [window,level] = l_RangeToWindowLevel(range_in)
	
	% Calculate window
	window = diff(range_in)/diff(Dat.OrigClim(Dat.DataInd,:))*100;
	window = min(max(window,0),100);
	
	% Calculate level
	level_abs = range_in(1)+diff(range_in)/2;
	level = (level_abs-Dat.OrigClim(Dat.DataInd,1))/diff(Dat.OrigClim(Dat.DataInd,:))*100;
	level = min(max(level,0),100);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Set window/level automatically
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_SetAutoWindowLevel(h,evd)
try
Clim = zeros(1,2);

tmp_data = DATA{Dat.DataInd}.FTDATA(:,:,:,Dat.CurrentVol);
max_val = max(tmp_data(:));
  
if any(strcmpi(class(DATA{Dat.DataInd}.FTDATA),{'single','double'}))
  % Saturate 1% of the low and high values
  normClim = stretchlim(tmp_data./max_val,...
    [0.005 0.995]);
  normClim = max(normClim,[],2);
  if max_val<1
    Clim(1) = max_val*normClim(1);
    Clim(2) = max_val*normClim(2);
  else
    Clim(1) = fix(max_val*normClim(1)*100)/100;
    Clim(2) = fix(max_val*normClim(2)*100)/100;
  end
else
	% Saturate 1% of the low and high values
	high_cut = 1-0.005;
	low_cut = 0.005;
	tmp = sort(tmp_data(:));
	n = length(tmp);
	high_ind = floor(n*high_cut);
	low_ind = ceil(n*low_cut);
	Clim = double([tmp(low_ind) tmp(high_ind)]);
end

% Update Window/Level
l_SetWindowLevel(0,[],Clim,0)

catch
  aedes_errordump(lasterror);
end
end % function l_SetAutoWindowLevel(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Build the recent files menu
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_BuildRecentFilesMenu(h,evd,opt,added_file)

if nargin>2 && ~isempty(opt)
  FILEMENU_OPEN_RECENT = opt;
else
  FILEMENU_OPEN_RECENT = H.FILEMENU_OPEN_RECENT;
end

% Clear all children from the current menu
delete(get(FILEMENU_OPEN_RECENT,'children'))

% Build recent files menu ----------------------------
try
  Dat.RecentFiles = getpref('Aedes','RecentFiles');
catch
  Dat.RecentFiles = {};
end

% Add a new file in the menu
if nargin==4
  
  % Check if the file already exists in the menu
  ind = strcmpi(Dat.RecentFiles,added_file);
  if any(ind)
	Dat.RecentFiles(find(ind)) = [];
	Dat.RecentFiles{end+1} = added_file;
  else
	Dat.RecentFiles{end+1} = added_file;
  end
  if length(Dat.RecentFiles)>9
	Dat.RecentFiles = {Dat.RecentFiles{end-9+1:end}};
  end
  setpref('Aedes','RecentFiles',Dat.RecentFiles)
end

RecentFiles = fliplr(Dat.RecentFiles);
for ii=1:length(RecentFiles)

  % Limit the length of the path
  [fp,fn,fe] = fileparts(RecentFiles{ii});
  filename = [fn,fe];
  if length(RecentFiles{ii})>70
	if length(filename)>60
	  if ispc
		f_label = [fp(1:3),filesep,'...',filesep,filename];
	  else
		f_label = [filesep,'...',filesep,filename];
	  end
	else
	  tmp=fliplr(RecentFiles{ii});
	  f_label = [fp(1:10),'...',filesep,'...',fliplr(tmp(1:50))];
	end
  else
	f_label = RecentFiles{ii};
  end
  f_label = [num2str(ii),'  ',f_label];
  recent_files_h = uimenu(FILEMENU_OPEN_RECENT,...
	'label',f_label,...
	'callback',{@l_OpenFile,'single',RecentFiles{ii}},...
	'userdata',RecentFiles{ii});
end

% Disable the whole menu if it's empty
if isempty(get(FILEMENU_OPEN_RECENT,'children'))
  set(FILEMENU_OPEN_RECENT,'enable','off')
else
  set(FILEMENU_OPEN_RECENT,'enable','on')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Build the plugins menu
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_BuildPluginsMenu(h,evd)
	
	% Clear the plugins menu
	try
		delete(get(H.UIPLUGINS,'children'));
	catch
		% pass
	end
	
	current_dir = pwd;
% Try to read plugins structure
try
  
	plugins=l_FindPlugins(Dat.PluginsFolder);
	counter=0;
	for ii=1:length(plugins.name)
		% Temporarily cd into the plugin directory
		cd(Dat.PluginsFolder)
		
		if plugins.isGroup(ii)
			group_h=uimenu(H.UIPLUGINS,...
				'Label',plugins.groupName{ii},...
				'enable','off');
			cd([Dat.PluginsFolder,plugins.groupFolder{ii}])
			for kk=1:length(plugins.fname{ii})
				set(group_h,'enable','on')
				func_h = str2func(plugins.fname{ii}{kk});
				tmp_h=uimenu(group_h,...
					'Label',plugins.name{ii}{kk},...
					'Userdata',func_h,...
					'callback',@l_ExecutePlugin);
			end
		else
			counter=counter+1;
			func_h = str2func(plugins.fname{ii});
			tmp_h=uimenu(H.UIPLUGINS,...
				'Label',plugins.name{ii},...
				'Userdata',func_h,...
				'callback',@l_ExecutePlugin);
			set(tmp_h,'separator','on')
			if counter==1
				set(tmp_h,'separator','on')
			else
				set(tmp_h,'separator','off')
			end
		end
	end
	cd(current_dir)
	% Add "rescan plugins folder"
	tmp_h=uimenu(H.UIPLUGINS,...
		'Label','Rescan plugins folder',...
		'callback',@l_BuildPluginsMenu);
	set(tmp_h,'separator','on')
catch
	cd(current_dir)
	errordlg(['Error: Could not initialize plugins menu.'],'Error initializing plugins','modal');
end
	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Invert Colormap
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_ColormapInvert(h,evd)
try
if get(h,'value')
  Dat.ColMapInverted = true;
else
  Dat.ColMapInverted = false;
end
Dat.ColMap = flipud(Dat.ColMap);
set(H.FIG,'colormap',Dat.ColMap)
  
catch
  aedes_errordump(lasterror);
end
end % function l_ColormapInvert(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Reset Window/Level sliders
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_ResetWindowLevel(h,evd)
try
% Look for a quick return
if ~strcmpi(get(H.FIG,'selectiontype'),'alt')
  return
end

% Determine which slider called the function
if h==H.WINDOW_SLIDER
  Dat.Window = 100;
  set(H.WINDOW_SLIDER,'value',Dat.Window)
elseif h==H.LEVEL_SLIDER
  Dat.Level = 50;
end

% Set Window/Level
l_SetWindowLevel(0,[])

catch
  aedes_errordump(lasterror);
end
end % function l_ResetWindowLevel(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Functions for ROI manipulation
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw/Erase ROI voxels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_RoiDrawErase(h,evd,opt,wbm,ax_ind,roi_ind)
try
  
if nargin<4
  wbm=false;
end

% Get axis handle
IMAXh = [H.IMAX1,H.IMAX2,H.IMAX3];
ax_h = IMAXh(ax_ind);

% Get current point from the axis
cp = get(ax_h,'currentpoint');
y=round(cp(1));
x=round(cp(3));
vol = Dat.CurrentVol;

%% Check that indices are within image limits
if ax_ind==1
  x=min(max(x,1),Dat.ImageDim(Dat.DataInd,1));
  y=min(max(y,1),Dat.ImageDim(Dat.DataInd,2));
  z=Dat.Slices(3);
  new_voxels = [x y z vol];
elseif ax_ind==2
  x=min(max(x,1),Dat.ImageDim(Dat.DataInd,1));
  y=min(max(y,1),Dat.ImageDim(Dat.DataInd,3));
  z=Dat.Slices(2);
  new_voxels = [x z y vol];
else
  x=min(max(x,1),Dat.ImageDim(Dat.DataInd,2));
  y=min(max(y,1),Dat.ImageDim(Dat.DataInd,3));
  z=Dat.Slices(1);
  new_voxels = [z x y vol];
end


%% Check if some points have been missed
if wbm
  x_old = Dat.x_old;
  y_old = Dat.y_old;
  dx = abs(x_old-x);
  dy = abs(y_old-y);
  if dx>1 | dy>1
    d = max(dx,dy);
    x_new=round(linspace(x_old,x,d+1))';
    y_new=round(linspace(y_old,y,d+1))';
    z_new=ones(d,1)*z;
    vol_new = ones(d,1)*vol;
    if ax_ind==1
      new_voxels=[x_new(2:end) y_new(2:end) z_new vol_new];
    elseif ax_ind==2
      new_voxels=[x_new(2:end) z_new y_new(2:end) vol_new];
    else
      new_voxels=[z_new x_new(2:end) y_new(2:end) vol_new];
    end
    Dat.x_old = x_new(end);
    Dat.y_old = y_new(end);
  end
else
  Dat.x_old = x;
  Dat.y_old = y;
  
  % Make sure that the ROI that is been drawn is also visible
  if ~any(Dat.RoiView==roi_ind)
    set(H.IMROI(:,roi_ind),'visible','on')
    set(H.ROIVIEW_LBOX,'value',[get(H.ROIVIEW_LBOX,'value') roi_ind])
    Dat.RoiView = get(H.ROIVIEW_LBOX,'value');
  end
end


% Convert subscripts to linear indices
if Dat.ImageDim(Dat.DataInd,3)==1 && Dat.ImageDim(Dat.DataInd,4)==1 
  ind=sub2ind(size(ROI(roi_ind).voxels{Dat.DataInd}),new_voxels(:,1),...
              new_voxels(:,2));
elseif Dat.ImageDim(Dat.DataInd,4)==1
  ind=sub2ind(size(ROI(roi_ind).voxels{Dat.DataInd}),new_voxels(:,1),...
              new_voxels(:,2),...
              new_voxels(:,3));
elseif Dat.ImageDim(Dat.DataInd,4)~=1
  ind=sub2ind(size(ROI(roi_ind).voxels{Dat.DataInd}),new_voxels(:,1),...
              new_voxels(:,2),...
              new_voxels(:,3),...
              new_voxels(:,4));
end

% Add/remove ROI voxels
if strcmpi(opt,'erase')
  ROI(roi_ind).voxels{Dat.DataInd}(ind)=false;
else
  ROI(roi_ind).voxels{Dat.DataInd}(ind)=true;
end

% Draw ROI
l_RefreshRoi([1 2 3],roi_ind)

% Update ROI Info
l_UpdateRoiInfoText([],[])

% Update ROI saved status
Dat.RoiSaved = false;

catch
  aedes_errordump(lasterror);
end
end % function l_DrawRemoveRoi(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add new ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_RoiAdd(h,evd,opt,comp_rois)
try
  
if nargin<3 || length(ROI)==0
  opt='new';
end

% Check if new ROI can be added
if length(ROI)==size(Dat.RoiColors,1)
  h=warndlg(['Maximum number of separate colors for ROIs reached. ',...
    'Individual ROIs are no longer visually separable.'],...
    'Maximum number of ROI colors reached','modal');
  uiwait(h)
end

% Ask a label for the new ROI
done=false;
while ~done
  resp=aedes_inputdlg('Please give a descriptive label for the ROI','ROI Label?', ...
                    ['ROI ' num2str(length(ROI)+1)]);
  if isempty(resp)
    return
  else
    resp = resp{1};
  end
  resp = strtrim(resp);
  if isempty(resp)
    return
  end
  if ~isempty(ROI)
    roi_labels = {ROI(:).label};
    if any(strcmp(resp,roi_labels))
      h=warndlg({['A ROI with the label "' resp '" already exists.'],...
                 'Choose another label.'},'Label already exists','modal');
      uiwait(h)
    else
      done=true;
    end
  else
    done=true;
  end
end

le=length(ROI);
roi_ind = le+1;
currentRoi = get(H.ROI_EDIT,'value');

% Add new ROI to the ROI structure
ROI(roi_ind).voxels = cell(1,Dat.DataL);
ROI(roi_ind).fpath = cell(1,Dat.DataL);
ROI(roi_ind).fname = cell(1,Dat.DataL);

for ii=1:Dat.DataL
  if strcmpi(opt,'new') % Add new empty ROI
    ROI(roi_ind).voxels{ii}=false(Dat.ImageDim(ii,1),...
                                  Dat.ImageDim(ii,2),...
                                  Dat.ImageDim(ii,3),...
                                  Dat.ImageDim(ii,4));
  elseif strcmpi(opt,'copy') % Add new copied ROI
    ROI(roi_ind).voxels{ii}=ROI(currentRoi).voxels{ii};
  elseif strcmpi(opt,'union') % Add union of two ROIs
    ROI(roi_ind).voxels{ii} = false(Dat.ImageDim(ii,1),...
                                  Dat.ImageDim(ii,2),...
                                  Dat.ImageDim(ii,3),...
                                  Dat.ImageDim(ii,4));
    ind_a=find(ROI(comp_rois(1)).voxels{ii});
    ind_b=find(ROI(comp_rois(2)).voxels{ii});
    ROI(roi_ind).voxels{ii}(ind_a)=true;
    ROI(roi_ind).voxels{ii}(ind_b)=true;
  elseif strcmpi(opt,'intersect') % Add intersection of two ROIs
    ROI(roi_ind).voxels{ii} = false(Dat.ImageDim(ii,1),...
                                  Dat.ImageDim(ii,2),...
                                  Dat.ImageDim(ii,3),...
                                  Dat.ImageDim(ii,4));
    ind_a=find(ROI(comp_rois(1)).voxels{ii});
    ind_b=find(ROI(comp_rois(2)).voxels{ii});
    intersect_ind = intersect(ind_a,ind_b);
    ROI(roi_ind).voxels{ii}(intersect_ind)=true;
  elseif strcmpi(opt,'subtract') % Add subtraction of two ROIs
     ROI(roi_ind).voxels{ii} = false(Dat.ImageDim(ii,1),...
                                  Dat.ImageDim(ii,2),...
                                  Dat.ImageDim(ii,3),...
                                  Dat.ImageDim(ii,4));
    ind_a=find(ROI(comp_rois(1)).voxels{ii});
    ind_b=find(ROI(comp_rois(2)).voxels{ii});
    setdiff_ind = setdiff(ind_a,ind_b);
    ROI(roi_ind).voxels{ii}(setdiff_ind)=true;
  elseif strcmpi(opt,'xor') % Add exclusive disjunction of two ROIs
    ROI(roi_ind).voxels{ii} = false(Dat.ImageDim(ii,1),...
                                  Dat.ImageDim(ii,2),...
                                  Dat.ImageDim(ii,3),...
                                  Dat.ImageDim(ii,4));
    ind_a=find(ROI(comp_rois(1)).voxels{ii});
    ind_b=find(ROI(comp_rois(2)).voxels{ii});
    setxor_ind = setxor(ind_a,ind_b);
    ROI(roi_ind).voxels{ii}(setxor_ind)=true;
  end
  ROI(roi_ind).fpath{ii} = DATA{ii}.HDR.fpath;
  ROI(roi_ind).fname{ii} = DATA{ii}.HDR.fname;
end
ROI(roi_ind).label = resp;

% Set the color of the new ROI to the first
% free color available
if le>0
  free_colors=Dat.RoiColors(find(ismember(Dat.RoiColors,...
    reshape([ROI(:).color],3,le)', ...
    'rows')==0),:);
  if isempty(free_colors)
    col_ind = roi_ind-floor(roi_ind/size(Dat.RoiColors,1))*size(Dat.RoiColors,1);
    if col_ind==0
      col_ind = size(Dat.RoiColors,1);
    end
    ROI(roi_ind).color = Dat.RoiColors(col_ind,:);
  else
    ROI(roi_ind).color = free_colors(1,:);
  end
else
  ROI(roi_ind).color = Dat.RoiColors(1,:);
end

  
% Create axes and image for the ROI...
for ii=1:3
  H.ROIAX(ii,roi_ind) = axes('parent',H.FIG,...
                             'units',eval(['get(H.IMAX' sprintf('%d',ii) ',''units'')']),...
                             'Position',...
                             eval(['get(H.IMAX' num2str(ii) ',''position'')']),...
                             'Ylim',eval(['get(H.IMAX' sprintf('%d',ii) ',''ylim'')']), ...
                             'xlim',eval(['get(H.IMAX' sprintf('%d',ii) ',''xlim'')']),...
                             'Yticklabel',{},...
                             'Xticklabel',{},...
                             'ydir','reverse',...
                             'Box','on', ...
                             'visible','off',...
                             'hittest','off',...
                             'alimmode','manual');
  H.IMROI(ii,roi_ind)=image('parent',H.ROIAX(ii,roi_ind),...
                            'cdata',[],'visible','off',...
                            'hittest','off');
	if Dat.HG2graphics
		set(H.IMROI(ii,roi_ind),'hittest','on',...
			'buttondownfcn',@l_SetMouseGestures)
	end
end


% Set New ROI label to the View listbox and
% Roi Edit popupmenu
str=get(H.ROIVIEW_LBOX,'string');
val = get(H.ROIVIEW_LBOX,'value');
if isempty(str)
  set(H.ROIVIEW_LBOX,'string',{resp},...
                    'value',[roi_ind])
  set(H.ROI_EDIT,'string',{resp},...
                 'value',roi_ind)
else
  set(H.ROIVIEW_LBOX,'string',{str{:} resp},...
                    'value',[val roi_ind])
  set(H.ROI_EDIT,'string',{str{:} resp},...
                 'value',roi_ind)
end

% Enable delete buttons, edit popup, and ROI save uimenu
set(H.UICH_ROIENABLED,'enable','on')

% Set viewed ROI(s)
Dat.RoiView = get(H.ROIVIEW_LBOX,'value');

% Set mouse for drawing ROIs
set(get(H.UITOGGLE_DRAWROI,'userdata'),'state','off')
set(H.UITOGGLE_DRAWROI,'state','on','enable','on')

l_DisplayData([],[])

% Draw crossbars
l_ShowHideCrossbars([],[])

% Update ROI Info
l_UpdateRoiInfoText([],[])

% ROI changed
Dat.RoiSaved=false;

catch
  aedes_errordump(lasterror);
end
end % function l_RoiAdd(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delete ROI(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_RoiDelete(h,evd,opt)
try
%% Delete selected or delete all?
if strcmpi(opt,'selected') % Delete selected
  
  % Get indices to the ROIs to be removed
  str = get(H.ROIVIEW_LBOX,'string');
  val = get(H.ROIVIEW_LBOX,'value');
  if isempty(val)
    return
  end
  SelRoiLabels = {str{val}};
  
  resp=questdlg({'The ROI(s) with the following labels will be deleted:','',...
                 SelRoiLabels{:},'','Are you sure?'},...
                'Delete selected ROI(s)?','Yes','No','No');
  if isempty(resp) || strcmpi(resp,'No')
    return
  end
  
  
  % Check if we are removing all ROIs
  if length(val)==length(str)
    rm_all=true;
  else
    rm_all=false;
  end
else
  if strcmpi(opt,'all')
    resp=questdlg('Delete ALL ROI(s). Are you sure?',...
                  'Delete ALL ROI(s)','Yes','No','No');
    if isempty(resp) || strcmpi(resp,'No')
      return
    end
  end
  rm_all=true;
end

% Remove all
if rm_all
  ROI = [];
  delete(H.ROIAX)
  H.ROIAX=[];
  H.IMROI=[];
  set(H.ROIVIEW_LBOX,'value',[],'string','')
  set(H.UICH_ROIENABLED,'enable','off')
  set(H.ROI_EDIT,'value',1,'string',{''},'enable','off')
  
  % Clear undo buffer
  l_RoiUndoBuffer('clear')
  
  
  
  % Set mouse for changing view
  set(get(H.UITOGGLE_ARROW,'userdata'),'state','off')
  set(H.UITOGGLE_ARROW,'state','on')
  
  % No ROI(s) exist, so they don't have to be saved
  Dat.RoiSaved=true;
  
  % Because image transparency is not needed anymore, 
  % the figure renderer can be set back to painters (faster).
  %set(H.FIG,'renderer','painters')
else % Remove some ROIs (but not all of them)
  new_str=str;
  new_str(val)=[];
  delete(H.ROIAX(:,val))
  H.ROIAX(:,val) = [];
  H.IMROI(:,val) = [];
  ROI(val) = [];
  
  set(H.ROIVIEW_LBOX,'value',[],'string',new_str)
  tmp_str=get(H.ROI_EDIT,'string');
  tmp_val=get(H.ROI_EDIT,'value');
  tmp=tmp_str(tmp_val);
  if any(strcmp(tmp,new_str))
    edit_val=find(strcmpi(tmp,new_str));
  else
    edit_val = 1;
  end
  set(H.ROI_EDIT,'value',edit_val,'string',new_str)
  
  % Correct ROI UNDO indices
  if ~isempty(Dat.RoiUndo)
    ind=find(ismember({Dat.RoiUndo(:).label},{str{val}}));
    l_RoiUndoBuffer('remove',ind)
  end
  
  
end
Dat.RoiView=[];

% Redraw crosshairs
try
	delete(H.CROSSBAR_LN)
	H = rmfield(H,'CROSSBAR_LN');
end
l_ShowHideCrossbars([],[]);

% Refresh ROI(s)
l_RefreshRoi;

% Update ROI Info
l_UpdateRoiInfoText([],[])

% ROI changed
Dat.RoiSaved=false;



catch
  aedes_errordump(lasterror);
end
end % function l_RoiDelete(h,


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fill/fill-erase ROIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_RoiFill(h,evd,opt,ax_ind,roi_ind)
try
  
% Get axis handle
IMAXh = [H.IMAX1,H.IMAX2,H.IMAX3];
ax_h = IMAXh(ax_ind);

% Get current point from the axis
cp = get(ax_h,'currentpoint');
y=round(cp(1));
x=round(cp(3));
%vol = Dat.CurrentVol;

%% Check that indices are within image limits
if ax_ind==1
  x=min(max(x,1),Dat.ImageDim(Dat.DataInd,1));
  y=min(max(y,1),Dat.ImageDim(Dat.DataInd,2));
  I = ROI(roi_ind).voxels{Dat.DataInd}(:,:,Dat.Slices(3),Dat.CurrentVol);
elseif ax_ind==2
  x=min(max(x,1),Dat.ImageDim(Dat.DataInd,1));
  y=min(max(y,1),Dat.ImageDim(Dat.DataInd,3));
  I = reshape(ROI(roi_ind).voxels{Dat.DataInd}(:,Dat.Slices(2),:,Dat.CurrentVol),...
              Dat.ImageDim(Dat.DataInd,1),Dat.ImageDim(Dat.DataInd,3));
else
  x=min(max(x,1),Dat.ImageDim(Dat.DataInd,2));
  y=min(max(y,1),Dat.ImageDim(Dat.DataInd,3));
  I = reshape(ROI(roi_ind).voxels{Dat.DataInd}(Dat.Slices(1),:,:,Dat.CurrentVol),...
              Dat.ImageDim(Dat.DataInd,2),Dat.ImageDim(Dat.DataInd,3));
end

% Flood fill binary (logical) image starting from
% mouse click point... 
if strcmp(opt,'fill')
  %I=logical(I);
elseif strcmp(opt,'fill_remove')
  I=I==false;
  %I=logical(I);
end

if I(x,y)
	% Nothing to fill...
	return
end

% Add undo information to buffer
l_RoiUndoBuffer('add',roi_ind,ax_ind);

% Check if Image Processing toolbox is installed and we if could use imfill
% function which is much faster than aedes_roifill for doing the flood fill
% operation.
try
  if Dat.isImageProc
    I=imfill(I,[x y],4);
  else
    I=aedes_roifill(I,[x y]);
  end
catch
  % Fallback to aedes_roifill...
  I=aedes_roifill(I,[x y]);
end

ind=find(I);
tmp=false(size(I));
tmp(ind)=true;
if strcmp(opt,'fill_remove')
  tmp=tmp==false;
end
if ax_ind==1
  ROI(roi_ind).voxels{Dat.DataInd}(:,:,Dat.Slices(3),Dat.CurrentVol)=tmp;
elseif ax_ind==2
  ROI(roi_ind).voxels{Dat.DataInd}(:,Dat.Slices(2),:,Dat.CurrentVol)=tmp;
else
  ROI(roi_ind).voxels{Dat.DataInd}(Dat.Slices(1),:,:,Dat.CurrentVol)=tmp;
end

% Make sure that the ROI that is been drawn is also visible
if ~any(Dat.RoiView==roi_ind)
  set(H.IMROI(:,roi_ind),'visible','on')
  set(H.ROIVIEW_LBOX,'value',[get(H.ROIVIEW_LBOX,'value') roi_ind])
  Dat.RoiView = get(H.ROIVIEW_LBOX,'value');
end

% Refresh ROI images
%l_RefreshRoi(ax_ind,roi_ind);
l_RefreshRoi();

% Update ROI Info
l_UpdateRoiInfoText([],[])

% ROI Edges
if Dat.ShowRoiEdges
  l_ShowRoiEdges([],[],'')
end

% ROI changed
Dat.RoiSaved=false;

catch
  aedes_errordump(lasterror);
end
end % function l_RoiFill(h,

%%%%%%%%%%%%%%%%%%%%%%%%%
% View selected ROI(s)
%%%%%%%%%%%%%%%%%%%%%%%%%
function l_RoiView(h,evd)
try
  
% Get selected ROI(s) i.e. listbox values
Dat.RoiView = get(H.ROIVIEW_LBOX,'value');

% Set visible off for unselected ROI(s)
if isempty(Dat.RoiView)
  set(H.IMROI,'visible','off')
else
  tmp=1:length(ROI);
  set(H.IMROI(:,ismember(tmp,Dat.RoiView)==0),'visible','off')
  set(H.IMROI(:,Dat.RoiView),'visible','on')
end

% Refresh ROI(s)
l_RefreshRoi;

% Draw ROI edges
if Dat.ShowRoiEdges
  l_ShowRoiEdges([],[],'all')
end

% If only one ROI is selected from the listbox, change the value of the
% edit/current ROI
if length(Dat.RoiView)==1
  set(H.ROI_EDIT,'value',Dat.RoiView)
end

return

catch
  aedes_errordump(lasterror);
end
end % function l_RoiView(h,

%%%%%%%%%%%%%%%%%%%%%%%%
% Undo last draw action
%%%%%%%%%%%%%%%%%%%%%%%%
function l_RoiUndo(h,evd)
try
  
if isempty(Dat.RoiUndo)
  return
end

str = get(H.ROIVIEW_LBOX,'string');
ax_ind = Dat.RoiUndo(end).ax_ind;
roi_ind = find(strcmp(Dat.RoiUndo(end).label,str));
slice = Dat.RoiUndo(end).slice;
vol = Dat.RoiUndo(end).vol;
dataind = Dat.RoiUndo(end).DataInd;

% Check if the undo buffer contains indices or slices
if ax_ind==0
  for ii=1:length(Dat.RoiUndo(end).voxels)
    ROI(roi_ind).voxels{ii}(:,:,:,vol)=false;
    ROI(roi_ind).voxels{ii}(Dat.RoiUndo(end).voxels{ii}...
                            +prod(Dat.ImageDim(ii,1:3))* ...
                            (vol-1))=true;
  end
else
  if ax_ind==1
    ROI(roi_ind).voxels{dataind}(:,:,slice,vol)=Dat.RoiUndo(end).voxels;
  elseif ax_ind==2
    ROI(roi_ind).voxels{dataind}(:,slice,:,vol)=Dat.RoiUndo(end).voxels;
  else
    ROI(roi_ind).voxels{dataind}(slice,:,:,vol)=Dat.RoiUndo(end).voxels;
  end
end

% Remove entry from the undo buffer
l_RoiUndoBuffer('remove',length(Dat.RoiUndo))


% ROI changed
Dat.RoiSaved=0;

% Refresh ROI(s)
l_RefreshRoi([1 2 3],roi_ind)

% Update ROI Info
l_UpdateRoiInfoText([],[])

% ROI Edges
if Dat.ShowRoiEdges
  l_ShowRoiEdges([],[],'all')
end

% ROI changed
Dat.RoiSaved=false;

catch
  aedes_errordump(lasterror);
end
end % function l_RoiUndo(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Move ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_RoiMove(h,evd,opt,wbm,ax_ind,roi_ind)
try
  
% Get axis handle
IMAXh = [H.IMAX1,H.IMAX2,H.IMAX3];
ax_h = IMAXh(ax_ind);

% Get initial point
if wbm==0
  cp=round(get(ax_h,'currentpoint'));
  y=cp(1);
  x=cp(3);
  vol = Dat.CurrentVol;
  
  Dat.RoiMoveInitialPoint = [x y];
  Dat.RoiMoveOrigSlice = [];
  if strcmpi(opt,'move2d') || Dat.isDataMixed
    if ax_ind==1
      if strcmpi(opt,'move3d') && Dat.isDataMixed
        %Dat.RoiMoveOrigSlice = ROI(roi_ind).voxels{Dat.DataInd}(:,:,:,vol);
        tmp_J = [];
        tmp_I = [];
        count=1;
        for ii=1:length(ROI(roi_ind).voxels)
          [I,J]=find(ROI(roi_ind).voxels{ii});
          if ~isempty(I) & ~isempty(J)
            tmp_J(count,:) = [-min(J)+1,Dat.ImageDim(ii,2)-max(J)];
            tmp_I(count,:) = [-min(I)+1,Dat.ImageDim(ii,1)- ...
                              max(I)];
            count=count+1;
          end
          Dat.RoiMoveOrigSlice{ii} = ROI(roi_ind).voxels{ii};
        end
        Dat.RoiMoveLimits = [max(tmp_J(:,1)),min(tmp_J(:,2));...
                            max(tmp_I(:,1)),min(tmp_I(:,2))];
      else
          
        Dat.RoiMoveOrigSlice = ROI(roi_ind).voxels{Dat.DataInd}(:,:, ...
                                                          Dat.Slices(3),vol);
        [I,J]=find(Dat.RoiMoveOrigSlice);
        Dat.RoiMoveLimits = [-min(J)+1,Dat.ImageDim(Dat.DataInd,2)-max(J);...
                            -min(I)+1,Dat.ImageDim(Dat.DataInd,1)- ...
                            max(I)];
      end
    elseif ax_ind==2
      Dat.RoiMoveOrigSlice = reshape(ROI(roi_ind).voxels{Dat.DataInd}(:,Dat.Slices(2), ...
                                                        :,vol),...
                                     Dat.ImageDim(Dat.DataInd,1), ...
                                     Dat.ImageDim(Dat.DataInd,3));
      [I,J]=find(Dat.RoiMoveOrigSlice);
      Dat.RoiMoveLimits = [-min(J)+1,Dat.ImageDim(Dat.DataInd,3)-max(J);...
                          -min(I)+1,Dat.ImageDim(Dat.DataInd,1)-max(I)];
    else
      Dat.RoiMoveOrigSlice = reshape(ROI(roi_ind).voxels{Dat.DataInd}(Dat.Slices(1),:,:,vol),...
                                     Dat.ImageDim(Dat.DataInd,2),Dat.ImageDim(Dat.DataInd,3));
      [I,J]=find(Dat.RoiMoveOrigSlice);
      Dat.RoiMoveLimits = [-min(J)+1,Dat.ImageDim(Dat.DataInd,3)-max(J);...
                          -min(I)+1,Dat.ImageDim(Dat.DataInd,2)-max(I)];
    end
  else %% Move 3D
    Dat.RoiMoveOrigSlice = ROI(roi_ind).voxels{Dat.DataInd}(:,:,:,vol);
    tmp_J = [];
    tmp_I = [];
    count=1;
    if ax_ind==1
      for ii=1:size(ROI(roi_ind).voxels{Dat.DataInd},3)
        [I,J]=find(ROI(roi_ind).voxels{Dat.DataInd}(:,:,ii,vol));
        if ~isempty(I) & ~isempty(J)
          tmp_J(count,:) = [-min(J)+1,Dat.ImageDim(Dat.DataInd,2)-max(J)];
          tmp_I(count,:) = [-min(I)+1,Dat.ImageDim(Dat.DataInd,1)- ...
                            max(I)];
          count=count+1;
        end
      end
      Dat.RoiMoveLimits = [max(tmp_J(:,1)),min(tmp_J(:,2));...
                          max(tmp_I(:,1)),min(tmp_I(:,2))];
    elseif ax_ind==2
      for ii=1:size(ROI(roi_ind).voxels{Dat.DataInd},2)
        [I,J]=find(squeeze(ROI(roi_ind).voxels{Dat.DataInd}(:,ii,:, ...
                                                          vol)));
        if ~isempty(I) & ~isempty(J)
          tmp_J(count,:) = [-min(J)+1,Dat.ImageDim(Dat.DataInd,3)-max(J)];
          tmp_I(count,:) = [-min(I)+1,Dat.ImageDim(Dat.DataInd,1)- ...
                            max(I)];
          count=count+1;
        end
      end
      Dat.RoiMoveLimits = [max(tmp_J(:,1)),min(tmp_J(:,2));...
                          max(tmp_I(:,1)),min(tmp_I(:,2))];
    else
      for ii=1:size(ROI(roi_ind).voxels{Dat.DataInd},1)
        [I,J]=find(squeeze(ROI(roi_ind).voxels{Dat.DataInd}(ii,:,:, ...
                                                          vol)));
        if ~isempty(I) & ~isempty(J)
          tmp_J(count,:) = [-min(J)+1,Dat.ImageDim(Dat.DataInd,3)-max(J)];
          tmp_I(count,:) = [-min(I)+1,Dat.ImageDim(Dat.DataInd,2)-max(I)];
          count=count+1;
        end
      end
      Dat.RoiMoveLimits = [max(tmp_J(:,1)),min(tmp_J(:,2));...
                          max(tmp_I(:,1)),min(tmp_I(:,2))];
    end
  end
  set(H.FIG,'pointer','fleur')
  return
end

% Check that moved ROI is not empty
if ~iscell(Dat.RoiMoveOrigSlice)
  if isempty(find(Dat.RoiMoveOrigSlice,1,'first'))
    return
  end
else
  if isempty(find(Dat.RoiMoveOrigSlice{Dat.DataInd},1,'first'))
    return
  end
end

% Get current point from the axis
cp = get(ax_h,'currentpoint');
y=round(cp(1));
x=round(cp(3));
vol = Dat.CurrentVol;

% Compare to the initial point
dfx = diff([Dat.RoiMoveInitialPoint(1) x]);
dfy = diff([Dat.RoiMoveInitialPoint(2) y]);
if dfy<Dat.RoiMoveLimits(1,1)
  dfy=Dat.RoiMoveLimits(1,1);
elseif dfy>Dat.RoiMoveLimits(1,2)
  dfy=Dat.RoiMoveLimits(1,2);
end
if dfx<Dat.RoiMoveLimits(2,1)
  dfx=Dat.RoiMoveLimits(2,1);
elseif dfx>Dat.RoiMoveLimits(2,2)
  dfx=Dat.RoiMoveLimits(2,2);
end

% Move ROI in the current slice by dfx and dfy
% pixels
if strcmp(opt,'move2d') || Dat.isDataMixed
  if ax_ind==1
    if strcmp(opt,'move3d') && Dat.isDataMixed
      for ii=1:Dat.DataL
        ROI(roi_ind).voxels{ii}(:,:,:,vol)=...
            circshift(Dat.RoiMoveOrigSlice{ii},[dfx,dfy]);
      end
    else
      ROI(roi_ind).voxels{Dat.DataInd}(:,:,Dat.Slices(3),vol)=...
          circshift(Dat.RoiMoveOrigSlice,[dfx,dfy]);
    end
  elseif ax_ind==2
    ROI(roi_ind).voxels{Dat.DataInd}(:,Dat.Slices(2),:,vol)=...
        circshift(Dat.RoiMoveOrigSlice,[dfx,dfy]);
  else
    ROI(roi_ind).voxels{Dat.DataInd}(Dat.Slices(1),:,:,vol)=...
        circshift(Dat.RoiMoveOrigSlice,[dfx,dfy]);
  end
  
else
  % Move whole ROI in 3D by dfx and dfy
  % pixels
  if ax_ind==1
    ROI(roi_ind).voxels{Dat.DataInd}(:,:,:,vol)=...
        circshift(Dat.RoiMoveOrigSlice,[dfx,dfy,0]);
  elseif ax_ind==2
    ROI(roi_ind).voxels{Dat.DataInd}(:,:,:,vol)=...
        circshift(Dat.RoiMoveOrigSlice,[dfx,0,dfy]);
  else
    ROI(roi_ind).voxels{Dat.DataInd}(:,:,:,vol)=...
        circshift(Dat.RoiMoveOrigSlice,[0,dfx,dfy]);
  end  
end

% Make sure that the ROI that is been drawn is also visible
if ~any(Dat.RoiView==roi_ind)
  set(H.IMROI(:,roi_ind),'visible','on')
  set(H.ROIVIEW_LBOX,'value',[get(H.ROIVIEW_LBOX,'value') roi_ind])
  Dat.RoiView = get(H.ROIVIEW_LBOX,'value');
end

% ROI changed
Dat.RoiSaved=false;

% Refresh ROI
l_RefreshRoi([1 2 3],roi_ind)

% Update ROI Info
l_UpdateRoiInfoText([],[])

% ROI changed
Dat.RoiSaved=false;

catch
  aedes_errordump(lasterror);
end
end % function l_MoveRoi(h,

%%%%%%%%%%%%%%%%%%%
% Copy ROI
%%%%%%%%%%%%%%%%%%%
function l_RoiCopy(h,evd)
try
  
roi_ind = get(H.ROI_EDIT,'value');
[slice_ind,dir_ind,roilabel,copytype]=aedes_roi_copy_gui(ROI,roi_ind);

if isempty(slice_ind)
  %% Action possibly canceled
  return
end
[h,txh]=aedes_calc_wait('Copying ROI(s)...');

for ii=1:length(roilabel)
  % Add current ROI to undo buffer
  l_RoiUndoBuffer('add_indices',find(strcmp({ROI(:).label},roilabel{ii})),[]);
  ROI = aedes_copy_roi(ROI,roilabel{ii},dir_ind,slice_ind,Dat.CurrentVol, ...
                 Dat.DataInd,copytype);
end

% Refresh ROI(s)
l_RefreshRoi;
delete(h)

% Update ROI Info
l_UpdateRoiInfoText([],[])

% ROI changed
Dat.RoiSaved=false;

catch
  aedes_errordump(lasterror);
end
end % function l_RoiCopy(h,


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copy current ROI slice
%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function l_RoiCopyCurrent(h,evd)
		
		if isempty(ROI)
			return
		end
		
		% Get current ROI
		currentRoi = get(H.ROI_EDIT,'value');
		
		if Dat.AxView~=0
			CopyDir = Dat.AxView;
		else
			resp = questdlg('Select source slice direction to copy from','ROI copy: select source direction',...
				'X','Y','Z','Z');
			if isempty(resp)
				return
			end
			if strcmpi(resp,'X')
				CopyDir = 1;
			elseif strcmpi(resp,'Y')
				CopyDir = 2;
			else
				CopyDir = 3;
			end
		end

		if CopyDir==3
			CurrentROISlice = ROI(currentRoi).voxels{Dat.DataInd}(:,:,Dat.Slices(3),Dat.CurrentVol);
		elseif CopyDir==2
			CurrentROISlice = reshape(ROI(currentRoi).voxels{Dat.DataInd}(:,Dat.Slices(2),:,Dat.CurrentVol),...
				Dat.ImageDim(Dat.DataInd,1),Dat.ImageDim(Dat.DataInd,3));
		else
			CurrentROISlice = reshape(ROI(currentRoi).voxels{Dat.DataInd}(Dat.Slices(1),:,:,Dat.CurrentVol),...
				Dat.ImageDim(Dat.DataInd,2),Dat.ImageDim(Dat.DataInd,3));
		end
		Dat.CopiedROI = CurrentROISlice;
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paste copied ROI slice
%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function l_RoiPasteCopied(h,evd)
		
		if isempty(ROI) || isempty(Dat.CopiedROI)
			return
		end
		
		if Dat.AxView==0
			resp = questdlg('Select target slice direction to paste ROI into.',...
				'ROI paste: select target direction',...
				'X','Y','Z','Z');
			if isempty(resp)
				return
			end
			if strcmpi(resp,'X')
				PasteDir = 1;
			elseif strcmpi(resp,'Y')
				PasteDir = 2;
			else
				PasteDir = 3;
			end
		else
			PasteDir = Dat.AxView;
		end
		
		% Get current ROI
		currentRoi = get(H.ROI_EDIT,'value');
		
		% Get target size
		if PasteDir==1
			target_sz = [Dat.ImageDim(2),...
				Dat.ImageDim(3)];
		elseif PasteDir==2
			target_sz = [Dat.ImageDim(1),...
				Dat.ImageDim(3)];
		else
			target_sz = [Dat.ImageDim(1),...
				Dat.ImageDim(2)];
		end
		
		% Check that target and source sizes match
		source_sz = size(Dat.CopiedROI);
		if ~all(source_sz==target_sz)
			if Dat.AxView==0
				errordlg('Size mismatch! Cannot paste copied ROI slice to selected slice direction.',...
					'Size mismatch.');
			end
			return
		end
		
		% Add current state to undo buffer
		ind = [3 2 1];
		l_RoiUndoBuffer('add',currentRoi,ind(PasteDir));
		
		% Copy ROI to current position
		if Dat.AxView==1
			ROI(currentRoi).voxels{Dat.DataInd}(Dat.Slices(1),:,:,Dat.CurrentVol)=Dat.CopiedROI;
		elseif Dat.AxView==2
			ROI(currentRoi).voxels{Dat.DataInd}(:,Dat.Slices(2),:,Dat.CurrentVol)=Dat.CopiedROI;
		else
			ROI(currentRoi).voxels{Dat.DataInd}(:,:,Dat.Slices(3),Dat.CurrentVol)=Dat.CopiedROI;
		end
		
		% Refresh ROIs
		l_RefreshRoi;
		
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save ROI(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_RoiSave(h,evd,opt)
try
  
if ~exist('opt','var')
  SaveTemplate = false;
else
  SaveTemplate = true;
end

% Default name for file
[fp,fn,fe]=fileparts([Dat.HDR{1}.fpath,Dat.HDR{1}.fname]);
filename = fn;

% If the file is a VNMR file, extract the filename from the path name
if strcmpi(fn,'fid')
  Ind=find(fp==filesep);
  if isempty(Ind)
    if ispc
      Ind=find(fp=='/');
    else
      Ind=find(fp=='\');
    end
  end
  filename = fp(Ind(end)+1:end);
  filename = filename(1:end-4);
end

% Default path for ROI file
if SaveTemplate
  try
    filepath = getpref('Aedes','PutRoiFileDir');
  catch
    filepath = '';
  end
  filename = [filename,'_template'];
  DialogTitle = 'Save ROI Template File As';
	
	% Prompt to save either current slice or volume
	if ~Dat.isDataMixed
		resp = questdlg({'ROI templates can be loaded on data whose size differs',...
			'from than of the ROI. 2D templates (current slice) can be',...
			'loaded on any data whereas 3D templates (current volume)',...
			'can be loaded on any 4D data with the same 3D resolution.',...
			'',...
			'Do you want to save 2D or 3D template?'},...
			'Save 2D or 3D template?',...
			'2D (current slice)','3D (current volume)','Cancel',...
			'2D (current slice)');
		if isempty(resp) || strcmpi(resp,'Cancel')
			return
		elseif strcmpi(resp,'2D (current slice)')
			Save3Dtemplate = false;
		else
			Save3Dtemplate = true;
		end
	else
		Save3Dtemplate = false;
	end
else
  try
    filepath = getpref('Aedes','PutRoiFileDir');
  catch
    filepath = '';
  end
  DialogTitle = 'Save ROI(s) As';
end

ok=false;
while ~ok
  % Ask file name
  [f_name,f_path,f_index]=uiputfile({'*.roi;*.ROI',...
                      'Aedes ROI-Files (*.roi)';...
                      '*.*','All Files (*.*)'},...
                                    DialogTitle,...
                                    [filepath,filename]);
  
  % Return if cancel is pressed
  if ( all(f_name==0) | all(f_path==0) )
    return
  end
  savefilename = [f_path,f_name];
  
  % Save ROI(s)
  if SaveTemplate
    if Dat.isDataMixed
      tmp_data = {DATA{Dat.DataInd}};
      for ii=1:length(ROI)
        tmp_roi(ii).voxels = {ROI(ii).voxels{Dat.DataInd}};
        tmp_roi(ii).fpath = {ROI(ii).fpath{Dat.DataInd}};
        tmp_roi(ii).fname = {ROI(ii).fname{Dat.DataInd}};
        tmp_roi(ii).label = ROI(ii).label;
        tmp_roi(ii).color = ROI(ii).color;
      end
		else
			if Save3Dtemplate
				tmp_data.FTDATA = DATA{1}.FTDATA(:,:,:,Dat.CurrentVol);
				tmp_data.HDR =  DATA{1}.HDR;
			else
				tmp_data.FTDATA = DATA{1}.FTDATA(:,:,Dat.Slices(3),Dat.CurrentVol);
				tmp_data.HDR =  DATA{1}.HDR;
			end
			for ii=1:length(ROI)
				if Save3Dtemplate
					tmp_roi(ii).voxels = {ROI(ii).voxels{Dat.DataInd}(:,:,:,Dat.CurrentVol)};
				else
					tmp_roi(ii).voxels = {ROI(ii).voxels{Dat.DataInd}(:,:,Dat.Slices(3),Dat.CurrentVol)};
				end
        tmp_roi(ii).fpath = {ROI(ii).fpath{Dat.DataInd}};
        tmp_roi(ii).fname = {ROI(ii).fname{Dat.DataInd}};
        tmp_roi(ii).label = ROI(ii).label;
        tmp_roi(ii).color = ROI(ii).color;
      end
		end
    [done,msg]=aedes_saveres(tmp_data,tmp_roi,savefilename,'SaveType','roi',...
                       'waitbar',true,'ConfirmOverwrite',false);
  else
	if Dat.isDataMixed
	  rotateflip = [Dat.DataRotation;Dat.DataFlip];
	else
	  rotateflip = Dat.RotateFlip3d;
	end
    [done,msg]=aedes_saveres(DATA,ROI,savefilename,'SaveType','roi',...
                       'waitbar',true,'rotateflip',rotateflip,...
                       'ConfirmOverwrite',false);
  end
  if ~done
    if strcmpi(msg,'Overwrite cancel')
      ok=false;
    else
      h=errordlg(msg,'Could not save ROI(s)','modal');
      return
    end
  else
    ok=true;
  end
end

% Store new directory in preferences
setpref('Aedes','PutRoiFileDir',f_path)

% Set saved state for ROIs
Dat.RoiSaved = true;



catch
  aedes_errordump(lasterror);
end
end % function l_RoiSave(h,


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load ROI(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function l_RoiLoad(h,evd,opt)
    try

      hCalcWait = [];
      resp_append = [];


      if strcmpi(opt,'interactive')

        % Default path for ROI file
        try
          filepath = getpref('Aedes','GetRoiFileDir');
        catch
          filepath = '';
        end

        % Ask for the ROI file
        [f_name,f_path,f_index]=uigetfile({'*.roi;*.ROI;*.nii;*.NII;*.hdr;*.HDR',...
          'ROI-Files (*.roi,*.nii,*,hdr)';...
          '*.*','All Files (*.*)'},...
          'Load ROI(s)',...
          filepath);

        % Return if cancel is pressed
        if ( all(f_name==0) | all(f_path==0) )
          return
        end

        % Display warning if ROIs already exist
        %   if ~isempty(ROI)
        %     resp=questdlg({'All existing ROIs have to be deleted before loading ROI(s) from file.',...
        %                    '',...
        %                    'By choosing "Yes" all existing ROI(s) will be deleted.',...
        %                    'Choose "No" if you want to abort.'},...
        %                   'Delete existing ROI(s) and continue loading?','Yes','No','No');
        %     if strcmpi(resp,'No')
        %       return
        %     end
        %   end
        if ~isempty(ROI)
          resp_append=questdlg({'Do you want to append the loaded ROI(s) to current ROI(s)',...'
            'or replace all the existing ROI(s)?',...
            '',...
            'By choosing "Replace" all existing ROI(s) will be deleted first.',...
            'By choosing "Append" the loaded ROI(s) will be appended to existing ROI(s).',...
            'Choose "Cancel" if you want to abort.'},...
            'Replace ROI(s) or append to current?','Replace','Append','Cancel','Append');
          if isempty(resp_append) || strcmpi(resp_append,'Cancel')
            return
          end
        end

        % Show aedes_calc_wait
        [hCalcWait,txh]=aedes_calc_wait({'Loading ROI(s) from',['"',f_path,f_name,'"']});

        % Store new directory in preferences
        setpref('Aedes','GetRoiFileDir',f_path)
        
        % Get file parts
        [fp,fn,fe]=fileparts([f_path,f_name]);
        
        % Load ROI(s)
        try
          if strcmpi(fe,'.roi')
            tmp=load([f_path,f_name],'-mat');
          else
            [tmp_data,msg]=aedes_read_nifti([f_path,f_name]);
            if isempty(tmp_data)
              error(msg)
            else
              if ~strcmpi(class(tmp_data.FTDATA),'uint8') || ...
                  min(tmp_data.FTDATA(:))~=0 || max(tmp_data.FTDATA(:))~=1
                error('Not a valid NIfTI/Analyze75 ROI file.')
              end
            end
            
            % Construct a ROI structure
            tmp.ROI = [];
            %tmp.FileInfo = [];
            tmp.ROI.voxels{1} = logical(tmp_data.FTDATA);
            tmp.ROI.fpath{1} = [fp,filesep];
            tmp.ROI.fname{1} = f_name;
            if ~isempty(fn)
              tmp.ROI.label = fn;
            else
              tmp.ROI.label = 'NIfTI ROI';
            end
            tmp.ROI.color = [255 0 0];
            %tmp.FileInfo.DataFileName{1} = '';
            %tmp.FileInfo.DataPathName{1} = '';
            clear tmp_data
          end
        catch
          h=errordlg({'Could not load ROI file',['"',f_path,f_name,'"'],...
            '','The following error was returned:',...
            lasterr},'Could not load ROI(s)','modal');
          delete(hCalcWait)
          return
        end

        % Check that the MAT file includes necessary fields
        if ~isempty(tmp) && isstruct(tmp)
          fldnames=fieldnames(tmp);
          if ~any(strcmp(fldnames,'ROI'))%length(find(ismember(fldnames,{'ROI','FileInfo'})))<2
            h=errordlg({'Required fields were not found from the file',...
              [f_path,f_name]},'Required field not found!','modal');
            clear tmp
            delete(hCalcWait)
            return
          end
        else
          h=errordlg({'Required fields were not found from the file',...
            ['"',f_path,f_name,'"'],'',...
            'The file does not seem to be a valid ROI file'},...
            'Required field not found!','modal');
          clear tmp
          delete(hCalcWait)
          return
        end

        % Delete all existing ROI(s) if replaced
        if strcmpi(resp_append,'Replace')
          l_RoiDelete([],[],'all_dont_ask')
        end

        % If the ROI-structure contains only one slice, do not check dimensions
        if length(tmp.ROI(1).voxels)==1 && ...
            length(size(tmp.ROI(1).voxels{1}))==2 && ...
            ( length(size(DATA{1}.FTDATA))>2 | length(DATA)>1 )

          if Dat.isDataMixed
            % Resize first ROI slice
            if size(tmp.ROI(1).voxels{1})~=size(DATA{Dat.DataInd}.FTDATA)
              fprintf(1,'Warning: Resizing ROIs to current data size...\n')
            end
            roi_slice = {};
            for ii=1:length(tmp.ROI)
              roi_slice{ii} = imresize(tmp.ROI(ii).voxels{1},...
                size(DATA{Dat.DataInd}.FTDATA),...
                'nearest');
            end

            % Add empty masks into the ROI structure
            for ii=1:Dat.DataL
              if ii==Dat.DataInd
                for kk=1:length(tmp.ROI)
                  tmp.ROI(kk).voxels{ii} = roi_slice{kk};
                end
              else
                for kk=1:length(tmp.ROI)
                  tmp.ROI(kk).voxels{ii} = false(size(DATA{ii}.FTDATA));
                end
              end
            end
            clear('roi_slices')
          else
            if ~isequal(size(tmp.ROI(1).voxels{1}),...
                size(DATA{Dat.DataInd}.FTDATA(:,:,Dat.Slices(3),Dat.CurrentVol)))
              fprintf(1,'Warning: Resizing ROIs to current data size...\n')
            end
            roi_slices = {};
            for ii=1:length(tmp.ROI)
              roi_slices{ii} = imresize(tmp.ROI(ii).voxels{1},...
                size(DATA{1}.FTDATA(:,:,1,Dat.CurrentVol)),...
                'nearest');
            end

            % Add empty masks into the ROI structure
            for kk=1:length(tmp.ROI)
              tmp.ROI(kk).voxels{1} = false(size(DATA{1}.FTDATA));
              tmp.ROI(kk).voxels{1}(:,:,Dat.Slices(3),Dat.CurrentVol) = roi_slices{kk};
            end
            clear('roi_slices')
          end


        else %% ----------------------------------------------------------------

          % Check that ROI(s) are valid in size
          equalFiles = true;
          resizeNeeded = false;
          isRotationNeeded = false;
          for ii=1:Dat.DataL
            roi_sz=[size(tmp.ROI(1).voxels{ii},1),...
							size(tmp.ROI(1).voxels{ii},2),...
							size(tmp.ROI(1).voxels{ii},3),...
							size(tmp.ROI(1).voxels{ii},4)];
            im_sz = [size(DATA{ii}.FTDATA,1),...
							size(DATA{ii}.FTDATA,2),...
							size(DATA{ii}.FTDATA,3),...
							size(DATA{ii}.FTDATA,4)];
            
            % Allow loading a 3D ROI over a 4D data if 3D dimensions
            % match...
            if ~Dat.isDataMixed && ndims(DATA{ii}.FTDATA)==4 && ...
                ndims(tmp.ROI(1).voxels{ii})==3 && ...
                ismember(roi_sz(1:3),im_sz(1:3),'rows')
              
              tmp2=tmp.ROI;
              for kk=1:length(tmp2)
                tmp.ROI(kk).voxels{ii} = false(size(DATA{ii}.FTDATA));
                tmp.ROI(kk).voxels{ii}(:,:,:,Dat.CurrentVol)=tmp2(kk).voxels{ii};
              end
              %roi_sz=size(tmp.ROI(1).voxels{ii});
              %im_sz = size(DATA{ii}.FTDATA);
							roi_sz(4) = im_sz(4);
              clear tmp2
            end
            
            %% Abort loading of ROIs if relative dimensions are not the same
            if length(tmp.ROI(1).voxels)~=length(DATA) || ...
                ( length(roi_sz)~=length(im_sz) ) || ...
								~ismember(roi_sz,im_sz,'rows')
							if all(roi_sz(3:4)==1)
								resizeNeeded = true;
							else
								h=errordlg('Cannot continue loading! ROI and DATA sizes do not match!',...
									'ROI and DATA sizes do not match!','modal');
								clear tmp
								delete(hCalcWait)
								return
							end
%             elseif ~all(roi_sz==im_sz)
%               %% Warn if resizing of ROIs is needed
%               resizeNeeded = true;
            end

            %if ( length(roi_sz)~=length(im_sz) ) || ...
            %      ~all(roi_sz==im_sz) || length(tmp.ROI(1).voxels)~=length(DATA)
            %  % Herja
            %  h=errordlg('Cannot continue loading! ROI and DATA sizes do not match!',...
            %             'ROI and DATA sizes do not match!','modal');
            %  clear tmp
            %  delete(hCalcWait)
            %  return
            %end
            if isfield(tmp,'FileInfo') && ~isempty(tmp.FileInfo.DataFileName{ii})
              if ~strcmp(tmp.FileInfo.DataFileName{ii},Dat.HDR{ii}.fname) || ...
                  ~strcmp(tmp.FileInfo.DataPathName{ii},Dat.HDR{ii}.fpath)
                equalFiles = false;
              end
            end
          end

          %% Warn about file names
          if ~equalFiles
            resp=questdlg({'The file name(s) in the ROI do not match with the current data file name(s)','',...
              'Continue loading ROI(s) anyway?'},...
              'File name mismatch!',...
              'Continue anyway','Abort','Continue anyway');
            if isempty(resp) || strcmpi(resp,'Abort')
              clear tmp
              delete(hCalcWait)
              return
            end
          end

          %% Resize ROIs if needed
          if resizeNeeded
            fprintf(1,'Warning: Resizing ROIs to current data size...\n')
            if Dat.isDataMixed
              for ii=1:Dat.DataL
                im_sz = size(DATA{ii}.FTDATA);
                for kk=1:length(tmp.ROI)
                  tmp.ROI(kk).voxels{ii}=imresize(tmp.ROI(kk).voxels{ii},im_sz,'nearest');
                end
              end
            else
              for ii=1:size(DATA{1}.FTDATA,3)
                im_sz = size(DATA{1}.FTDATA(:,:,ii));
                for kk=1:length(tmp.ROI)
                  tmp.ROI(kk).voxels{1}(:,:,ii)=imresize(tmp.ROI(kk).voxels{1}(:,:,ii),...
                    im_sz,'nearest');
                end
              end
            end
          end



          % Check data rotation and flipping
          if isfield(tmp,'RotateFlip')

            % Check if rotation or flipping is needed
            if isstruct(tmp.RotateFlip)
              rotNeeded = ~ismember(Dat.DataRotation,tmp.RotateFlip.Rotate,'rows','legacy');  % Added 'legacy' -flag  mjn 09/2018
              flipNeeded = ~ismember(Dat.DataFlip,tmp.RotateFlip.Flip,'rows','legacy');       % Added 'legacy' -flag  mjn 09/2018
              isRotationNeeded = (rotNeeded | flipNeeded);
            else
              if ~isempty(tmp.RotateFlip)
                isRotationNeeded = true;
              end
            end
          end

          if isRotationNeeded
            resp=questdlg({['The ROI file suggests that some images should ' ...
              'be rotated and/or flipped before loading the saved ROI(s).'],'',...
              ['Do you want to rotate/flip images according to the',...
              ' orientation information saved in the ROI file?']},...
              'Rotate/flip images before loading ROI(s)?',...
              'Yes','No','Abort',...
              'Abort');
            if isempty(resp) || strcmpi(resp,'Abort')
              clear tmp
              delete(hCalcWait)
              return
            elseif strcmpi(resp,'Yes')


              if isstruct(tmp.RotateFlip)
                % Reset images before rotating/flipping
                l_RotateFlip([],[],'reset')

                Dat.DataRotation = tmp.RotateFlip.Rotate;
                Dat.DataFlip = tmp.RotateFlip.Flip;
                for ii=1:Dat.DataL
                  if Dat.DataRotation(ii)~=0
                    DATA{ii}.FTDATA = rot90(DATA{ii}.FTDATA,Dat.DataRotation(ii));
                    if ~isempty(ROI)
                      for kk=1:length(ROI)
                        ROI(kk).voxels{ii} = rot90(ROI(kk).voxels{ii},Dat.DataRotation(ii));
                      end
                    end
                  end

                  if Dat.DataFlip(ii)~=0
                    if Dat.DataFlip(ii)==1
                      DATA{ii}.FTDATA = flipud(DATA{ii}.FTDATA);
                      if ~isempty(ROI)
                        for kk=1:length(ROI)
                          ROI(kk).voxels{ii} = flipud(ROI(kk).voxels{ii});
                        end
                      end
                    elseif Dat.DataFlip(ii)==2
                      DATA{ii}.FTDATA = fliplr(DATA{ii}.FTDATA);
                      if ~isempty(ROI)
                        for kk=1:length(ROI)
                          ROI(kk).voxels{ii} = fliplr(ROI(kk).voxels{ii});
                        end
                      end
                    end

                  end
                end

              else
                % Rotate 3d data

                % To be written ...

                clear tmp
                delete(hCalcWait)
                h=errordlg('This feature has not been implemented yet.',...
                  'Feature not implemented','modal');
                return
              end



            end
          end
        end % if length(tmp.

        % Resolve possible conflicts if the ROI(s) will be appended
        if strcmpi(resp_append,'Append')

          % First check that we don't break the 18 ROI limit
          if (length(ROI)+length(tmp.ROI))>18
            h=errordlg({'Could not load ROI file',['"',f_path,f_name,'"'],...
              'The maximum number of ROI(s) (18) is exceeded!'},...
              'Could not load ROI(s)','modal');
            delete(hCalcWait)
            delete tmp
            return
          end

          % Change ROI colors
          le=length(ROI);
          free_colors=Dat.RoiColors(find(ismember(Dat.RoiColors,...
            reshape([ROI(:).color],3,le)', ...
            'rows')==0),:);

          for ii=1:length(tmp.ROI)
            tmp.ROI(ii).color = free_colors(ii,:);
          end

          % Change conflicting labels
          ind_label = ismember({tmp.ROI(:).label},{ROI(:).label});
          if any(ind_label)
            for ii=find(ind_label)
              done = false;
              counter = 1;
              while ~done
                if counter==1
                  new_label = [tmp.ROI(ii).label,'_appended'];
                else
                  new_label = [tmp.ROI(ii).label,'_appended',num2str(counter)];
                end

                if ~any(strcmpi(new_label,{ROI(:).label}))
                  tmp.ROI(ii).label = new_label;
                  done = true;
                else
                  counter = counter+1;
                end
              end
            end
          end

          % Append ROI(s) to existing
          tmp2 = tmp;
          tmp.ROI=ROI;
          for ii=1:length(tmp2.ROI)
            tmp.ROI(end+1) = tmp2.ROI(ii);
          end

          % Clear now the existing ROI(s)
          l_RoiDelete([],[],'all_dont_ask')
        end

        %% Update file names and paths
        for ii=1:length(tmp.ROI)
          tmp.ROI(ii).fpath = {};
          tmp.ROI(ii).fname = {};
          for kk=1:length(DATA)
            tmp.ROI(ii).fpath{kk}=DATA{kk}.HDR.fpath;
            tmp.ROI(ii).fname{kk}=DATA{kk}.HDR.fname;
          end
        end

        % Set loaded ROI(s)
        ROI = tmp.ROI;
        clear('tmp')
			end

      for kk=1:length(ROI)
        % Create axes and image for the ROI...
        for ii=1:3
          H.ROIAX(ii,kk) = axes('parent',H.FIG,...
            'units',...
            eval(['get(H.IMAX' sprintf('%d',ii) ',''units'')']),...
            'Position',...
            eval(['get(H.IMAX' num2str(ii) ',''position'')']),...
            'Ylim',...
            eval(['get(H.IMAX' sprintf('%d',ii) ',''ylim'')']), ...
            'xlim',...
            eval(['get(H.IMAX' sprintf('%d',ii) ',''xlim'')']),...
            'Yticklabel',{},...
            'Xticklabel',{},...
            'ydir','reverse',...
            'Box','on', ...
            'visible','off',...
            'hittest','off',...
            'alimmode','manual');
          H.IMROI(ii,kk)=image('parent',H.ROIAX(ii,kk),...
            'cdata',[],'visible','off',...
            'hittest','off');
					if Dat.HG2graphics
						set(H.IMROI(ii,kk),'hittest','on',...
							'buttondownfcn',@l_SetMouseGestures)
          end
        end
      end

      % Set New ROI label to the View listbox and
      % Roi Edit popupmenu
      set(H.ROIVIEW_LBOX,'string',{ROI(:).label},...
        'value',1:length(ROI))
      set(H.ROI_EDIT,'string',{ROI(:).label},...
        'value',1)

      % Enable delete buttons and edit popup
      set(H.UICH_ROIENABLED,'enable','on')

      % Set viewed ROI(s)
      Dat.RoiView = get(H.ROIVIEW_LBOX,'value');

      % Refresh ROIs
      l_RefreshRoi;
      
      % Check if aedes_calc_wait window still exists
      if ~isempty(hCalcWait) && ishandle(hCalcWait)
        delete(hCalcWait)
      end
      
      % Draw crossbars
      l_ShowHideCrossbars([],[])

      % Update ROI Info
      l_DisplayData([],[])


    catch
      aedes_errordump(lasterror);
    end
  end % function l_RoiLoad()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% View/Export ROI Statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_RoiViewStats(h,evd)
try
  
% Calculate results
Res = aedes_roi_stats(DATA,ROI);
if isempty(Res)
  h=errordlg('Unknown error: Statistics calculation failed!',...
             'Unknown error','modal');
  return
end

% Open results in aedes_resviewer
aedes_resviewer(Res)

catch
  aedes_errordump(lasterror);
end
end % function l_RoiViewStats(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotate/Flip ROI(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_RoiFlip(h,evd,opt)
try
if isempty(ROI)
  return
end


% Get current ROI
currentRoi = get(H.ROI_EDIT,'value');

% Add current ROI to undo buffer
l_RoiUndoBuffer('add_indices',currentRoi,[]);
  
% Loop over slices
if Dat.isDataMixed
  
  % Add current ROI to undo buffer
  l_RoiUndoBuffer('add_indices',currentRoi,[]);
  
  for ii=1:length(ROI(currentRoi).voxels)
    if strcmpi(opt,'flipud')
      % FlipUD current ROI
      ROI(currentRoi).voxels{ii} = flipud(ROI(currentRoi).voxels{ii});
    elseif strcmpi(opt,'fliplr')
      % FlipLR current ROI
      ROI(currentRoi).voxels{ii} = fliplr(ROI(currentRoi).voxels{ii});
    end
  end
else
  
  % Prompt for flip direction
  if Dat.AxView==0
    resp = questdlg('Select fliping plane', ...
      'Flipping plane?', ...
      'X', 'Y', 'Z', 'X');
    if isempty(resp)
      % Cancelled
      return
    end
    FlipPlane = resp;
  elseif Dat.AxView==1
    if strcmpi(opt,'flipud')
      FlipPlane='Y';
    elseif strcmpi(opt,'fliplr')
      FlipPlane='X';
    end
  elseif Dat.AxView==2
    FlipPlane='Y';
  else
    FlipPlane='Z';
  end
  
  if FlipPlane=='X'
    if strcmpi(opt,'flipud')
      ROI(currentRoi).voxels{1} = flipdim(ROI(currentRoi).voxels{1},1);
    elseif strcmpi(opt,'fliplr')
      ROI(currentRoi).voxels{1} = flipdim(ROI(currentRoi).voxels{1},2);
    end
  elseif FlipPlane=='Y'
    if strcmpi(opt,'flipud')
      ROI(currentRoi).voxels{1} = flipdim(ROI(currentRoi).voxels{1},1);
    elseif strcmpi(opt,'fliplr')
      ROI(currentRoi).voxels{1} = flipdim(ROI(currentRoi).voxels{1},3);
    end
  else
    if strcmpi(opt,'flipud')
      ROI(currentRoi).voxels{1} = flipdim(ROI(currentRoi).voxels{1},2);
    elseif strcmpi(opt,'fliplr')
      ROI(currentRoi).voxels{1} = flipdim(ROI(currentRoi).voxels{1},3);
    end
  end
end

% Refresh ROI(s)
l_RefreshRoi;

return

catch
  aedes_errordump(lasterror);
end
end % function l_RoiRotateFlip(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate boolean operations for ROIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_RoiBooleanOperations(h,evd)
try
if length(ROI)<2
  h=errordlg(['At least two (2) ROIs have to be defined in order to compare ',...
              'ROIs'],'Error','modal');
  return
end

cancel = true;

%% Load default font and colors
GD=aedes_gui_defaults;

% ROI Labels
roi_labels = {ROI(:).label};

% Overlay controls figure
fig_h = 220;
fig_w = 245;
fig_location = aedes_dialoglocation([fig_w,fig_h]);
fig_pos = [fig_location(1) fig_location(2) fig_w fig_h];

H.ROIBOOLEAN_FIG = figure('units','pixel',...
                             'position',...
                             fig_pos,...
                             'Name','ROI Boolean Operations',...
                             'numbertitle','off',...
                             'Toolbar','none',...
                             'Menubar','none',...
                             'DockControls','off',...
                             'renderer','painters',...
                             'resize','off',...
                             'CloseRequestFcn','delete(gcbf)',...
                             'Handlevisibility','off',...
                             'windowstyle','modal');
if ~GD.HG2graphics
	set(H.ROIBOOLEAN_FIG,'DoubleBuffer','on')
end
uipanel_h = uipanel('parent',H.ROIBOOLEAN_FIG,...
                    'units','pixel',...
                    'position',[5 40 fig_w-10 fig_h-45]);
set(H.ROIBOOLEAN_FIG,'color',get(uipanel_h,'backgroundcolor'))
% ROI A
roia_tx = uicontrol('parent',uipanel_h,...
                    'units','pixel',...
                    'position',[10 140 150 15],...
                    'style','text',...
                    'string','ROI "A"',...
                    'horizontalalign','left');
tmp = get(roia_tx,'position');
roia_popup = uicontrol('parent',uipanel_h,...
                       'units','pixel',...
                       'position',[tmp(1) tmp(2)-20 200 20],...
                       'style','popup',...
                       'string',roi_labels,...
                       'value',1,...
                      'callback',['if get(get(gcbo,''userdata''),''value'')==get(gcbo,''value''),',...
                    'if get(gcbo,''value'')==1,set(get(gcbo,''userdata''),''value'',get(gcbo,''value'')+1),',...
                    'elseif get(gcbo,''value'')==length(get(gcbo,''string'')),',...
                    'set(get(gcbo,''userdata''),''value'',get(gcbo,''value'')-1),',...
                    'else,set(get(gcbo,''userdata''),''value'',get(gcbo,''value'')-1),end,end'],...
                       'backgroundcolor','w');
tmp = get(roia_popup,'position');

% ROI B
roib_tx = uicontrol('parent',uipanel_h,...
                    'units','pixel',...
                    'position',[tmp(1) tmp(2)-30 150 15],...
                    'style','text',...
                    'string','ROI "B"',...
                    'horizontalalign','left');
tmp = get(roib_tx,'position');
roib_popup = uicontrol('parent',uipanel_h,...
                       'units','pixel',...
                       'position',[tmp(1) tmp(2)-20 200 20],...
                       'style','popup',...
                       'string',roi_labels,...
                       'value',2,...
                       'callback',['if get(get(gcbo,''userdata''),''value'')==get(gcbo,''value''),',...
                    'if get(gcbo,''value'')==1,set(get(gcbo,''userdata''),''value'',get(gcbo,''value'')+1),',...
                    'elseif get(gcbo,''value'')==length(get(gcbo,''string'')),',...
                    'set(get(gcbo,''userdata''),''value'',get(gcbo,''value'')-1),',...
                    'else,set(get(gcbo,''userdata''),''value'',get(gcbo,''value'')-1),end,end'],...
                       'backgroundcolor','w');
tmp = get(roib_popup,'position');
set(roia_popup,'userdata',roib_popup)
set(roib_popup,'userdata',roia_popup)

% Comparisons
comp_tx = uicontrol('parent',uipanel_h,...
                    'units','pixel',...
                    'position',[tmp(1) tmp(2)-30 150 15],...
                    'style','text',...
                    'string','Boolean Operations',...
                    'horizontalalign','left');
tmp = get(comp_tx,'position');
comp_popup = uicontrol('parent',uipanel_h,...
                       'units','pixel',...
                       'position',[tmp(1) tmp(2)-20 200 20],...
                       'style','popup',...
                       'value',1,...
                       'string',{'Union (either A or B)',...
                    'Intersect (mutual A and B)',...
                    'Subtraction (A - B)',...
                    'Exclusive disjunction (A xor B)'},...
                       'backgroundcolor','w');

% OK and Cancel buttons
tmp = get(uipanel_h,'position');
cancel_btn = uicontrol('parent',H.ROIBOOLEAN_FIG,...
                       'position',[tmp(1)+tmp(3)-65 5 65 30],...
                       'string','Cancel',...
                       'callback','delete(gcbf)');
tmp = get(cancel_btn,'position');
ok_btn = uicontrol('parent',H.ROIBOOLEAN_FIG,...
                   'position',[tmp(1)-65-5 tmp(2) tmp(3) tmp(4)],...
                   'string','OK',...
                   'callback','uiresume(gcbf)');

% Wait for exit
uiwait(H.ROIBOOLEAN_FIG)
if ishandle(H.ROIBOOLEAN_FIG)
  
  % Get values from the popup menus
  roi_a = get(roia_popup,'value');
  roi_b = get(roib_popup,'value');
  comp_ind = get(comp_popup,'value');
  
  % Close figure
  delete(H.ROIBOOLEAN_FIG)
  
  switch comp_ind
   case 1
    str = 'union';
   case 2
    str = 'intersect';
   case 3
    str = 'subtract';
   case 4
    str = 'xor';
  end
  
  l_RoiAdd([],[],str,[roi_a,roi_b])
end
return

catch
  aedes_errordump(lasterror);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Rename current ROI
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_RoiRename(h,evd)
try
  
  if isempty(ROI)
    return
  end
  
  % Get current ROI
  roi_ind=get(H.ROI_EDIT,'value');
  roilabel = ROI(roi_ind).label;
  
  % Ask for a new label
  done=false;
  while ~done
    new_label=aedes_inputdlg('Give a new label for the ROI','New ROI Label?', ...
                           roilabel);
    if isempty(new_label)
      return
    else
      new_label = new_label{1};
    end
    new_label = strtrim(new_label);
    if isempty(new_label)
      return
    end
    
    % Check if ROI with the same label already exists
    roi_labels = {ROI(:).label};
    if any(strcmp(new_label,roi_labels))
      h=warndlg({['A ROI with the label "' new_label '" already exists.'],...
                 'Choose another label.'},'Label already exists','modal');
      uiwait(h)
    else
      done=true;
    end
  end
  
  % Apply new label to the ROI
  ROI(roi_ind).label = new_label;
  
  % Apply changes to ROI listbox and popupmenu
  popupstr=get(H.ROI_EDIT,'string');
  popupval=get(H.ROI_EDIT,'value');
  if iscell(popupstr)
    popupstr{roi_ind}=new_label;
  else
    popupstr=new_label;
  end
  set(H.ROI_EDIT,'string',popupstr,...
                 'value',popupval)
  
  lboxstr = get(H.ROIVIEW_LBOX,'string');
  lboxval = get(H.ROIVIEW_LBOX,'value');
  if iscell(lboxstr)
    lboxstr{roi_ind}=new_label;
  else
    lboxstr=new_label;
  end
  set(H.ROIVIEW_LBOX,'string',lboxstr,...
                    'value',lboxval)
  
catch
  aedes_errordump(lasterror);
end
end % function l_RoiRename(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Show/Hide ROI edges
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_UiRoiEdgesCB(h,evd)
try
  
if strcmpi(get(h,'checked'),'off')
  set(h,'checked','on')
  l_ShowRoiEdges([],[],'enable')
else
  set(h,'checked','off')
  l_ShowRoiEdges([],[],'disable')
end

catch
  aedes_errordump(lasterror);
end
end % function l_UiRoiEdgesCB(h,
function l_ShowRoiEdges(h,evd,opt)
try
  
if strcmpi(opt,'enable')
  Dat.ShowRoiEdges = true;
  l_ShowRoiEdges([],[],'all')
  return
elseif strcmpi(opt,'disable')
  if Dat.isDataMixed
    for ii=1:length(ROI)
      OldLines=findobj(H.ROIAX(1,ii),'tag','roi_edgeline');
      if ~isempty(OldLines)
        delete(OldLines);
      end
    end
  else
    for kk=1:3
      for ii=1:length(ROI)
        OldLines=findobj(H.ROIAX(kk,ii),'tag','roi_edgeline');
        if ~isempty(OldLines)
          delete(OldLines);
        end
      end
    end
  end
  Dat.ShowRoiEdges = false;
  return
end

if isempty(ROI) || ~Dat.ShowRoiEdges %|| ~Dat.isDataMixed
  return
end

% Update all
if strcmpi(opt,'all')
  UpdatedRois=1:length(ROI);
else
  % Get current ROI
  CurrentRoi = get(H.ROI_EDIT,'value');
  UpdatedRois=CurrentRoi;
end

% Get Orientations to update
if Dat.AxView==0
  update_ind = [1 2 3];
else
  update_ind = Dat.AxView;
end

ind = [3 2 1];

for ii=update_ind
  for kk=UpdatedRois
    
    % Check if ROI is empty
    if ii==3
      CurrentSlice = ROI(kk).voxels{Dat.DataInd}(:,:,Dat.Slices(3),Dat.CurrentVol);
    elseif ii==2
      CurrentSlice = reshape(ROI(kk).voxels{Dat.DataInd}(:,Dat.Slices(2),:,Dat.CurrentVol),...
                             Dat.ImageDim(Dat.DataInd,1),Dat.ImageDim(Dat.DataInd,3));
    else
      CurrentSlice = reshape(ROI(kk).voxels{Dat.DataInd}(Dat.Slices(1),:,:,Dat.CurrentVol),...
                             Dat.ImageDim(Dat.DataInd,2),Dat.ImageDim(Dat.DataInd,3));
    end
    
    tmp=find(CurrentSlice,1,'first');
    if isempty(tmp) || ~any(kk==Dat.RoiView)
      OldLines=findobj(H.ROIAX(ind(ii),kk),'tag','roi_edgeline');
      if ~isempty(OldLines)
        delete(OldLines);
      end
      continue
    end
    % Calculate ROI edges
    B=bwboundaries(CurrentSlice,4,'holes');
    
    
    if isempty(B)
      OldLines=findobj(H.ROIAX(ind(ii),kk),'tag','roi_edgeline');
      if ~isempty(OldLines)
        delete(OldLines);
      end
      continue
    end
    
    OldLines=findobj(H.ROIAX(ind(ii),kk),'tag','roi_edgeline');
    if ~isempty(OldLines)
      delete(OldLines);
    end
    
    for jj=1:length(B)
      boundary = B{jj};
      line('parent',H.ROIAX(ind(ii),kk),...
           'xdata',boundary(:,2),'ydata',boundary(:,1),...
           'color',ROI(kk).color./255,...
           'tag','roi_edgeline',...
           'linewidth',1.5,...
           'linestyle','-',...
           'hittest','off');
    end
  end
end
drawnow
  
catch
  aedes_errordump(lasterror);
end
end % function l_ShowRoiEdges(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Set ROI color
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function l_RoiSetColor(h,evd)
    try
      if isempty(ROI)
        return
      end
      
      roi_ind = get(H.ROI_EDIT,'value');
      str = get(H.ROI_EDIT,'string');
      str = str{roi_ind};
      old_color = ROI(roi_ind).color./255;
      
      % Prompt for color
      C=uisetcolor(old_color,['ROI: ',str]);
      
      % Set the new color
      ROI(roi_ind).color = round(C*255);
      
      % Redraw ROIs
      l_RefreshRoi;
    catch
      aedes_errordump(lasterror);
    end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Draw ROI(s) on top of the slices
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_RefreshRoi(update_ind,roi_ind)
try

% Look for a quick return
if isempty(H.ROIAX)
  return
end


if nargin<2
  roi_ind = Dat.RoiView;
  %roi_ind=1:length(ROI);
end
ind = [3 2 1];
% Refresh ROI masks. Draw only the masks that
% are needed at the moment.

% 3D-view. This is the slowest case because all
% three of the ROI masks have to be possibly drawn...
if Dat.AxView==0
  if nargin<1
    update_ind=[1 2 3];
  end
  for ii=update_ind
    for kk=roi_ind
      if ii==3
        tmp = ROI(kk).voxels{Dat.DataInd}(:,:,Dat.Slices(3),Dat.CurrentVol);
      elseif ii==2
        tmp = reshape(ROI(kk).voxels{Dat.DataInd}(:,Dat.Slices(2),:,Dat.CurrentVol),...
                      Dat.ImageDim(Dat.DataInd,1),Dat.ImageDim(Dat.DataInd,3));
      else
        tmp = reshape(ROI(kk).voxels{Dat.DataInd}(Dat.Slices(1),:,:,Dat.CurrentVol),...
                      Dat.ImageDim(Dat.DataInd,2),Dat.ImageDim(Dat.DataInd,3));
      end
      if isempty(find(tmp,1,'first'))
        set(H.IMROI(ind(ii),kk),'visible','off')
        continue
      end
      tmp=uint8(tmp);
      
      cdata = zeros([size(tmp) 3],'uint8');
      if ROI(kk).color(1)~=0
        cdata(:,:,1) = tmp*ROI(kk).color(1);
      end
      if ROI(kk).color(2)~=0
        cdata(:,:,2) = tmp*ROI(kk).color(2);
      end
      if ROI(kk).color(3)~=0
        cdata(:,:,3) = tmp*ROI(kk).color(3);
      end
      alphadata=double(tmp)*Dat.RoiTransp;
      
      set(H.IMROI(ind(ii),kk),'cdata',cdata,...
                        'AlphaData',alphadata,'visible','on')
      %if ~any(any(tmp))
      %  set(H.IMROI(ii,kk),'visible','off')
      %  continue
      %end
    end
  end
  
% Single slice view (X, Y, or Z). Draw only the ROI(s)
% visible in this view...
else
  for kk=roi_ind
    if Dat.AxView==3
      tmp = ROI(kk).voxels{Dat.DataInd}(:,:,Dat.Slices(3),Dat.CurrentVol);
    elseif Dat.AxView==2
      tmp = reshape(ROI(kk).voxels{Dat.DataInd}(:,Dat.Slices(2),:,Dat.CurrentVol),...
                    Dat.ImageDim(Dat.DataInd,1),Dat.ImageDim(Dat.DataInd,3));
    else
      tmp = reshape(ROI(kk).voxels{Dat.DataInd}(Dat.Slices(1),:,:,Dat.CurrentVol),...
                    Dat.ImageDim(Dat.DataInd,2),Dat.ImageDim(Dat.DataInd,3));
    end
    if isempty(find(tmp,1,'first'))
      set(H.IMROI([1 2 3]~=ind(Dat.AxView),kk),'visible','off')
      set(H.IMROI(ind(Dat.AxView),kk),'visible','off')
      continue
    end
    tmp=uint8(tmp);
    
    cdata = zeros([size(tmp) 3],'uint8');
    if ROI(kk).color(1)~=0
      cdata(:,:,1) = tmp*ROI(kk).color(1);
    end
    if ROI(kk).color(2)~=0
      cdata(:,:,2) = tmp*ROI(kk).color(2);
    end
    if ROI(kk).color(3)~=0
      cdata(:,:,3) = tmp*ROI(kk).color(3);
    end
    alphadata=double(tmp)*Dat.RoiTransp;
    
    set(H.IMROI([1 2 3]~=ind(Dat.AxView),kk),'visible','off')
    set(H.IMROI(ind(Dat.AxView),kk),'cdata',cdata,...
                      'AlphaDataMapping','none',...
                      'AlphaData',alphadata,...
                      'visible','on')
    
  end
end

catch
  aedes_errordump(lasterror);
end
end % function l_RefreshRoi(H,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Set transparencies for ROI(s)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_SetRoiTransparency(h,evd,opt)
try

	if ~isempty(h) && h==H.ROI_TRANSP_SLIDER
		opt = 'slider';
	end
switch opt
 case 'slider'
  s_val = get(H.ROI_TRANSP_SLIDER,'value');
  s_val = fix(s_val*100)/100;
  set(H.ROI_TRANSP_EDIT,'string',sprintf('%.2f',s_val))
  Dat.RoiTransp = s_val;
  set(H.ROI_TRANSP_EDIT,'userdata',Dat.RoiTransp)
  
 case 'edit'
  val=str2num(get(H.ROI_TRANSP_EDIT,'string'));
  % Check input values
  if ~isreal(val) || isempty(val) || val<0 || val>1
    set(H.ROI_TRANSP_EDIT,'string',...
                      num2str(get(H.ROI_TRANSP_EDIT,'userdata')))
  else
    % Round value to two decimals
    Dat.RoiTransp = fix(val*100)/100;
    set(H.ROI_TRANSP_EDIT,'userdata',Dat.RoiTransp)
  end
  set(H.ROI_TRANSP_SLIDER,'value',Dat.RoiTransp)
  
end

% Update preferences
setpref('Aedes','RoiTransp',Dat.RoiTransp)

% Update ROI(s)
l_RefreshRoi;

catch
  aedes_errordump(lasterror);
end
end % function l_SetRoiTransparency(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Set mouse to control either ROI
%% functions or viewing
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_SetMouseGestures(h,evd,opt)
try
% Check which mouse button was pressed
seltype = find(strcmpi(get(H.FIG,'selectiontype'),{'normal','extend', ...
                    'alt','open'}));

% Get modifier button
CurrentModifier=get(H.FIG,'CurrentModifier');

%% Draw ROI
if strcmpi(get(H.UITOGGLE_DRAWROI,'state'),'on') || ...
      strcmpi(get(H.UITOGGLE_ERASEROI,'state'),'on') || ...
			strcmpi(get(H.UITOGGLE_FILLROI,'state'),'on')
  
  if isempty(ROI)
    return
  end
  
  % ROI changed
  Dat.RoiSaved=0;
  
  % Check current ROI
  roi_ind = get(H.ROI_EDIT,'value');
  str = get(H.ROI_EDIT,'string');
  str = str(roi_ind);
  if isempty(str)
    return
	end
	if any(h==H.IM)
		ax_ind = find(H.IM==h);
	elseif any(h==H.IMROI(:))
		[ax_ind,J] = find(H.IMROI==h);
	elseif any(h==H.IMOVERLAY)
		ax_ind = find(H.IMOVERLAY==h);
	else
		ax_ind = find(H.IM==get(h,'userdata'));
	end
  
  %%%%%%%%%%%%%%%%%%%%%%%
  % Left mouse button
  %%%%%%%%%%%%%%%%%%%%%%%
  if seltype==1
    
    if length(CurrentModifier)==2
      if all(ismember({'shift','control'},CurrentModifier))
				% Add undo information to buffer
				l_RoiUndoBuffer('add',roi_ind,ax_ind);
				
        l_RoiMove([],[],'move2d',false,ax_ind,roi_ind)
        set(H.FIG,'windowbuttonmotionfcn',...
                  {@l_RoiMove,'move2d',true,ax_ind,roi_ind},...
                  'windowbuttonupfcn',...
                  @l_WindowButtonUpFcn)
      else
        return
      end
    else
      % Draw ROI -----------------------
      if strcmpi(get(H.UITOGGLE_DRAWROI,'state'),'on')
				% Add undo information to buffer
				l_RoiUndoBuffer('add',roi_ind,ax_ind);

        set(H.FIG,'windowbuttonmotionfcn',...
                  {@l_RoiDrawErase,'draw',true,ax_ind,roi_ind},...
                  'windowbuttonupfcn',...
                  @l_WindowButtonUpFcn)
        l_RoiDrawErase([],[],'draw',false,ax_ind,roi_ind)
        
        % Erase ROI ----------------------
      elseif strcmpi(get(H.UITOGGLE_ERASEROI,'state'),'on')
				% Add undo information to buffer
				l_RoiUndoBuffer('add',roi_ind,ax_ind);
				
        set(H.FIG,'windowbuttonmotionfcn',...
                  {@l_RoiDrawErase,'erase',true,ax_ind,roi_ind},...
                  'windowbuttonupfcn',...
                  @l_WindowButtonUpFcn)
        l_RoiDrawErase([],[],'erase',false,ax_ind,roi_ind)

				% Fill ROI -----------------------
			elseif strcmpi(get(H.UITOGGLE_FILLROI,'state'),'on')
				l_RoiFill([],[],'fill',ax_ind,roi_ind)
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Middle mouse button
    %%%%%%%%%%%%%%%%%%%%%%%%
  elseif seltype==2
    if all(ismember({'shift','control'},CurrentModifier))
      % Drag view
      cp = get(H.FIG,'currentpoint');
      
      % Get image slider values
      val1 = get(H.IMSLIDER(1),'value');
      val2 = get(H.IMSLIDER(2),'value');
      setptr(H.FIG,'closedhand');
      set(H.FIG,'windowbuttonmotionfcn',...
                {@l_DragView,cp,val1,val2},...
                'windowbuttonupfcn',...
                @l_WindowButtonUpFcn)

		else
      
      % Fill ROI
      if strcmpi(get(H.UITOGGLE_DRAWROI,'state'),'on')
        l_RoiFill([],[],'fill',ax_ind,roi_ind)
        
        % Erase Fill ROI
      elseif strcmpi(get(H.UITOGGLE_ERASEROI,'state'),'on')
        l_RoiFill([],[],'fill_remove',ax_ind,roi_ind)
			elseif strcmpi(get(H.UITOGGLE_FILLROI,'state'),'on')
				set(H.FIG,'windowbuttonmotionfcn',...
                  {@l_RoiDrawErase,'draw',true,ax_ind,roi_ind},...
                  'windowbuttonupfcn',...
                  @l_WindowButtonUpFcn)
        l_RoiDrawErase([],[],'draw',false,ax_ind,roi_ind)
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Right mouse button
    %%%%%%%%%%%%%%%%%%%%%%%%%
  elseif seltype==3
    if length(CurrentModifier)==2
      if all(ismember({'shift','control'},CurrentModifier))
				 % Add undo information to buffer
				 l_RoiUndoBuffer('add',roi_ind,ax_ind);
				 
        if strcmpi(get(H.UITOGGLE_DRAWROI,'state'),'on')
          l_RoiFill([],[],'fill_remove',ax_ind,roi_ind)
        elseif strcmpi(get(H.UITOGGLE_ERASEROI,'state'),'on')
          l_RoiFill([],[],'fill',ax_ind,roi_ind)
        end
      else
        return
      end
    else
      if strcmpi(get(H.UITOGGLE_DRAWROI,'state'),'on')
				 % Add undo information to buffer
				 l_RoiUndoBuffer('add',roi_ind,ax_ind);
				 
        set(H.FIG,'windowbuttonmotionfcn',...
                  {@l_RoiDrawErase,'erase',true,ax_ind,roi_ind},...
                  'windowbuttonupfcn',...
                  @l_WindowButtonUpFcn)
        l_RoiDrawErase([],[],'erase',false,ax_ind,roi_ind)
      elseif strcmpi(get(H.UITOGGLE_ERASEROI,'state'),'on')
				 % Add undo information to buffer
				 l_RoiUndoBuffer('add',roi_ind,ax_ind);
				 
        set(H.FIG,'windowbuttonmotionfcn',...
                  {@l_RoiDrawErase,'draw',true,ax_ind,roi_ind},...
                  'windowbuttonupfcn',...
                  @l_WindowButtonUpFcn)
        l_RoiDrawErase([],[],'draw',false,ax_ind,roi_ind)
				elseif strcmpi(get(H.UITOGGLE_FILLROI,'state'),'on')
					l_RoiFill([],[],'fill_remove',ax_ind,roi_ind)
      end
    end
  end
  
elseif strcmpi(get(H.UITOGGLE_ARROW,'state'),'on') %% Slice/window/level viewing
  if seltype==1
    % Check which image was pressed
		if any(h==H.IM)
			ax_ind = find(H.IM==h);
		elseif any(h==H.IMROI(:))
			[ax_ind,J] = find(H.IMROI==h);
		elseif any(h==H.IMOVERLAY)
			ax_ind = find(H.IMOVERLAY==h);
		else
			ax_ind = find(H.IM==get(h,'userdata'));
		end
		ind = [3 2 1];
    l_ChangeSlice([],[],'mouse',ind(ax_ind))
    set(H.FIG,'windowbuttonmotionfcn',...
              {@l_ChangeSlice,'mouse',ind(ax_ind)},...
              'windowbuttonupfcn',...
              @l_WindowButtonUpFcn)
    
  elseif seltype==2
    % Drag view
    cp = get(H.FIG,'currentpoint');
    
    % Get image slider values
    val1 = get(H.IMSLIDER(1),'value');
    val2 = get(H.IMSLIDER(2),'value');
    setptr(H.FIG,'closedhand');
    set(H.FIG,'windowbuttonmotionfcn',...
              {@l_DragView,cp,val1,val2},...
              'windowbuttonupfcn',...
              @l_WindowButtonUpFcn)
    
  elseif seltype==3
    
    % Load mouse custom pointer from cdata.mat
    try
      tmp = H.btn_cdata;
      set(H.FIG,'pointer','custom',...
                'PointerShapeCData',tmp.cdata.contrastPointer.PointerShapeCData,...
                'PointerShapeHotSpot',tmp.cdata.contrastPointer.PointerShapeHotSpot)
    catch
      set(H.FIG,'pointer','arrow')
    end
    
    % Use the mouse to set window and level
    cp = get(H.FIG,'currentpoint');
    Dat.tmp_point = cp;
    
    set(H.FIG,'windowbuttonmotionfcn',...
              {@l_SetWindowLevel,[],1},...
              'windowbuttonupfcn',...
              @l_WindowButtonUpFcn)
  end
  
elseif strcmpi(get(H.UITOGGLE_PAN,'state'),'on')
  % Drag view
  cp = get(H.FIG,'currentpoint');
  
  % Get image slider values
  val1 = get(H.IMSLIDER(1),'value');
  val2 = get(H.IMSLIDER(2),'value');
  setptr(H.FIG,'closedhand');
  set(H.FIG,'windowbuttonmotionfcn',...
            {@l_DragView,cp,val1,val2},...
            'windowbuttonupfcn',...
            @l_WindowButtonUpFcn)
elseif strcmpi(get(H.UITOGGLE_MOVEROI,'state'),'on') || ...
    strcmpi(get(H.UITOGGLE_MOVEROI3D,'state'),'on')
  
  if isempty(ROI)
    return
  end
  
  % Check if we want to move in 2D or 3D
  if strcmpi(get(H.UITOGGLE_MOVEROI3D,'state'),'on')
    isMove3D = true;
  else
    isMove3D = false;
  end
  
  % ROI changed
  Dat.RoiSaved=0;
  
  % Check current ROI
  roi_ind = get(H.ROI_EDIT,'value');
  str = get(H.ROI_EDIT,'string');
  str = str(roi_ind);
  if isempty(str)
    return
	end
	if any(h==H.IM)
		ax_ind = find(H.IM==h);
	elseif any(h==H.IMROI(:))
		[ax_ind,J] = find(H.IMROI==h);
	elseif any(h==H.IMOVERLAY)
		ax_ind = find(H.IMOVERLAY==h);
	else
		ax_ind = find(H.IM==get(h,'userdata'));
	end
  
  % Add undo information to buffer
  if isMove3D
    l_RoiUndoBuffer('add_indices',roi_ind,[]);
    
     % Move ROI
     l_RoiMove([],[],'move3d',false,ax_ind,roi_ind)
     set(H.FIG,'windowbuttonmotionfcn',...
               {@l_RoiMove,'move3d',true,ax_ind,roi_ind},...
               'windowbuttonupfcn',...
               @l_WindowButtonUpFcn)
  else
    l_RoiUndoBuffer('add',roi_ind,ax_ind);
    
     % Move ROI
     l_RoiMove([],[],'move2d',false,ax_ind,roi_ind)
     set(H.FIG,'windowbuttonmotionfcn',...
               {@l_RoiMove,'move2d',true,ax_ind,roi_ind},...
               'windowbuttonupfcn',...
               @l_WindowButtonUpFcn)
  end
end

catch
  aedes_errordump(lasterror);
end
end % function l_SetMouseGestures(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Mouse wheel callback
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function l_MouseWheelFcn(h,evd)

    % Get the point where mouse cursor was when scroll wheel event was
		% fired
		p = get(0,'PointerLocation');
		cp = get(H.FIG,'Currentpoint');
		fig_pos = get(H.FIG,'position');
		x = (p(1)-fig_pos(1))/fig_pos(3);
		y = (p(2)-fig_pos(2))/fig_pos(4);
		pos = [232 0 fig_pos(3) fig_pos(4)-232];
		coord = [cp(1) cp(2)];
		coord(1) = coord(1)-pos(1);
		coord(2) = abs(coord(2)-pos(2)-pos(4));
		if coord(1) < 0 || coord(2) < 0 || ...
				coord(1) > pos(3) || coord(2) > pos(4)
			return
		end
		
		% Get axes absolute positions
		dim = 0;
		ax3_pos = get(H.IMAX1,'position');
		ax3_pos(1:2) = ax3_pos(1:2)+fig_pos(1:2);
		ax2_pos = get(H.IMAX2,'position');
		ax2_pos(1:2) = ax2_pos(1:2)+fig_pos(1:2);
		ax1_pos = get(H.IMAX3,'position');
		ax1_pos(1:2) = ax1_pos(1:2)+fig_pos(1:2);
		
		% Determine on which axes the pointer is currently on
		if Dat.AxView == 0 || Dat.AxView == 3
			if p(1)>=ax3_pos(1) && p(2)>=ax3_pos(2) && ...
					p(1)<=(ax3_pos(1)+ax3_pos(3)) && p(2)<=(ax3_pos(2)+ax3_pos(4))
				dim = 3;
			end
		end
		if Dat.AxView == 0 || Dat.AxView == 2
			if p(1)>=ax2_pos(1) && p(2)>=ax2_pos(2) && ...
					p(1)<=(ax2_pos(1)+ax2_pos(3)) && p(2)<=(ax2_pos(2)+ax2_pos(4))
				dim = 2;
			end
		end
		if Dat.AxView == 0 || Dat.AxView == 1
			if p(1)>=ax1_pos(1) && p(2)>=ax1_pos(2) && ...
					p(1)<=(ax1_pos(1)+ax1_pos(3)) && p(2)<=(ax1_pos(2)+ax1_pos(4))
				dim = 1;
			end
		end

    %cp=get(H.FIG,'CurrentPoint')
		modifier = get(H.FIG,'CurrentModifier');
		
		if length(modifier)>1
			return
		end

		if isempty(modifier)
			if dim==0
				% Not over image...
				return
			end
			
			% Change slice
			l_ChangeSlice([],dim,'wheel',evd.VerticalScrollCount);
		elseif strcmpi(modifier,'control')
			% Change zoom level
			%if Dat.ZoomLevel==0
			%	l_Zoom([],[],'normalize')
			%end
			
			if evd.VerticalScrollCount<0
				l_Zoom([],[],'+')
			else
				l_Zoom([],[],'-')
			end
			l_AxesPositions([],[],Dat.AxView,Dat.ZoomLevel)
		elseif strcmpi(modifier,'shift')
			% Fast scroll. Change 3 slices at a time.
			if evd.VerticalScrollCount<0
				cnt = evd.VerticalScrollCount-3;
			elseif evd.VerticalScrollCount>0
				cnt = evd.VerticalScrollCount+3;
			else
				return
			end
			l_ChangeSlice([],dim,'wheel',cnt);
		end
		
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Add to ROI UNDO buffer
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_RoiUndoBuffer(opt,roi_ind,ax_ind)
try
% This local function is used to manage the undo
% information in ROI drawing.


switch opt
  %%%%%%%%%%%%%%%%%%%%%
  % Add to undo buffer
  %%%%%%%%%%%%%%%%%%%%%
 case 'add'
  % The Undo buffer retains a copy of 10 last modified slices. The reason
  % that it does store indices (which would use in many cases less
  % memory) is that finding the indices using FIND command is relatively
  % slow (can be even few seconds with large datas). Nevertheless,
  % indices are stored when multiple slices are altered simultaneously
  % (e.g. when copying ROI).
  
  le=length(Dat.RoiUndo);
  if ax_ind==1
    Dat.RoiUndo(le+1).voxels = ROI(roi_ind).voxels{Dat.DataInd}(:,:,Dat.Slices(3),Dat.CurrentVol);
    Dat.RoiUndo(le+1).ax_ind = ax_ind;
    Dat.RoiUndo(le+1).label = ROI(roi_ind).label;
    Dat.RoiUndo(le+1).slice = Dat.Slices(3);
    Dat.RoiUndo(le+1).vol = Dat.CurrentVol;
    Dat.RoiUndo(le+1).DataInd = Dat.DataInd;
  elseif ax_ind==2
    Dat.RoiUndo(le+1).voxels = reshape(ROI(roi_ind).voxels{Dat.DataInd}(:,Dat.Slices(2),:,Dat.CurrentVol),...
                                       Dat.ImageDim(Dat.DataInd,1),Dat.ImageDim(Dat.DataInd,3));
    Dat.RoiUndo(le+1).ax_ind = ax_ind;
    Dat.RoiUndo(le+1).label = ROI(roi_ind).label;
    Dat.RoiUndo(le+1).slice = Dat.Slices(2);
    Dat.RoiUndo(le+1).vol = Dat.CurrentVol;
    Dat.RoiUndo(le+1).DataInd = Dat.DataInd;
  else
    Dat.RoiUndo(le+1).voxels = reshape(ROI(roi_ind).voxels{Dat.DataInd}(Dat.Slices(1),:,:,Dat.CurrentVol),...
                                       Dat.ImageDim(Dat.DataInd,2),Dat.ImageDim(Dat.DataInd,3));
    Dat.RoiUndo(le+1).ax_ind = ax_ind;
    Dat.RoiUndo(le+1).label = ROI(roi_ind).label;
    Dat.RoiUndo(le+1).slice = Dat.Slices(1);
    Dat.RoiUndo(le+1).vol = Dat.CurrentVol;
    Dat.RoiUndo(le+1).DataInd = Dat.DataInd;
  end
  undo_length=length(Dat.RoiUndo);
  if undo_length==11
    Dat.RoiUndo(1)=[];
  end
  set(H.ROIBTN_UNDO,'enable','on','string', num2str(length(Dat.RoiUndo)))
  set(H.UIROITOOLS_UNDO,'enable','on','label',['Undo draw (' num2str(length(Dat.RoiUndo)) ')'])
  set(H.UIPUSH_UNDOROI,'enable','on')
  %Dat.RoiUndo.voxels = ROI(roi_ind).voxels;
  %Dat.RoiUndo.ind = roi_ind;
  
 case 'add_indices'
  % Store current ROI indices. Can be slow so don't hold your breath...
  
  le=length(Dat.RoiUndo);
  for ii=1:length(ROI(roi_ind).voxels)
    Dat.RoiUndo(le+1).voxels{ii} = ...
        find(ROI(roi_ind).voxels{ii}(:,:,:,Dat.CurrentVol));
  end
  Dat.RoiUndo(le+1).ax_ind = 0;
  Dat.RoiUndo(le+1).label = ROI(roi_ind).label;
  Dat.RoiUndo(le+1).slice = 0;
  Dat.RoiUndo(le+1).vol = Dat.CurrentVol;
  Dat.RoiUndo(le+1).DataInd = Dat.DataInd;
  
  undo_length=length(Dat.RoiUndo);
  if undo_length==11
    Dat.RoiUndo(1)=[];
  end
  set(H.ROIBTN_UNDO,...
      'enable','on','string',...
      num2str(length(Dat.RoiUndo)))
  set(H.UIROITOOLS_UNDO,...
      'enable','on','label',...
      ['Undo draw (' num2str(length(Dat.RoiUndo)) ')'])
  set(H.UIPUSH_UNDOROI,'enable','on')
  %Dat.RoiUndo.voxels = ROI(roi_ind).voxels;
  %Dat.RoiUndo.ind = roi_ind;

  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Remove from undo buffer
  %%%%%%%%%%%%%%%%%%%%%%%%%%
 case 'remove'
  
  % Correct ROI UNDO indices
  %ind=find(ismember({Dat.RoiUndo(:).label},{str{val}}));
  Dat.RoiUndo(roi_ind)=[];
  %setappdata(H.FIG,'Dat',Dat)
  if isempty(Dat.RoiUndo)
    set(H.ROIBTN_UNDO,'enable','off',...
                      'string','0')
    set(H.UIROITOOLS_UNDO,'enable','off',...
                      'label','Undo draw (0)')
    set(H.UIPUSH_UNDOROI,'enable','off')
  else
    set(H.ROIBTN_UNDO,'string',num2str(length(Dat.RoiUndo)))
    set(H.UIROITOOLS_UNDO,'label',['Undo draw (' num2str(length(Dat.RoiUndo)) ...
                        ')'])
    set(H.UIPUSH_UNDOROI,'enable','on')
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%
  % Clear the undo buffer
  %%%%%%%%%%%%%%%%%%%%%%%%
 case 'clear'
  
  Dat.RoiUndo = [];
  set(H.ROIBTN_UNDO,'enable','off','string','0')
  set(H.UIROITOOLS_UNDO,'enable','off',...
                      'label','Undo draw (0)')
  set(H.UIPUSH_UNDOROI,'enable','off')
end


catch
  aedes_errordump(lasterror);
end
end % function l_RoiUndoBuffer(opt,


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Change figure colormap
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_ChangeColormap(h,evd,opt)
try
  
if isempty(opt)
  str = lower(popupstr(H.COLMAP_POPUP));
else
  str = lower(opt);
  
  % Set popupmenu value
  tmp=get(H.COLMAP_POPUP,'string');
  ind=find(strcmpi(tmp,str));
  set(H.COLMAP_POPUP,'value',ind)
end

% Set new colormap
if strcmpi(str,'Custom')
  % Try to get the colormap from preferences
  if ispref('Aedes','CustomColorMap')
	% Check that the size of the custom colormap is ok
	custom_cmap = getpref('Aedes','CustomColorMap');
	if ~isequal(size(custom_cmap),[256 3]) || ...
		~isnumeric(custom_cmap)
	  h=errordlg(['The custom colormap in preferences is not valid.',...
		' It should be a 256x3 numeric matrix.'],...
		'Invalid custom colormap!','modal');
	  
	  % Set popupmenu value
	  if ispref('Aedes','ColorMap')
		str=getpref('Aedes','ColorMap');
	  else
		str = 'gray';
	  end
	  tmp=get(H.COLMAP_POPUP,'string');
	  ind=find(strcmpi(tmp,str));
	  set(H.COLMAP_POPUP,'value',ind)
	  return
	end
	Dat.ColMap = custom_cmap;
	if Dat.ColMapInverted
	  Dat.ColMap=flipud(Dat.ColMap);
	end
	set(H.FIG,'colormap',Dat.ColMap)
	setpref('Aedes','ColorMap',str)
  else
	h=errordlg({['Custom colormap not defined!',...
	  ' You can store a custom 256x3 colormap matrix CMAP',...
	  ' into Aedes preferences using the following command:'],...
	  '','setpref(''Aedes'',''CustomColorMap'',CMAP)'},...
	  'Custom colormap not defined!','modal');
	
	% Set popupmenu value
	if ispref('Aedes','ColorMap')
	  str=getpref('Aedes','ColorMap');
	else
	  str = 'gray';
	end
	tmp=get(H.COLMAP_POPUP,'string');
	ind=find(strcmpi(tmp,str));
	set(H.COLMAP_POPUP,'value',ind)
	return
  end
else
  eval(['Dat.ColMap=' str '(256);'])
  if Dat.ColMapInverted
	Dat.ColMap=flipud(Dat.ColMap);
  end
  set(H.FIG,'colormap',Dat.ColMap)
  
  setpref('Aedes','ColorMap',str)
end

catch
  aedes_errordump(lasterror);
end
end % function l_ChangeColormap(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Export Images
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_ExportImages(h,evd)
try

  % 0 = clim unlocked
  % 1 = clim locked global
  % 2 = clim locked slice
  % 3 = clim unlocked and auto-balanced
if Dat.isDataMixed
  if Dat.LockClim==2
    clim = Dat.SliceClim;
  elseif Dat.LockClim==0
    clim = Dat.OrigClim;
  elseif Dat.LockClim==3
    clim = Dat.OrigClim;
  else
    clim = Dat.Clim;
  end
else
  clim = Dat.Clim;
end
if ~Dat.isImageOverlayLoaded || Dat.isImageOverlayHidden
  overlay = [];
else
  overlay.ImageOverlay = Dat.ImageOverlay;
  overlay.isOverlayRGB = Dat.isOverlayRGB;
  overlay.ImOverlayMax = Dat.ImOverlayMax;
  overlay.ImOverlayMin = Dat.ImOverlayMin;
  overlay.ImageOverlayThold = Dat.ImageOverlayThold;
  overlay.ImageOverlayClim = Dat.ImageOverlayClim;
  overlay.ImageOverlayAlpha = Dat.ImageOverlayAlpha;
  overlay.ImageOverlayTholdDirPos = Dat.ImageOverlayTholdDirPos;
  overlay.ImageOverlayCmapStr = Dat.ImageOverlayCmapStr;
  overlay.ImageOverlayCmap = Dat.ImageOverlayCmap;
end

aedes_export_gui(DATA,ROI,clim,Dat.ColMap,Dat.RoiTransp,Dat.ShowRoiEdges,overlay,Dat.CurrentVol)

catch
  aedes_errordump(lasterror);
end
end % function l_ExportImages(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Unfold data
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_UnfoldData(h,evd)
try

%% Ask for confirmation
if Dat.isDataMixed
	resp=questdlg({'Choose unfolding direction:','',...
		['Left/Right (current Y-position = ',num2str(Dat.Slices(2)),')'],...
		'or',...
		['Up/Down (current X-position = ',num2str(Dat.Slices(1)),')']},...
		'Choose direction',...
		'Left/Right','Up/Down','Cancel','Left/Right');
else
	resp=questdlg({'Choose unfolding direction:','',...
		['Left/Right (current Y-position = ',num2str(Dat.Slices(2)),')'],...
		'',...
		['Front/Back (current Z-position = ',num2str(Dat.Slices(3)),')'],...
		'',...
		['Up/Down (current X-position = ',num2str(Dat.Slices(1)),')'],...
		'','To cancel: close the dialog.',''},...
		'Choose direction',...
		'Left/Right','Up/Down','Front/Back','Left/Right');
end
if isempty(resp) || strcmpi(resp,'Cancel')
  return
elseif strcmpi(resp,'Left/Right')
	if Dat.isDataMixed
		FoldInd = [0 -Dat.Slices(2)];
	else
		FoldInd = [0 -Dat.Slices(2) 0];
	end
elseif strcmpi(resp,'Up/Down')
	if Dat.isDataMixed
		FoldInd = [-Dat.Slices(1) 0];
	else
		FoldInd = [-Dat.Slices(1) 0 0];
	end
elseif strcmpi(resp,'Front/Back')
	FoldInd = [0 0 -Dat.Slices(3)];
end

%% Unfold data
if Dat.isDataMixed
  for ii=1:Dat.DataL
    DATA{ii}.FTDATA = circshift(DATA{ii}.FTDATA,FoldInd);
  end
else
  DATA{1}.FTDATA = circshift(DATA{1}.FTDATA,FoldInd);
end

% Unfold ROI(s)
if ~isempty(ROI)
	for kk = 1:length(ROI)
		if Dat.isDataMixed
			for ii=1:Dat.DataL
				ROI(kk).voxels{ii} = circshift(ROI(kk).voxels{ii},FoldInd);
			end
		else
			ROI(kk).voxels{1} = circshift(ROI(kk).voxels{1},FoldInd);
		end
	end
end


% Refresh data
l_DisplayData([],[])

% Refresh ROIs
l_RefreshRoi([],[])


catch
  aedes_errordump(lasterror);
end
end % function l_UnfoldData(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% KeyPressFcn - Handle keyboard shortcuts
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_KeyPressFcn(h,evd)
try
  
% If data is not loaded, return
if strcmpi(get(H.BLANK_FRAME,'visible'),'on')
  return
end

CurrentKey=lower(evd.Key);
CurrentModifier = evd.Modifier;

switch CurrentKey
  %% View X direction ----------------------------
 case 'f1'
  if Dat.ImageDim(Dat.DataInd,3)~=1 
    l_ChangeView([],[],1)
  end
  
  %% View Y direction ---------------------------
 case 'f2'
  if Dat.ImageDim(Dat.DataInd,3)~=1 
    l_ChangeView([],[],2)
  end
  
  %% View Z direction ---------------------------
 case 'f3'
  if Dat.ImageDim(Dat.DataInd,3)~=1 
    l_ChangeView([],[],3)
  end
  
  %% View 3D -------------------------
 case 'f4'
  if Dat.ImageDim(Dat.DataInd,3)~=1 
    l_ChangeView([],[],0)
  end
    
  %% Change + slice -----------------------------
 case {'rightarrow','pageup'}
  
  % Check that control is also pressed
  if strcmpi(CurrentKey,'rightarrow') && ~(length(CurrentModifier)==1 && ...
			strcmp(CurrentModifier{1},'control'))
    return
  end
  
  % Move Slice Slider by one notch
  sl_val = get(H.SL_SLIDER,'value');
  sl_min = get(H.SL_SLIDER,'min');
  sl_max = get(H.SL_SLIDER,'max'); 
  new_val = min(max(sl_val+1,sl_min),sl_max);
  set(H.SL_SLIDER,'value',new_val)
  l_ChangeSlice(H.SL_SLIDER,[],'slider')
  return
  
  %% Change - slice -----------------------------
 case {'leftarrow','pagedown'}
  
  % Check that control is also pressed
  if strcmpi(CurrentKey,'leftarrow') && ~( length(CurrentModifier)==1 && ...
			strcmp(CurrentModifier{1},'control') )
    return
  end
  
  % Move Slice Slider by one notch
  sl_val = get(H.SL_SLIDER,'value');
  sl_min = get(H.SL_SLIDER,'min');
  sl_max = get(H.SL_SLIDER,'max'); 
  new_val = min(max(sl_val-1,sl_min),sl_max);
  set(H.SL_SLIDER,'value',new_val)
  l_ChangeSlice(H.SL_SLIDER,[],'slider')
  return
  
  
  %% Zoom view in ----------------------------
  % Ctrl - downarrow
 case 'downarrow'
  
  % Check that control is also pressed
  if ~(length(CurrentModifier)==1 && strcmp(CurrentModifier{1},'control'))
    return
  end
  
  l_Zoom([],[],'+')
  
  %% Zoom view out -----------------------------
  % Ctrl - uparrow
 case 'uparrow'
  
  % Check that control is also pressed
  if ~(length(CurrentModifier)==1 && strcmp(CurrentModifier{1},'control'))
    return
  end
  
  l_Zoom([],[],'-')
  
  
  %% Toggle zoom view normalized
  % Ctrl - Return
 case 'return'
  
  % Check that control is also pressed
  if ~(length(CurrentModifier)==1 && strcmp(CurrentModifier{1},'control'))
    return
  end
  
  % Call zoom function
  l_Zoom([],[],'normalize')
  
	case 'c'
		
		if strcmpi(CurrentModifier,'control')
			if ~isempty(ROI) || Dat.AxView~=0
				l_RoiCopyCurrent;
			end
		end
		
		
	case 'v'
		
		if strcmpi(CurrentModifier,'control')
			l_RoiPasteCopied;
		end
 otherwise
  return
end

% Refresh
%drawnow

catch
  aedes_errordump(lasterror);
end
end % function l_KeyPressFcn(h,


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Save Results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_SaveResults(h,evd,filename)
try
  
%% If there isn't anything to save...
if isempty(ROI)
  h=errordlg('Cannot save results. ROI(s) not defined.','Cannot save results','modal');
  return
end

% Default name for file
[fp,fn,fe]=fileparts([Dat.HDR{1}.fpath,Dat.HDR{1}.fname]);
fname = fn;

% If the file is a VNMR file, extract the filename from the path name
if strcmpi(fn,'fid')
  Ind=find(fp==filesep);
  if isempty(Ind)
    if ispc
      Ind=find(fp=='/');
    else
      Ind=find(fp=='\');
    end
  end
  fname = fp(Ind(end)+1:end);
  fname = fname(1:end-4);
end
  
% Default path for ROI file
try
  fpath = getpref('Aedes','PutResFileDir');
catch
  fpath = '';
end
defaultfile = [fpath,fname];

% Loop while filename is not exepted
ok=false;
while ~ok
  
  %% Ask for a file name
  if isempty(filename)
    figure(H.FIG)
    [fn,fp,fi]=uiputfile({'*.roi;*.ROI;*.res;*.RES','Save ROI(s) and Statistics (*.roi,*.res)';...
                        '*.res;*.RES','Save Statistics (*.res)';...
                        '*.roi;*.ROI','Save ROI(s) (*.roi)'},'Save Results',...
                         defaultfile);
    if isequal(fn,0) || isequal(fp,0)
      return
    end
    filename = [fp,fn];
    if fi==1
      filetype='all';
    elseif fi==2
      filetype='res';
    else
      filetype='roi';
    end
  end
  if Dat.isDataMixed
	rotateflip = [Dat.DataRotation;Dat.DataFlip];
  else
	rotateflip = Dat.RotateFlip3d;
  end
  [done,msg]=aedes_saveres(DATA,ROI,filename,'SaveType',filetype,...
                     'waitbar',true,'rotateflip',rotateflip);
  if ~done && iscell(msg)
    h=errordlg(msg,'Error while saving Results/ROIs','modal');
    return
  elseif ~done && strcmpi(msg,'Overwrite cancel')
    ok=false;
    filename=[];
  elseif ~done
    h=errordlg(msg,'Error while saving Results/ROIs','modal');
    return
  else
    ok=true;
  end
end

% Update preferences
setpref('Aedes','PutResFileDir',fp)

catch
  aedes_errordump(lasterror);
end
end % function l_SaveResults(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate world coordinates from voxel coordinates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function coord=l_WorldCoord(h,evd)
	
	coord = [];
	if ~isfield(DATA{Dat.DataInd},'DataFormat')
	  return
	end
	
	% Do the coordinate transformation for different file formats
	switch lower(DATA{Dat.DataInd}.DataFormat)
	  case {'nifti(1)','nifti(2)'}
		
	  case 'vnmr'
		
	  otherwise
		
	end
	
	  
  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FIND PLUGINS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function plugins = l_FindPlugins(plugin_path)

	plugins.name = {};
	plugins.groupName = {};
	plugins.isGroup = [];
	plugins.fname = {};
	plugins.groupFolder = {};
	
	dir_struct=dir(plugin_path);
	plugin_folders = {dir_struct([dir_struct(:).isdir]).name};
	plugin_files = {dir_struct(~[dir_struct(:).isdir]).name};
	
	% Remove folders that begin with a dot
	ind = regexpi(plugin_folders,'^\.');
	plugin_folders = {plugin_folders{cellfun(@isempty,ind)}};
  
  % Remove files that start with a dot (hidden files)
  ind = regexpi(plugin_files,'^\.');
	plugin_files = {plugin_files{cellfun(@isempty,ind)}};
	
	% Remove files that don't have .m extension
	ind = regexpi(plugin_files,'\.m$');
	plugin_files = {plugin_files{~cellfun(@isempty,ind)}};
  
	
	% Go through the subfolders
	for kk=1:length(plugin_folders)
	  tmp_dir = dir([plugin_path,plugin_folders{kk}]);
	  tmp_files = {tmp_dir(~[tmp_dir(:).isdir]).name};
	  ind = regexpi(tmp_files,'\.m$');
	  tmp_files = {tmp_files{~cellfun(@isempty,ind)}};
    ind = regexpi(tmp_files,'^\.');
    tmp_files = {tmp_files{cellfun(@isempty,ind)}};
	  
	  plugins.groupFolder{end+1} = plugin_folders{kk};
	  plugins.groupName{end+1} = strrep(plugin_folders{kk},'_',' ');
	  plugins.isGroup(end+1) = true;
	  tmp_files = regexprep(tmp_files,'\.m','','ignorecase');
	  plugins.fname{end+1} = tmp_files;
	  tmp_files = cellfun(@(x)[upper(x(1)),x(2:end)],...
		tmp_files,'UniformOutput',false);
	  plugins.name{end+1} = regexprep(tmp_files,'_',' ');
	end
	
	% Go through the M-files in the plugins directory
	for ii=1:length(plugin_files)
	  plugins.groupFolder{end+1} = '';
	  plugins.groupName{end+1} = '';
	  plugins.isGroup(end+1) = false;
	  tmp_file = regexprep(plugin_files{ii},'\.m','','ignorecase');
	  plugins.fname{end+1} = tmp_file;
	  tmp_file = [upper(tmp_file(1)),tmp_file(2:end)];
	  plugins.name{end+1} = regexprep(tmp_file,'_',' ');
	end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Execute plugins
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_ExecutePlugin(h,evd)
try
% Get handle to the executed plugin function
fhandle = get(gcbo,'userdata');

% Additional information to be passed on to the plugin
AddInfo.hFigure = H.FIG;
AddInfo.hAxes = [H.IMAX1,H.IMAX2,H.IMAX3];
if Dat.isDataMixed
  AddInfo.CurrentSlice = Dat.DataInd;
else
  AddInfo.CurrentSlice = Dat.Slices;
end
if ~isempty(ROI)
  AddInfo.CurrentROI = get(H.ROI_EDIT,'value');
else
  AddInfo.CurrentROI = [];
end
AddInfo.CurrentVol = Dat.CurrentVol;
AddInfo.Clim = Dat.Clim;
AddInfo.isDataMixed = Dat.isDataMixed;
AddInfo.AxView = Dat.AxView;

% Add Overlay info
if ~isempty(Dat.ImageOverlay)
	AddInfo.ImageOverlay = Dat.ImageOverlay;
	AddInfo.ImageOverlayThold = Dat.ImageOverlayThold;
	AddInfo.ImOverlayMin = Dat.ImOverlayMin;
	AddInfo.ImOverlayMax = Dat.ImOverlayMax;
	AddInfo.ImageOverlayClim = Dat.ImageOverlayClim;
	AddInfo.ImageOverlayTholdDirPos = Dat.ImageOverlayTholdDirPos;
end


% Execute plugin
try
  feval(fhandle,DATA,ROI,AddInfo);
catch
  tmp=lasterror;
  h=errordlg({'An error occurred while executing plugin file',...
              ['"',tmp.stack(1).file,'"'],'',...
             'The following error was returned:',['"',tmp.message,'"'],...
             '','Please, fix the plugin'},'Error executing plugin!');
end

catch
  aedes_errordump(lasterror);
end
end % function l_ExecutePlugin(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Launch Header Browser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function l_LaunchHeaderBrowser(fh,evd)
   try
		 if ( isfield(DATA{Dat.DataInd},'PROCPAR') && ~isempty(DATA{Dat.DataInd}.PROCPAR) ) || ...
				 ( isfield(DATA{Dat.DataInd},'HDR') && isfield(DATA{Dat.DataInd}.HDR,'FileHeader') && ...
				 ~isempty(DATA{Dat.DataInd}.HDR.FileHeader) )
			 aedes_headerbrowser(DATA{Dat.DataInd});
		 else
			 errordlg('Cannot open file header browser. File header not found!','File Header Not Found!','modal');
		 end
	 catch
		 aedes_errordump(lasterror);
	 end
	end

function l_PrintLicense(h,evd)

  fprintf(1,'****************************************************************************\n');
  fprintf(1,'* Aedes - A graphical tool for analyzing medical images\n');
  fprintf(1,'*\n')
  fprintf(1,'* Copyright (C) 2006 Juha-Pekka Niskanen <juhapekka.niskanen@gmail.com>\n');
  fprintf(1,'*\n')
  fprintf(1,'* Department of Physics, Department of Neurobiology\n');
  fprintf(1,'* University of Eastern Finland, Kuopio, FINLAND\n');
  fprintf(1,'*\n')
  fprintf(1,'* This program may be used under the terms of the GNU General Public\n');
  fprintf(1,'* License version 2.0 as published by the Free Software Foundation\n');
  fprintf(1,'* and appearing in the file LICENSE.TXT included in the packaging of\n');
  fprintf(1,'* this program.\n');
  fprintf(1,'*\n')
  fprintf(1,'* This program is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE\n');
  fprintf(1,'* WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.\n');
  fprintf(1,'*\n')
  fprintf(1,'*\n')
  if Dat.MatlabVersion>=7.05
    fprintf(1,'* NOTE: You can suppress this license notification by ');
    fprintf(1,'<a href="matlab:setpref(''Aedes'',''ShowLicenseAtStartUp'',false)">clicking here</a>\n');
  else
    fprintf(1,'* NOTE: You can suppress this license notification by typing the following\n')
    fprintf(1,'* in the Matlab command window:\n');
    fprintf(1,'* setpref(''Aedes'',''ShowLicenseAtStartUp'',false)\n');
  end
  fprintf(1,'****************************************************************************\n');
  
end

end % function Res = Aedes(DATA)


