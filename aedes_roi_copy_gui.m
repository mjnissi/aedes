function [slice_ind,dir_ind,roilabel,copytype] = aedes_roi_copy_gui(ROI,roi_ind)
% AEDES_ROI_COPY_GUI - 
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


% Initialize public variables
H = [];
slice_ind=[];
roilabel = [];
copytype=[];
dir_ind = [];
cancel=0;

%% Determine if ROI is for normal or mixed type data
if length(ROI(1).voxels)>1
  isMixedType = true;
else
  isMixedType = false;
end

H = l_DrawGUI(ROI,roi_ind);
fig_h = H.FH;

%% Wait for exit
waitfor(fig_h)
if cancel
  slice_ind=[];
  roilabel=[];
  dir_ind=[];
end
clear H
return

%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw GUI
%%%%%%%%%%%%%%%%%%%%%%%%%
function H=l_DrawGUI(ROI,roi_ind)

%% Parse ROI structure
roi_labels = {ROI(:).label};


%% Load default font and colors
FigColor=get(0,'DefaultUicontrolBackgroundcolor');
GD=aedes_gui_defaults;
GD.col.mainfig = FigColor;%[236 233 216]./255;
fig_h = 415;%395;
fig_w = 250;
fig_location = aedes_dialoglocation([fig_w,fig_h]);
fig_pos = [fig_location(1) fig_location(2) fig_w fig_h];


%% The main figure
H.FH = figure('position',fig_pos,...
              'Units','Pixel', ...
              'Name','Copy ROI', ...
              'Numbertitle','off', ...
              'Tag','roi_copy_fig', ...
              'Color',FigColor,...%[236 233 216]./255, ...
              'Toolbar','none', ...
              'Menubar','none', ...
              'DoubleBuffer','on', ... 
              'DockControls','off',...
              'renderer','painters',...
              'KeyPressFcn','',...
              'resize','off',...
              'windowstyle','modal');
fh = H.FH;

%% Select ROI text
h1 = uicontrol('parent',fh,...
               'units','pixels',...
               'position',[10 fig_h-10-15 170 15],...
               'style','text',...
               'horizontalalign','left',...
               'backgroundcolor',FigColor,...[236 233 216]./255,...
               'string','Select ROI(s) to be copied');

%% ROI label popupmenu
tmp=get(h1,'extent');
H.ROILBOX = uicontrol('parent',fh,...
                      'units','pixels',...
                      'position',[10 fig_h-25-65 230 65],...
                      'style','listbox',...
                      'string',roi_labels,...
                      'value',roi_ind,...
                      'min',0,'max',2,...
                      'backgroundcolor','w');


tmp=get(H.ROILBOX,'position');
%% Source select popup
% $$$ H.SOURCEPOPUP = uicontrol('parent',fh,...
% $$$                           'units','pixels',...
% $$$                           'position',[tmp(1) tmp(2)-tmp(4)-5 tmp(3) tmp(4)],...
% $$$                           'style','popup',...
% $$$                           'string',,...
% $$$                           'value',1,...
% $$$                           'backgroundcolor','w');

%% Direction uipanel
% $$$ h1 = uipanel('parent',fh,...
% $$$              'units','pixels',...
% $$$              'position',[5 tmp(2)-90-15 240 90],...
% $$$              'title','Direction',...
% $$$              'backgroundcolor',GD.col.mainfig);
% $$$ tmp=get(h1,'position');

%% Button group for direction radiobuttons
H.DIRBTNGRP = uibuttongroup('parent',fh,...
                      'units','pixels',...
                      'position',[5 tmp(2)-105-15 240 105],...
                      'title','Copy Direction',...
                      'backgroundcolor',GD.col.mainfig);
H.XRADIO = uicontrol('parent',H.DIRBTNGRP,...
                     'units','pixel',...
                     'position',[15 70 150 15],...
                     'style','radio',...
                     'string','X direction (row)',...'value',1,...
                     'backgroundcolor',GD.col.mainfig);
tmp=get(H.XRADIO,'position');
H.YRADIO = uicontrol('parent',H.DIRBTNGRP,...
                     'units','pixel',...
                     'position',[15 tmp(2)-20 150 15],...
                     'style','radio',...
                     'string','Y direction (column)',...'value',0,...
                     'backgroundcolor',GD.col.mainfig);
tmp=get(H.YRADIO,'position');
H.ZRADIO = uicontrol('parent',H.DIRBTNGRP,...
                     'units','pixel',...
                     'position',[15 tmp(2)-20 150 15],...
                     'style','radio',...
                     'string','Z direction (slice)',...'value',0,...
                     'backgroundcolor',GD.col.mainfig);
tmp=get(H.ZRADIO,'position');
H.VRADIO = uicontrol('parent',H.DIRBTNGRP,...
                     'units','pixel',...
                     'position',[15 tmp(2)-20 150 15],...
                     'style','radio',...
                     'string','V direction (volume)',...'value',0,...
                     'backgroundcolor',GD.col.mainfig,...
					 'enable','off');
if size(ROI(1).voxels{1},4)>1
  set(H.VRADIO,'enable','on')
  set(H.DIRBTNGRP,'SelectedObject',H.VRADIO)
else
  set(H.DIRBTNGRP,'SelectedObject',H.ZRADIO)
end
tmp=get(H.DIRBTNGRP,'position');


%% Copy type btngroup
H.CPTYPEGRP = uibuttongroup('parent',fh,...
                            'units','pixels',...
                            'position',[5 tmp(2)-70-10 240 70],...
                            'title','Copy type',...
                            'backgroundcolor',GD.col.mainfig);
H.APPENDRADIO = uicontrol('parent',H.CPTYPEGRP,...
                          'units','pixel',...
                          'position',[15 35 150 15],...
                          'style','radio',...
                          'string','Append',...
                          'backgroundcolor',GD.col.mainfig);
tmp=get(H.APPENDRADIO,'position');
H.OVERWRITERADIO = uicontrol('parent',H.CPTYPEGRP,...
                             'units','pixel',...
                             'position',[15 tmp(2)-25 80 15],...
                             'style','radio',...
                             'string','Overwrite',...
                             'backgroundcolor',GD.col.mainfig);
set(H.CPTYPEGRP,'SelectedObject',H.APPENDRADIO)
%if ~isMixedType
%  set([H.APPENDRADIO,H.OVERWRITERADIO],'enable','off')
%end
%% Slice selection btngroup
tmp=get(H.CPTYPEGRP,'position');
H.SLBTNGRP = uibuttongroup('parent',fh,...
                     'units','pixels',...
                     'position',[5 tmp(2)-70-10 240 70],...
                     'title','Copy ROI to...',...
                     'backgroundcolor',GD.col.mainfig);

H.ALLRADIO = uicontrol('parent',H.SLBTNGRP,...
                       'units','pixel',...
                       'position',[15 35 150 15],...
                       'style','radio',...
                       'string','All slices/volumes',...
                       'backgroundcolor',GD.col.mainfig);
tmp=get(H.ALLRADIO,'position');

H.CUSTOMRADIO = uicontrol('parent',H.SLBTNGRP,...
                          'units','pixel',...
                          'position',[15 tmp(2)-25 80 15],...
                          'style','radio',...
                          'string','Custom',...
                          'backgroundcolor',GD.col.mainfig);
tmp=get(H.CUSTOMRADIO,'position');
H.CUSTOMEDIT = uicontrol('parent',H.SLBTNGRP,...
                         'units','pixel',...
                         'position',[85 tmp(2) 145 18],...
                         'style','edit',...
                         'backgroundcolor','w',...
                         'horizontalalign','left',...
                         'string','1:end',...
                         'enable','off');
set(H.SLBTNGRP,'SelectedObject',H.ALLRADIO,...
               'SelectionChangeFcn',@l_SelectCustom)

tmp=get(H.SLBTNGRP,'position');
H.CANCELBTN = uicontrol('parent',fh,...
                        'units','pixel',...
                        'position',[175 tmp(2)-10-25 70 25],...
                        'string','Cancel',...
                        'callback',@l_Cancel);
H.OKBTN = uicontrol('parent',fh,...
                    'units','pixel',...
                    'position',[100 tmp(2)-10-25 70 25],...
                    'string','OK',...
                    'callback',@l_CheckValues);

%% If the ROI is for mixed type data, disable Y and Z buttons
if isMixedType
  set([H.YRADIO,H.XRADIO],'enable','off','value',0)
  set(H.ZRADIO,'value',1)
end

end % function H=l_DrawGUI()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select ALL/CUSTOM slices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function l_SelectCustom(h,evd)
  
  if evd.NewValue==H.CUSTOMRADIO
    set(H.CUSTOMEDIT,'enable','on')
  else
    set(H.CUSTOMEDIT,'enable','off')
  end
  
  end % function l_SelectCustom(h,


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check values when OK is pressed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function l_CheckValues(h,evd)
  
  %% Get handle to the selected radio button
  h=get(H.DIRBTNGRP,'SelectedObject');
  
  %% Get ROI label
  roi_vals=get(H.ROILBOX,'value');
  if isempty(roi_vals)
    warndlg('No ROI(s) selected! Please, select at least one ROI.',...
            'Error','modal');
    return
  end
  str=get(H.ROILBOX,'string');
  roilabel = {str{roi_vals}};
  
  %% ROI index
  %roi_ind=find(strcmp({ROI(:).label},str));
  roi_ind=1;
  
  %% Copy Type
  hCpType=get(H.CPTYPEGRP,'SelectedObject');
  if hCpType==H.APPENDRADIO
    copytype='append';
  else
    copytype='overwrite';
  end
  
  %% Determine ROI size in the selected direction
  if isMixedType
    sz=length(ROI(roi_ind).voxels);
  else  
    if h==H.ZRADIO
      sz=size(ROI(roi_ind).voxels{1},3);
      dir_ind=3;
    elseif h==H.YRADIO
      sz=size(ROI(roi_ind).voxels{1},2);
      dir_ind=2;
	elseif h==H.XRADIO
      sz=size(ROI(roi_ind).voxels{1},1);
      dir_ind=1;
	elseif h==H.VRADIO
	  sz=size(ROI(roi_ind).voxels{1},4);
      dir_ind=4;
    end
  end
  
  %% Check if all slices or custom radio button is selected
  h2=get(H.SLBTNGRP,'SelectedObject');
  
  if h2==H.ALLRADIO
    slice_ind = 1:sz;
  else %% Parse the custom string
    
    %% Get custom string
    custom_str = get(H.CUSTOMEDIT,'string');
    
    %% Check that the string contains only valid characters
    tmp_ind=ismember(custom_str,'0123456789end:,*+-/');
    if ~all(tmp_ind)
      h=warndlg(['The custom text "' custom_str '" contains invalid characters'],...
                'Invalid expression','modal');
      uiwait(h)
      return
    end
    
    % replace 'end' statement with ROI size
    custom_str=strrep(custom_str,'end','sz');
    
    % Try to evaluate the index string
    try
      eval(['slice_ind=[' custom_str '];']);
    catch
      h=warndlg({['Could not evaluate expression "' custom_str '".'],...
                 ['Custom string has to be a valid Matlab vector expression.']},...
                'Could not evaluate expression','modal');
      uiwait(h)
      return
    end
    
    % Sort indices and make sure that all indices are unique
    slice_ind=unique(slice_ind);
    
    % Make sure that min>=1 and max<=sz
    if max(slice_ind)>sz
      h=warndlg({['Could not evaluate expression "' custom_str '".'],...
                 'Maximum slice index exceeds the number of slices.'},...
                'Could not evaluate expression','modal');
      uiwait(h)
      return
    end
    if min(slice_ind)<1
      h=warndlg({['Could not evaluate expression "' custom_str '".'],...
                 ['Slice indices have to be larger than zero.']},...
                'Could not evaluate expression','modal');
      uiwait(h)
      return
    end
    
  end
  
  % If this point is reached, all is well and we can exit normally
  delete(H.FH)
  
  end % function l_CheckValues(h,


%%%%%%%%%%%%%%%%%%%%%%
% Cancel is pressed
%%%%%%%%%%%%%%%%%%%%%%
  function l_Cancel(h,evd)
  
  %% Set the cancel flag and resume with exit
  cancel=1;
  
  delete(H.FH)
  
  end % function l_Cancel(h,
end
