function copy_data_to_workspace(DATA,ROI,AddInfo)
% COPY_DATA_TO_WORKSPACE - Copy (assign) DATA and ROI structures to workspace
%   (this is an Aedes plugin)
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


%% Init figure size
FigColor=get(0,'DefaultUicontrolBackgroundcolor');
fig_h = 160;
fig_w = 245;
fig_location = aedes_dialoglocation([fig_w,fig_h]);
fig_pos = [fig_location(1) fig_location(2) fig_w fig_h];

fh = figure('units','pixel',...
            'position',...
            fig_pos,...
            'Name','Copy to workspace',...
            'numbertitle','off',...
            'Toolbar','none',...
            'Color',FigColor,...%[236 233 216]./255,...
            'Menubar','none',...
            'DoubleBuffer','on',... 
            'DockControls','off',...
            'renderer','painters',...
            'resize','off',...
            'Handlevisibility','off');

uipanel_h = uipanel('parent',fh,...
                    'units','pixel',...
                    'position',[10 40 fig_w-20 fig_h-50]);

data_tx = uicontrol('parent',uipanel_h,...
                    'units','pixel',...
                    'position',[10 80 205 18],...
                    'style','text',...
                    'backgroundcolor',FigColor,...%[236 233 216]./255,...
                    'horizontalalign','left',...
                    'string','Data structure variable name');
tmp=get(data_tx,'position');
data_edit = uicontrol('parent',uipanel_h,...
                      'units','pixel',...
                      'position',[tmp(1) tmp(2)-tmp(4) 205 20],...
                      'style','edit',...
                      'backgroundcolor','w',...
                      'horizontalalign','left',...
					  'tag','data_edit',...
                      'string','DATA',...
                      'userdata','DATA',...
					  'keypress',@l_VariableNameKeyPress);

tmp=get(data_edit,'position');
roi_tx = uicontrol('parent',uipanel_h,...
                   'units','pixel',...
                   'position',[tmp(1) tmp(2)-tmp(4)-10 205 18],...
                   'style','text',...
                   'backgroundcolor',FigColor,...%[236 233 216]./255,...
                   'horizontalalign','left',...
                   'string','ROI structure variable name');
tmp=get(roi_tx,'position');
roi_edit = uicontrol('parent',uipanel_h,...
                     'units','pixel',...
                     'position',[tmp(1) tmp(2)-tmp(4) tmp(3:4)],...
                     'style','edit',...
                     'backgroundcolor','w',...
                     'horizontalalign','left',...
					 'tag','roi_edit',...
                     'string','ROI',...
                     'userdata','ROI',...
                     'keypress',@l_VariableNameKeyPress);

%% OK and Cancel buttons
tmp=get(uipanel_h,'position');
ok_btn = uicontrol('parent',fh,...
                   'units','pixel',...
                   'position',[tmp(1)+79 tmp(2)-35 70 30],...
                   'style','pushbutton',...
                   'string','OK',...
                   'callback',@l_AssignVariables);
tmp=get(ok_btn,'position');
cancel_btn = uicontrol('parent',fh,...
                       'units','pixel',...
                       'position',[tmp(1)+tmp(3)+5 tmp(2:4)],...
                       'style','pushbutton',...
                       'string','Cancel',...
                       'callback','close(gcbf)');
dummy_btn = uicontrol('parent',fh,...
                      'units','pixel',...
                      'position',[-30 0 1 1],...
                      'style','pushbutton',...
					  'tag','dummy_btn',...
                      'string','',...
                      'callback','');
%set([data_edit,roi_edit],'userdata',dummy_btn)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable name checking
%
% NOTE: This doesn't work properly in newer Matlab versions...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_VariableNameKeyPress(h,evd)

% This function doesn't work in new Matlab versions so do nothing...
return
  
% Switch focus between uicontrols (a Matlab bug workaround)
uicontrol(dummy_btn)
uicontrol(h)

% Call drawnow twice to refresh the editbox (yet another Matlab bug
% workaround)
drawnow
drawnow

% Query the editbox string value
str=get(h,'string')

%set(h,'value',length(str))
%drawnow
%drawnow
if isempty(str)
  return
end

%% Check first letter
if ~any('abcdefghijklmnopqrstuvwxyz'==lower(str(1)))
  warndlg('The first character in the variable name has to be a letter.',...
          'Invalid character entered.','modal')
  set(h,'string',get(h,'userdata'))
  uicontrol(dummy_btn)
  uicontrol(h)
  drawnow
  drawnow
  return
end

%% Check other letters
if length(str)>1
  ind=ismember(lower(str(2:end)), ...
               '1234567890abcdefghijklmnopqrstuvwxyz_');
  if ~all(ind)
    warndlg({'Variable name can only contain the following characters:',...
             '','1234567890abcdefghijklmnopqrstuvwxyz_'},...
            'Invalid character entered.','modal')
    set(h,'string',get(h,'userdata'))
    uicontrol(dummy_btn)
    uicontrol(h)
    drawnow
    drawnow
    return
  end
end

set(h,'userdata',str)

end % function l_VariableNameKeyPress(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_AssignVariables(h,evd)
  
%% Check current workspace variable names
var_names = evalin('base','who');

%% Get DATA structure name
data_name = get(data_edit,'string');

%% Get ROI structure name
roi_name = get(roi_edit,'string');

%% Check that the variable names are not empty
if isempty(data_name)
  hh=warndlg('Data variable name is empty!',...
             'Empty variable name','modal');
  return
end
if isempty(roi_name)
  hh=warndlg('Roi variable name is empty!',...
             'Empty variable name','modal');
  return
end

%% Check variable names for invalid characters
str = {data_name,roi_name};
for ii=1:2
  %% Check first letter
  if ~any('abcdefghijklmnopqrstuvwxyz'==lower(str{ii}(1)))
	h_warn=warndlg('The first character in the variable name has to be a letter.',...
	  'Invalid character entered.','modal');
	uiwait(h_warn)
	uicontrol(dummy_btn)
	if ii==1
	  uicontrol(data_edit)
	else
	  uicontrol(roi_edit)
	end
	drawnow
	drawnow
	return
  end
  
  %% Check other letters
  if length(str{ii})>1
	ind=ismember(lower(str{ii}(2:end)), ...
	  '1234567890abcdefghijklmnopqrstuvwxyz_');
	if ~all(ind)
	  warndlg({'Variable name can only contain the following characters:',...
		'','1234567890abcdefghijklmnopqrstuvwxyz_'},...
		'Invalid character entered.','modal')
	  uicontrol(dummy_btn)
	  if ii==1
		uicontrol(data_edit)
	  else
		uicontrol(roi_edit)
	  end
	  drawnow
	  drawnow
	  return
	end
  end
end


%% Check if variables are being overwritten
if any(strcmpi(data_name,var_names))
  resp=questdlg({['A variable named "',data_name,...
                  '" already exists in the current workspace.'],'',...
                'Overwrite variable?'},...
                'Overwrite workspace variable?','Yes','No','No');
  if strcmpi(resp,'no')
    return
  end
end
if any(strcmpi(roi_name,var_names))
  resp=questdlg({['A variable named "',roi_name,...
                  '" already exists in the current workspace.'],'',...
                 'Overwrite variable?'},...
                'Overwrite workspace variable?','Yes','No','No');
  if strcmpi(resp,'no')
    return
  end
end

%% Assign variables into workspace
if length(DATA)==1
  assignin('base',data_name,DATA);
else
  assignin('base',data_name,DATA);
end
assignin('base',roi_name,ROI);

%% Close window
delete(fh)

end % function l_AssignVariables(h,

end
