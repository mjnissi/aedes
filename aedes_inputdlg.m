function Resp=aedes_inputdlg(prompt_str,title_str,default_answer)
% AEDES_INPUTDLG - Input dialog box
%   
%
% Synopsis: 
%        Resp=aedes_inputdlg(prompt_str,title_str,default_answer)
%
% Description:
%
%
% Examples:
%
% See also:
%        AEDES

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


if nargin==0
  prompt_str = 'Input:';
  title_str = '';
  default_answer = '';
elseif nargin==1
  title_str = '';
  default_answer = '';
elseif nargin==2
  default_answer = '';
elseif nargin>3
  error('Too many input arguments!')
end

if iscell(default_answer)
  default_answer = default_answer{1};
end

Resp = {};
callerFig = gcbf;

%% Dialog figure -------------------
GD = aedes_gui_defaults;
fig_w=215;
fig_h=85;
fig_location = aedes_dialoglocation([fig_w,fig_h]);
fig_pos = [fig_location(1) fig_location(2) fig_w fig_h];

fh=dialog('position',fig_pos,...
  'Visible','on', ...
  'KeyPressFcn',@l_KeyPressFcn, ...
  'Name',title_str, ...
  'Pointer','arrow', ...
  'Units','pixels', ...
  'HandleVisibility','callback', ...
  'Color',GD.col.frame, ...
  'WindowStyle','modal', ...
  'userdata',false,...
  'resize','off');
if ~GD.HG2graphics
	set(fh,'DoubleBuffer','on')
end


%% Title ---------------------------
title_h = uicontrol('parent',fh,...
  'units','pixels',...
  'position',[5 fig_h-20 fig_w-10 15],...
  'string',prompt_str,...
  'backgroundcolor',GD.col.frame,...
  'style','text',...
  'KeyPressFcn',@l_KeyPressFcn, ...
  'horizontalalign','left');


%% Edit box ------------------------
tmp = get(title_h,'position');
edit_h = uicontrol('parent',fh,...
  'units','pixels',...
  'position',[tmp(1) tmp(2)-25 tmp(3) 25],...
  'string',default_answer,...
  'style','edit',...
  'horizontalalign','left',...
  'KeyPressFcn',@l_KeyPressFcn, ...
  'backgroundcolor','w');

%% Cancel button -------------------
tmp = get(edit_h,'position');
cancelbtn_h = uicontrol('parent',fh,...
  'units','pixels',...
  'position',[tmp(1)+tmp(3)-60 tmp(2)-35 60 30],...
  'String','Cancel',...
  'KeyPressFcn',@l_KeyPressFcn, ...
  'Callback',@l_CancelButtonPressed);


%% OK-button -----------------------
tmp = get(cancelbtn_h,'position');
okbtn_h = uicontrol('parent',fh,...
  'units','pixels',...
  'position',[tmp(1)-60-5 tmp(2) tmp(3:4)],...
  'String','OK',...
  'KeyPressFcn',@l_KeyPressFcn, ...
  'Callback',@l_OKButtonPressed);


% Set focus to editbox
uicontrol(edit_h);

%% Wait for dialog figure
uiwait(fh);

if ishandle(fh)
    Resp={};
    if get(fh,'UserData')
	  Resp={get(edit_h,'string')};
    end
    delete(fh);
else
    Resp={};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_KeyPressFcn(h,evd)
switch(evd.Key)
 case {'return'}
   % Make sure that text visible in the editbox is registered
   uicontrol(findall(gcbf,'style','text'))
   uicontrol(findall(gcbf,'style','edit'))
   set(gcbf,'userdata',true);
   uiresume(gcbf);
 case {'escape'}
  delete(gcbf);
end


function l_OKButtonPressed(h,evd)

set(gcbf,'userdata',true)
uiresume(gcbf)

function l_CancelButtonPressed(h,evd)

delete(gcbf)


