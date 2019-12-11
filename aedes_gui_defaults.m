function Dat = aedes_gui_defaults;
% AEDES_GUI_DEFAULTS - Default fonts and colors for GUIs
%   
%
% Synopsis: 
%       Dat = aedes_gui_defaults;
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


% Default fonts
set(0,'defaultuicontrolfontname','helvetica')

% %% Load color for frame and figure
% if isunix
%   DefaultColor = [239 235 222]/255;
% else
%   DefaultColor=get(0,'DefaultUicontrolBackgroundcolor');
% end

if strcmpi(computer,'pcwin'),	% - Windows ------
  Dat.title_fs = 10;	        % FontSize
  Dat.text_fs = 8;
  Dat.map_fs = 8;
  Dat.ax_fs = 7;
  Dat.marker_fs = 5;
  Dat.boldline_w = 2;
else				% - Others -------
  Dat.title_fs = 14;	        % FontSize
  Dat.text_fs = 10;
  Dat.map_fs = 10;
  Dat.ax_fs = 8;
  Dat.marker_fs = 6;
  Dat.boldline_w = 3;
end				

% Colors for Windows
if ispc
  DefaultColor=get(0,'DefaultUicontrolBackgroundcolor');
  Dat.col.txt        = [0 0 0];
  Dat.col.ax         = [1 1 1];
  Dat.col.mainfig    = DefaultColor;
  %Dat.col.mainfig    = [0.95 0.95 0.95];
  %col.frame      = [0.9 0.9 0.9];
  %Dat.col.frame      = [0.86 0.86 0.86];
  Dat.col.frame      = DefaultColor;
  %Dat.col.frame      = [191 204 230]./255;
  Dat.col.frame_brd  = [0 0 0];
  Dat.col.subframe   = [0.759 0.759 0.8];
  Dat.col.shline	 = [0.2 0.2 0.6];
  Dat.col.shlinebg	 = [0.1 0.1 0.1];
  Dat.col.button     = [0.75 0.75 0.9];
  Dat.col.button2    = [0.7 0.7 0.8];
  Dat.col.radio_btn  = [0.72 0.72 0.82];
  Dat.col.toggle_btn = [0.72 0.72 0.82];
  Dat.col.edit       = [1 1 1];
  Dat.col.checkbox   = Dat.col.subframe; %[0.72 0.72 0.82];
  Dat.col.listbox    = [0.72 0.72 0.82];
  Dat.col.popup      = [1 1 1];
  %Dat.col.slider     = [0.398 0.398 0.578];
  Dat.col.slider     = [166 179 201]./255;
  Dat.col.normch     = Dat.col.mainfig;
  Dat.col.selch      = [0.35 0.85 0.35];
  Dat.col.rejch      = [0.85 0.35 0.35];
  
else
  DefaultColor = [239 235 222]/255;
  Dat.col.txt        = [0 0 0];
  Dat.col.ax         = [1 1 1];
  Dat.col.mainfig    = DefaultColor;
  %Dat.col.mainfig    = [0.95 0.95 0.95];
  %col.frame      = [0.9 0.9 0.9];
  %Dat.col.frame      = [0.86 0.86 0.86];
  Dat.col.frame      = DefaultColor;
  %Dat.col.frame      = [191 204 230]./255;
  Dat.col.frame_brd  = [0 0 0];
  Dat.col.subframe   = [0.759 0.759 0.8];
  Dat.col.shline	 = [0.2 0.2 0.6];
  Dat.col.shlinebg	 = [0.1 0.1 0.1];
  Dat.col.button     = [0.75 0.75 0.9];
  Dat.col.button2    = [0.7 0.7 0.8];
  Dat.col.radio_btn  = [0.72 0.72 0.82];
  Dat.col.toggle_btn = [0.72 0.72 0.82];
  %Dat.col.edit       = [0.92 0.92 1];
  Dat.col.edit       = [1 1 1];
  Dat.col.checkbox   = Dat.col.subframe; %[0.72 0.72 0.82];
  %Dat.col.listbox    = [0.72 0.72 0.82];
  Dat.col.listbox    = [0.95 0.95 0.95];
  Dat.col.popup      = [0.72 0.72 0.82];
  %Dat.col.slider     = [0.398 0.398 0.578];
  Dat.col.slider     = [166 179 201]./255;
  Dat.col.normch     = Dat.col.mainfig;
  Dat.col.selch      = [0.35 0.85 0.35];
  Dat.col.rejch      = [0.85 0.35 0.35];
  
end

% - Number for tmp-figure -
Dat.tmpfig = 99;
fig_handles = get(0,'children');
if ~isempty(fig_handles)
	while ismember(Dat.tmpfig,double(fig_handles)),
		Dat.tmpfig = Dat.tmpfig+1;
	end
end

% Detect if HG2 graphics are used
tmp = figure('units','pixels','position',[1 1 1 1],'visible','off');
if isa(tmp,'matlab.ui.Figure')
	Dat.HG2graphics = true;
else
	Dat.HG2graphics = false;
end
delete(tmp)

return

%% - EOF - %% 
