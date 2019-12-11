function [h,txh] = aedes_calc_wait(str)
% AEDES_CALC_WAIT - Wait message box
%
% Synopsis: 
%	[h,txh] = aedes_calc_wait(str)
% 
% Description:
%	Displays wait message box with text 'str'. Output argument 'h' is
%	handle to message box and 'txh' is handle to text string. 
% 
% Examples:
%	[h,txh] = aedes_calc_wait('Demo');
%	pause(1)
%	set(txh,'String','Almost ready.');
%	pause(1)
%	delete(h)
% 
% See also:
%       AEDES_WBAR, AEDES

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
%
% Special Thanks to Perttu Ranta-aho for this marvelous piece of code...
% :-)

h = msgbox(str,'Processing...','help');

% - Replace 'OK'-button with text 'Wait...' -
txh = findall(findobj(h,'type','axes'),'tag','MessageBox');
set(txh,'Fontsize',8);
btn_h = findobj(h,'style','pushbutton');
set(btn_h,'units','normal')
pos = get(btn_h,'position');
set(btn_h,...
    'String','Please wait...',...
    'Style','Text',...
    'position',[0 pos(2) 1 pos(4)], ...
    'Fontsize',12)

% - Set pointer -
set(h,'Pointer','watch','Units','normal');

% - Commit changes -
set(h,'visible','on')
drawnow

