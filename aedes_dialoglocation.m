function fig_pos=aedes_dialoglocation(figure_size)
% AEDES_DIALOGLOCATION - Returns a good location for a dialog window
%   
%
% Synopsis: 
%        fig_pos=aedes_dialoglocation([figure_width,figure_height]);
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
% Copyright (C) 2006 Juha-Pekka Niskanen <Juha-Pekka.Niskanen@uef.fi>
% 
% Department of Applied Physics, Department of Neurobiology
% University of Eastern Finland, FINLAND
%
% This program may be used under the terms of the GNU General Public
% License version 2.0 as published by the Free Software Foundation
% and appearing in the file LICENSE.TXT included in the packaging of
% this program.
%
% This program is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
% WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

if nargin < 1
	error('Too few input arguments!')
end

hFig = gcbf;

% Get the primary display location
try
	scrsz = get(0,'MonitorPositions');
	scrsz = scrsz(1,:);
catch
	scrsz = get(0,'ScreenSize');
end

if ~isempty(hFig)
	old_units = get(hFig,'Units');
	try
		set(hFig,'Units','pixels');
		contsz = get(hFig,'position');
		set(hFig,'Units',old_units);
	catch
		set(hFig,'Units',old_units);
		return
	end
else
	contsz = scrsz;
end
fig_pos(1) = contsz(1)  + 1/2*(contsz(3) - figure_size(1));
fig_pos(2) = contsz(2)  + 2/3*(contsz(4) - figure_size(2));
