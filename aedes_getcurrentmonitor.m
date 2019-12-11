function cMon=aedes_getcurrentmonitor(hFig)
% AEDES_GETCURRENTMONITOR - Returns the monitor, where the pointer
%                           currently resides
%   
%
% Synopsis: 
%        CurrentMonitor=aedes_getcurrentmonitor(hFig);
%
% Description:
%        If called without input arguments, returns the index for the 
%        monitor where the mouse pointer is currently located (1=primary, 
%        2=secondary). If a figure handle is given as an input argument the 
%        function returns the index of the monitor where the figure is
%        current located.
%
% Examples:
%        
%        % Get the screen size of the monitor where the pointer is located
%        scrsz = get(0,'MonitorPositions');
%        sz = scrsz(aedes_getcurrentmonitor,:)
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

% Default monitor
cMon = 1;

try
	mon_pos = get(0,'MonitorPositions');
catch
	% If the above code fails, Matlab version is probably old and we cannot
	% detect monitor positions at all...
	return
end

% If only one monitor is attached
nMonitors = size(mon_pos,1);
if nMonitors==1
	return
end

switch nargin
	case 0
		
		% Get current pointer location
		pl = get(0,'PointerLocation');
		
		ind = find(pl(1)>=mon_pos(:,1));
		tmp_pos = mon_pos(ind,1);
		[mn,ind2]=min(abs(tmp_pos-pl(1)));
		
		cMon = find(mon_pos(:,1)==tmp_pos(ind2));
		
		if isempty(cMon)
			cMon = 1;
		end
		
	case 1
		
		% Check if the input argument is a figure handle
		if ~ishandle(hFig) || ~strcmpi(get(hFig,'type'),'figure')
			error('Input argument has to be a valid figure handle.')
		end
		
		% Get figure position
		fig_pos = get(hFig,'position');
		
		% Calculate center of the figure in pixels
		fig_c = fig_pos(1)+fig_pos(3)/2;
		if fig_c < 1
			fig_c=1;
		end
		
		ind = find(fig_c>=mon_pos(:,1));
		tmp_pos = mon_pos(ind,1);
		[mn,ind2]=min(abs(tmp_pos-fig_c));
		
		cMon = find(mon_pos(:,1)==tmp_pos(ind2));
		
		if isempty(cMon)
			cMon = 1;
		end
		
		%cMon = min(max(1,length(find(~[mon_pos(:,1)<=fig_c],1))+1),nMonitors);
		
	otherwise
		error('Too many input arguments.')
end


















