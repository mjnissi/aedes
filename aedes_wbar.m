function H = aedes_wbar(varargin)
% AEDES_WBAR - Slightly modified Matlab Waitbar...
% Synopsis: 
%       See the WAITBAR function synopsis by typing 'help waitbar' into
%       the Matlab command prompt.
%
% Description:
%
% Examples:
%
% See also:
%        WAITBAR

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



if nargout==1
  h=waitbar(varargin{1},'');
  %h = waitbar(varargin{1},varargin{2:end});
  set(h,'Name',sprintf('Processing...%.0f%%',varargin{1}*100))
  set(get(findobj(h,'type','axes'),'Title'),'fontsize',8,...
	'interpreter','none',...
	'string',varargin{2})
  H=h;
elseif nargout==0
  waitbar(varargin{1:end})
  set(varargin{2},'Name',['Processing...' sprintf('%.0f',varargin{1}*100) '%'])
end

