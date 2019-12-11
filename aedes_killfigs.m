function aedes_killfigs(opt)
% AEDES_KILLFIGS - Close all Matlab figures (also with hidden handles) with brute force
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

if nargin<1
  opt = 'all';
end

switch opt
  case 'all'
    % Get handles to all figures (with hidden and normal handles)
    H = findall(0,'type','figure');
    fprintf(1,'Killing %d figure(s)...\n',length(H))
    delete(H);
    
  case 'aedes'
    % Find handles to aedes windows
    tags = {'aedes_main_fig','aedes_juigetfiles_main_fig',...
      'header_browser_fig','aedes_overlay_controls_fig',...
      'aedes_resview_fig','aedes_export_fig'};
    for ii=1:length(tags)
      h=findall(0,'tag',tags{ii});
      delete(h)
    end
    
  otherwise
   error('Valid input arguments are ''all'' or ''aedes''!') 
end

