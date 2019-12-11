function [version_number,isImProcToolbox] = aedes_getmatlabversion()
% AEDES_GETMATLABVERSION - Returns the numerical Matlab (main) version
%   
%
% Synopsis: 
%       [version_number,isImProcToolbox] = aedes_getmatlabversion;
%
% Description:
%       Returns the numerical Matlab (main) version. The number is converted
%       to a proper decimal format (e.g. Matlab 7.6 -> 7.06, Matlab 7.10 ->
%       7.1). The availability of Image Processing Toolbox is also
%       returned in isImProcToolbox.
%
% Examples:
%       [version_number,isImProcToolbox] = aedes_getmatlabversion;
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

version_str = version;
isImProcToolbox = false;
ind = find(version_str=='.');
major_ver = str2double(version_str(1:ind(1)-1));
minor_ver = str2double(version_str(ind(1)+1:ind(2)-1));
version_number = str2double(sprintf('%d.%02d',major_ver,minor_ver));
if nargout>1
  try
    s=ver;
    toolboxes = {s(:).Name};
    if any(strcmpi(toolboxes,'Image Processing Toolbox'))
      isImProcToolbox = true;
    end
  catch
    isImProcToolbox = false;
  end
end
