function ind = shiftind(data,shift)
% SHIFTIND - 
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

if length(data)==1
  ind = 1:data;
  dataLength = length(ind);
else
  ind = 1:length(data);
  dataLength = length(data);
end

if shift==0 % No shifting required
  return
elseif shift>0 % Shift right
  ind = [(shift+1):-1:2 1:dataLength-shift];
elseif shift<0 % Shift left
  shift = abs(shift);
  ind = [(shift+1):dataLength dataLength-1:-1:dataLength-shift];
end
