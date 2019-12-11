function C = aedes_cellsprintf(form_str,s)
% AEDES_CELLSPRINTF - CELLSTR with SPRINTF formalism
%   
%
% Synopsis: 
%
% Description:
%
% Examples:
%
% See also:
%       CELLSTR, SPRINTF

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


%% - Check input arguments -
if nargin<1
  error('Too few input arguments')
elseif nargin<2
  error('Second input argument not defined')
end

if ndims(s)~=2, error('S must be 2-D.'); end


if isempty(s)
  C = {''};
end

C = cell(size(s));
for ii=1:size(s,1)
  for jj=1:size(s,2)
    C{ii,jj} = sprintf(form_str,s(ii,jj));
  end
end


