function zz = aedes_trendest(z,alpha);
% AEDES_TRENDEST - Estimate trends from signal using smoothness priors
%
% Synopsis:
%       function zz = aedes_trendest(z,alpha)
%
% Description:
%       Estimate trend from signal z using smoothness priors. t is the time
%       scale for z. alpha is the smoothing parameter.
%
% Examples:
%
% See also:
%       
%

% This function is a part of Aedes - A graphical tool for analyzing 
% medical images
%
% Copyright (C) 2006 Juha-Pekka Niskanen <Juha-Pekka.Niskanen@uku.fi> and
% Mika Tarvainen.
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



%% Check dimensions of t and z -------------
if min(size(z)) == 1,
  z = z(:);
end

M = size(z,1);
if M<3
  zz=mean(z);
end

if nargin < 2, alpha = 1e4; end

% Create a sparse second difference matrix
e = ones(M,1);
D2 = spdiags([e -2*e e], 0:2, M-2, M);

% Create sparse identity matrix
H = speye(M,M);
zhat = (H+alpha^2*D2'*D2)\z;   % Smooth. priors estimate
zz = zhat;
