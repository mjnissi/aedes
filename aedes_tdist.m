function p = aedes_tdist(x,v)
% TDIST - Calculates Student's cumulative distribution function.
% 
% Synopsis: 
%       p = aedes_tdist(X,V)
%
% Description:
%       Calculates Student's cumulative distribution function with V
%       degrees of freedom at the values in X.
%
% Examples:
%
% See also:
%        AEDES_FMRI_ANALYSIS

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

if nargin < 2,
    error('Too few input arguments.');
end

if numel(v)==1
  v = ones(size(x),class(v)).*v;
end

% Initialize P.
if isa(x,'single') || isa(v,'single')
    p = NaN(size(x),'single');
else
    p = NaN(size(x));
end

nans = (isnan(x) | ~(0<v));

% Handle special cases of degrees of freedom
ind = (v == 1);
if any(ind(:))
  p(ind) = 0.5+(1/pi).*atan(x);
end
ind = (v == 2);
if any(ind(:))
  p(ind) = 0.5*(1+x./sqrt(2+x.^2));
end

ind = isnan(p);
if any(ind(:))
  % When degrees of freedom are other than 1 or 2, the cumulative
  % distribution function is calculated as an incomplete beta function
  p(ind) = betainc((x(ind)+sqrt(x(ind).^2+v(ind)))./(2*sqrt(x(ind).^2+v(ind))),v(ind)/2,v(ind)/2);  
end
p(x == 0 & ~nans) = 0.5;