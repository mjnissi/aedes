function data = fmri_global_norm(data)
% FMRI_GLOBAL_NORM - Globally normalize fMRI data
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

if nargin == 0
  error('Too few input arguments!');
elseif nargin > 1
  error('Too many input arguments!');
end

if isstruct(data) && isfield(data,'FTDATA')
  data = data.FTDATA;
elseif ~isnumeric(data)
  error('Invalid first input argument!')
end

% Return error if data is not 4D
if ndims(data)~=4
  error('fMRI data must be in a 4D matrix!')
end
  
% Make sure that data is in either single or double form...
if ~any(strcmpi(class(data),{'single','double'}))
  data = single(data);
end

% Do global normalization
for ii=1:size(data,4)
  vol = data(:,:,:,ii);
  data(:,:,:,ii) = vol-sum(vol(:))/numel(vol);
end
  