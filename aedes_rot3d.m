function B = aedes_rot3d(A,k,dim)
% AEDES_ROT3D - Rotate 3D (or 4D) matrix in 90 degree steps in 3 dimensions
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


if nargin<2
  error('AEDES_ROT3D: Too few input arguments!')
elseif nargin<3
  dim = 3; % Use dim 3 as a default
end

if ~( isnumeric(A) || islogical(A) )
  error('First input argument must be a numerical 3D (or 4D)-matrix!')
end

if ~isnumeric(k) || ~any(k==[0 1 2 3 4])
  error('The second input argument has to be a scalar 0,1,2,3, or 4!')
end

dim_inds = 1:ndims(A);

switch dim
  case 1
		dim_inds = [1 3 2 4];
		%dim_inds(1:3)=dim_inds(3:-1:1);
	if k==1 % Rotate 90 degrees along dim 1 (rows)
	  B = flipdim(permute(A,dim_inds),1);
	elseif k==2 % Rotate 180 degrees along dim 1 (rows)
	  B = flipdim(flipdim(A,3),1);
    elseif k==3 % Rotate 270 degrees along dim 1 (rows)
	  B = flipdim(permute(A,dim_inds),3);
	elseif k==0 || k==4
	  B=A;
	end
  case 2
		dim_inds = [3 2 1 4];
		%dim_inds(2:3)=dim_inds(3:-1:2);
	if k==1 % Rotate 90 degrees along dim 2 (cols)
      B = flipdim(permute(A,dim_inds),3);
    elseif k==2 % Rotate 180 degrees along dim 2 (cols)
      B = flipdim(flipdim(A,3),2);
    elseif k==3 % Rotate 270 degrees along dim 2 (cols)
      B = flipdim(permute(A,dim_inds),2);
    elseif k==0 || k==4
      B=A;
	end
  case 3
    dim_inds = [2 1 3 4];
		%dim_inds(1:2)=dim_inds(2:-1:1);
    if k==1 % Rotate 90 degrees along dim 3
      B = flipdim(permute(A,dim_inds),1);
    elseif k==2 % Rotate 180 degrees along dim 3
      B = flipdim(flipdim(A,2),1);
    elseif k==3 % Rotate 270 degrees along dim 3
      B = flipdim(permute(A,dim_inds),2);
    elseif k==0 || k==4
      B=A;
	end
  otherwise
    return
end
