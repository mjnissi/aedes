function bw_out = aedes_roifill(bw,seed)
% AEDES_ROIFILL - Do flood fill operation (4-conn) for binary image
%   
%
% Synopsis:
%        bw_out = aedes_roifill(bw_in,seed);
%
% Description:
%       Flood fill operation (4-connected) for binary 2D images.
%
% Examples:
%
% See also:
%       AEDES

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

% Flood fill algorithm to be used
method = 1;

% Look for a quick return
if bw(seed(1),seed(2))
	% Seed pixel is 1, no need to do anything...
	bw_out = bw;
	return
end

% matrix of linear indices
sz = size(bw);
bw = ~bw; % Image complement
lind = reshape(1:prod(sz),sz(1),sz(2));

% Matrix of row and column indices
[xind,yind] = ndgrid(1:sz(1),1:sz(2));

switch method
  case 1
    % This is the simplest and slowest flood fill method...
    Q = [];
    bw_out = false(size(bw));
    if ~bw(lind(seed(1),seed(2)))
      return
    end
    Q(end+1) = lind(seed(1),seed(2));
    while ~isempty(Q)
      n = Q(1);
      ni = [xind(n) yind(n)];
      if bw(n)
        bw_out(n) = true;
      end
      Q=Q(2:end);
      %Q(1)=[];
      
      % Go west
      if ni(2)~=1
        nn = lind(ni(1),ni(2)-1);
        if bw(nn) && ~bw_out(nn)
          bw_out(nn) = true;
          Q(end+1) = nn;
        end
      end
      
       % Go east
      if ni(2)~=sz(2)
        nn = lind(ni(1),ni(2)+1);
        if bw(nn) && ~bw_out(nn)
          bw_out(nn) = true;
          Q(end+1) = nn;
        end
      end
      
       % Go south
      if ni(1)~=sz(1)
        nn = lind(ni(1)+1,ni(2));
        if bw(nn) && ~bw_out(nn)
          bw_out(nn) = true;
          Q(end+1) = nn;
        end
      end
      
       % Go north
      if ni(1)~=1
        nn = lind(ni(1)-1,ni(2));
        if bw(nn) && ~bw_out(nn)
          bw_out(nn) = true;
          Q(end+1) = nn;
        end
      end
      
      
    end
    
  case 2
    % To be written...
    
end

bw_out = ~bw | bw_out;
