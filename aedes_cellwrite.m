function [done,msg] = aedes_cellwrite(incell,filename,varargin)
% AEDES_CELLWRITE - Write cell array to a text file
%   
%
% Synopsis: 
%        [done,msg] = aedes_cellwrite(incell,filename,varargin)
%
% Description:
%
% Examples:
%
% See also:
%        AEDES

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



% Defaults
done = false;
delim = ';';

% Check input arguments
if nargin<2
  error('Too few input arguments')
elseif ~iscell(incell)
  error('First input argument must be a cell array')
elseif ~ischar(filename)
  error('Filename must be of class char')
end

% Parse Input arguments
for ii=1:2:length(varargin)
  switch varargin{ii}
   case 'delimitter'
    if strcmpi(varargin{ii+1},'tab')
      delim = '\t';
    elseif strcmpi(varargin{ii+1},'space')
      delim = ' ';
    else
      delim=varargin{ii+1};
    end
   otherwise
    msg = sprintf('Unknown parameter "%s"!',varargin{ii});
    return
  end
end

% Open file for writing
fid = fopen(filename,'w');
if fid<0
  msg = sprintf('Could not open file "%s" for writing',filename);
  return
end

nRows = size(incell,1);
nCols = size(incell,2);

% Write cell to file
for ii=1:nRows
  for kk=1:nCols
    fprintf(fid,['%s',delim],incell{ii,kk});
  end
  fprintf(fid,'\r\n');
end

% Close file
fclose(fid);

% All went well...
done=true;
msg='';

% - EOF -
