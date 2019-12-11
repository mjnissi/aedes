function ind = aedes_check_file_exist(fnames,fpath)
% AEDES_CHECK_FILE_EXIST - Check if files exist in a specified directory
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

% This function is a part of Aedes 2.0 - A graphical tool for analyzing 
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


ind=[];

%% Check input arguments
if nargin<2
  error('Too few input arguments!');
elseif nargin>2
  error('Too many input erguments!')
end

if ~iscell(fnames)
  fnames={fnames};
end

%% Get files in the specified directory
if ~isdir(fpath)
  error('Directory "%s" does not exist!',fpath)
end
d=dir(fpath);
dir_files={d(~[d(:).isdir]).name};

if isempty(dir_files)
  ind=false(size(fnames));
  return
end

%% Check if any of the files exist
if ~isunix % Check file names, in windows ignore case
  ind=ismember(lower(fnames),lower(dir_files));
else
  ind=ismember(fnames,dir_files);
end

return
% - EOF -
