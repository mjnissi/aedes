function [rev,repo,workingcopy] = aedes_revision()
% AEDES_REVISION - Returns the current revision of Aedes
%   
%
% Synopsis: 
%        [rev,repo,workingcopy] = aedes_revision;
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


rev = [];
repo = '';

% Get path to the svn working copy
fname= mfilename('fullpath');
[fp,fn,fe] = fileparts(fname);
workingcopy = [fp,filesep];

% Check if working copy is under version control.
if ~(exist([workingcopy,'.svn'])==7)
  rev = '$Revision: 219 $';
  rev = str2num(rev(12:end-2));
  return
end

% Check the current repository and revision
if isunix
  s = [];
  w = [];
  [s,w] = unix(['svn info "',fp,'"']);
  if isempty(w)
	return
  end
elseif ispc
  s = [];
  w = [];
  [s,w] = dos(['svn info "',fp,'"']);
  if isempty(w)
	return
  end
else
  
end

% Scan the lines from the output
C=textscan(w,'%s','delimiter','\n');
if ~isempty(C)
  C=C{:};
else
  return
end

% Get URL line
url_ind = find(strncmpi(C,'URL:',4));
if ~isempty(url_ind)
  url_str = C{url_ind}(6:end);
else
  return
end
repo = url_str;

% Get revision line
rev_ind = find(strncmpi(C,'Revision:',4));
if ~isempty(rev_ind)
  rev_str = C{rev_ind};
else
  return
end
rev=str2num(rev_str(11:end));

% The next comment line is changed in every commit by the svncommit
% bash-script every time it is called so that this file "aedes_revision.m" is
% always in the list of committed files. DO NOT EDIT THE NEXT LINE!!!
% - Svn Hook -
