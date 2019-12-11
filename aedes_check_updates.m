function [isUpdateAvailable,HeadRev,WorkingCopyRev,error_msg]=aedes_check_updates()
% AEDES_CHECK_UPDATES - Check if Aedes updates are available
%
% Synopsis: 
%	[isUpdateAvailable,HeadRev,WorkingCopyRev] = aedes_check_updates()
% 
% Description:
%	Checks if there are updates available for Aedes. The
%	"isUpdateAvailable" output argument is true if the head revision number
%   in SVN repository is larger than working copy revision. If the
%   "HeadRev" output argument is -1, there is a problem with the network.
%   If "HeadRev" is empty, there is a problem running SVN commands.
% 
% Examples:
%	
% 
% See also:
%       AEDES, AEDES_UPDATE

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



isUpdateAvailable = false;
HeadRev = [];
error_msg = '';

% Get current revision and repository URL
[rev,repo,workingcopy] = aedes_revision;
if isempty(repo)
  WorkingCopyRev = rev;
  return
end
WorkingCopyRev = rev;

% Get head repository revision
[HeadRev,msg]=l_GetHeadRevision(repo);
if isempty(HeadRev)
  % Something wrong with running SVN commands. Perhaps SVN is not
  % installed? Network problems?
  error_msg = msg;
else
  if rev<HeadRev
    isUpdateAvailable = true;
  end
end
  
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Get Head Repository Revision
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [HeadRev,msg]=l_GetHeadRevision(repo_url)
  
  HeadRev = [];
  msg = '';
  
  % Try to determine the head revision in svn
  if isunix
    % Unix/Linux
    [s,w] = unix(['svn info ',repo_url]);
    if s~=0
      msg=w;
      return
    end
  elseif ispc
    % Windows
    [s,w] = dos(['svn info ',repo_url]);
    if s~=0
      msg=w;
      return
    end
  else
    % Mac OS X
    [s,w] = unix(['svn info ',repo_url]);
    if s~=0
      msg=w;
      return
    end
  end
  
  % Scan the lines from the output
  C=textscan(w,'%s','delimiter','\n');
  if ~isempty(C)
    C=C{:};
  else
    return
  end
  
  % Get revision line
  rev_ind = find(strncmpi(C,'Revision:',4));
  if ~isempty(rev_ind)
    rev_str = C{rev_ind};
  else
    return
  end
  HeadRev=str2num(rev_str(11:end));
  
  