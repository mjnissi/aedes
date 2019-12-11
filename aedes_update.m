function done=aedes_update(opt)
% AEDES_UPDATE - Aedes update tool
%
% Synopsis: 
%	done=aedes_update(opt)
% 
% Description:
%	
% 
% Examples:
%	
% 
% See also:
%       AEDES_CHECK_UPDATES, AEDES_REVISION, AEDES

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


done = false;
ShowChangeLog = false;

if nargin==0
  opt = 'semiprompt';
end

% Check if updates are available
[isUpdateAvailable,HeadRev,WorkingCopyRev,error_msg]=aedes_check_updates();

if isUpdateAvailable
  if strcmpi(opt,'noprompt')
    fprintf(1,'Updating from revision %d to %d\n',...
      WorkingCopyRev,HeadRev)
  elseif any(strcmpi(opt,{'prompt','semiprompt'}))
    % Prompt for updates
    dialog_str = {'Updates for Aedes are available.',...
      '',...
      ['Your current revision: ',num2str(WorkingCopyRev)],...
      ['The latest revision: ',num2str(HeadRev)],...
      '',...
      'Do you want to install updates now?'};
    resp=questdlg(dialog_str,'Updates available. Install updates now?',...
      'Yes','No','Yes');
    if strcmpi(resp,'No')
      done = [];
      return
    end
    [cwh,txh]=aedes_calc_wait('Installing updates...');
  end
elseif isempty(HeadRev)
  % Something wrong with SVN
  if strcmpi(opt,'prompt')
    h=warndlg({'Could not check updates. Do you have hetwork access and is SVN properly installed?','',...
      'The following error was encountered:','',...
      error_msg},...
      'Could not check updates','modal');
    uiwait(h);
    return
  elseif any(strcmpi(opt,{'noprompt','semiprompt'}))
    fprintf(1,'Could not check updates.\n');
    fprintf(1,'The following error was encountered:\n');
    fprintf(1,error_msg);
    return
  end
else
  % No updates available, Aedes is up-to-date. Don't display anything if
  % semiprompt mode is selected
  if strcmpi(opt,'prompt')
    h=helpdlg('No updates available. Aedes is up-to-date.',...
      'Aedes up-to-date');
    uiwait(h);
  elseif strcmpi(opt,'noprompt')
    fprintf(1,'No updates available. Aedes is up-to-date.\n');
  end
  return
end


% Update to the latest revision
[fp,fn,fe]=fileparts(mfilename('fullpath'));
if isunix
  [s,w]=unix(['svn up "',fp,'"']);
  C=textscan(w,'%s','delimiter','\n');
  C=C{:};
elseif ispc
  [s,w]=dos(['svn up "',fp,'"']);
  C=textscan(w,'%s','delimiter','\n');
  C=C{:};
else
  % Mac OSX
  [s,w]=unix(['svn up "',fp,'"']);
  C=textscan(w,'%s','delimiter','\n');
  C=C{:};
end
if s~=0
  try
    delete(cwh)
  catch
  end
  if any(strcmpi(opt,{'prompt','semiprompt'}))
    h=errordlg({'Update failed because of following error:',...
      C{:}},...
      'Update failed','modal');
    uiwait(h);
    return
  elseif strcmpi(opt,'noprompt')
    fprintf(1,'\n********************* ERROR ***********************\n')
    fprintf(1,'Update failed because of following error:\n');
    fprintf(1,'%s\n',C{:});
    fprintf(1,'***************************************************\n')
    return
  end
end

try
  set(txh,'String','Installing updates...done')
  pause(0.5)
  delete(cwh)
catch
end


% Check for conflict or merge
c_ind = [];
g_ind = [];
for ii=1:length(C)
  if strncmpi(C{ii},'C  ',3)
    % Conflict
    c_ind(end+1)=ii;
  elseif strncmpi(C{ii},'G  ',3)
    % Merge
    g_ind(end+1)=ii;
  end
end

%%%%%%%%%%%%%%%%%%
%c_ind = [];
%g_ind = [];
%%%%%%%%%%%%%%%%%%

%% Everything went ok ------------------------
if isempty(c_ind) && isempty(g_ind)
  % Update successful
  if any(strcmpi(opt,{'prompt','semiprompt'}))
    resp=questdlg({['Updated successfully to revision ',num2str(HeadRev),'.'],...
      'Please close all Aedes windows and restart Aedes.','',...
      ['You can use "View Details..." the view the changelog.']},...
      'Update successful','OK','View Details...','View Details...');
    if strcmpi(resp,'View Details...')
      ShowChangeLog = true;
    else
      ShowChangeLog = false;
    end
  elseif strcmpi(opt,'noprompt')
    fprintf(1,'Updated successfully to revision %d\n',HeadRev);
    
    % Echo the update information to the command window
    fprintf(1,'%s\n',C{:});
    fprintf(1,'Update complete.\n');
  end
  
  % Clear functions just in case...
  clear functions
  
  done = true;
  
  % Show the Change log and information
  if ShowChangeLog
    if isunix
      [s,w] = unix(sprintf('svn log "%s" -r%d:%d',...
        fp,HeadRev,WorkingCopyRev+1));
      if s~=0
        
      end
    elseif ispc
      [s,w] = dos(sprintf('svn log "%s" -r%d:%d',...
        fp,HeadRev,WorkingCopyRev+1));
      if s~=0
        
      end
    else
      
    end
    W = textscan(w,'%s','delimiter','\n');
    W=W{:};
    web(['text://',...
      '<title>Update information and changelog</title>',...
      sprintf('<h1>Updated Aedes successfully to revision %d</h1>',...
      HeadRev),...
      sprintf('<h2>Updated files from revision %d to %d</h2>',...
      WorkingCopyRev,HeadRev),...
      sprintf('%s<br>',C{:}),...
      '<br>',...
      sprintf('<h2>Changelog from revision %d to %d</h2>',...
      WorkingCopyRev,HeadRev),...
      sprintf('%s<br>',W{:})]);
  end
  
  % Merges and/or conflicts detected ----------
elseif ~isempty(c_ind) || ~isempty(g_ind)
  
  % Echo the update information to the command window
  fprintf(1,'%s\n',C{:});
  fprintf(1,'Update complete.\n');
  
  % List of files with conflicts
  c_files = C(c_ind);
  
  % List of files with merges
  g_files = C(g_ind);
  
  % Update successful but with conflicts and/or merges
  if isempty(c_ind)
    % Only merges
    fprintf(1,'\n********************* WARNING *********************\n')
    fprintf(1,['Updated to revision %d but with errors.\n',...
      'The following files contain MERGES:\n'],HeadRev)
    fprintf(1,'%s\n',g_files{:})
    fprintf(1,'***************************************************\n')
    if any(strcmpi(opt,{'prompt','semiprompt'}))
      h = warndlg({['Updated to revision ',num2str(HeadRev),' but with errors.'],...
        '',...
        'The following files contain MERGES:',...
        '',...
        g_files{:},...
        ''},...
        'Updated with errors','modal');
      uiwait(h);
    end
    
  elseif isempty(g_ind)
    % Only conflicts
    fprintf(1,'\n********************* WARNING *********************\n')
    fprintf(1,['Updated to revision %d but with errors.\n',...
      'The following files contain CONFLICTS:\n'],HeadRev)
    fprintf(1,'%s\n',g_files{:})
    fprintf(1,'***************************************************\n')
    if any(strcmpi(opt,{'prompt','semiprompt'}))
      h = warndlg({['Updated to revision ',num2str(HeadRev),' but with errors.'],...
        '',...
        'The following files contain CONFLICTS:',...
        '',...
        c_files{:},...
        ''},...
        'Updated with errors','modal');
      uiwait(h);
    end
    
  else
    % Both conflicts and merges
    fprintf(1,'\n********************* WARNING *********************\n')
    fprintf(1,['Updated to revision %d but with errors.\n',...
      'The following files contain CONFLICTS:\n'],HeadRev)
    fprintf(1,'%s\n',c_files{:})
    fprintf(1,'\nThe following files contain MERGES:\n')
    fprintf(1,'%s\n',g_files{:})
    fprintf(1,'***************************************************\n')
    if any(strcmpi(opt,{'prompt','semiprompt'}))
      h = warndlg({['Updated to revision ',num2str(HeadRev),' but with errors.'],...
        '',...
        'The following files contain CONFLICTS:',...
        '',...
        c_files{:},...
        '',...
        'The following files contain MERGES:',...
        '',...
        g_files{:},...
        ''},...
        'Updated with errors','modal');
      uiwait(h);
    end
  end
  
  % Clear functions just in case...
  clear functions
  
  done = true;
  
end
