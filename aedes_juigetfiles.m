function [filename,filepath,filterindex] = aedes_juigetfiles(filefilter,ftitle,defaultdir)
% AEDES_JUIGETFILES - A Java based file browser for selecting multiple files in
% non-sorted order from different folders.
%   
%
% Synopsis: 
%       [FILENAME, PATHNAME, FILTERINDEX] = AEDES_JUIGETFILES(FILTERSPEC,TITLE, DEFAULTPATH)
%
% Description:
%       The behavior of AEDES_JUIGETFILES is very similar to UIGETFILE with a few
%       exceptions. The UIGETFILE properties 'location' and 'multiselect'
%       are not available. Output arguments FILENAME and PATHNAME will
%       always be cell arrays, if cancel is not pressed. The last input
%       argument is not the default file but the default
%       directory. FILTERSPEC must always be a cell array.
%
%       NOTE: Utilizes heavily the Matlab Java interface and undocumented
%       and unsupported functions. This is due to change from version to
%       version and, thus, this function will probably be broken in future
%       Matlab releases. At the moment, however, AEDES_JUIGETFILES should
%       work with Matlab version from R14SP2 to R2007b.
%
% Examples:
%       [filename, pathname, filterindex] = aedes_juigetfiles( ...
%       {'*.m;*.fig;*.mat;*.mdl', 'All MATLAB Files (*.m, *.fig, *.mat, *.mdl)';
%        '*.m',  'M-files (*.m)'; ...
%        '*.fig','Figures (*.fig)'; ...
%        '*.mat','MAT-files (*.mat)'; ...
%        '*.mdl','Models (*.mdl)'; ...
%        '*.*',  'All Files (*.*)'}, ...
%        'Pick a file(s)',pwd);
%
% See also:
%       UIGETFILE

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


%% Parse input arguments
if nargin==0
  filefilter = {'*.*','All Files (*.*)'};
  ftitle = 'Select File(s)';
  defaultdir = 0;
elseif nargin==1
  ftitle = 'Select File(s)';
  defaultdir = 0;
elseif nargin==2
  defaultdir = 0;
elseif nargin>3
  error('Too many input arguments')
end
%keyboard
%% Default colors and fonts for gui
%FigColor=get(0,'DefaultUicontrolBackgroundcolor');
GD=aedes_gui_defaults;
FigColor = GD.col.mainfig;

fig_w = 870;%780;
fig_h = 545;%535;
fig_location = aedes_dialoglocation([fig_w,fig_h]);
fig_pos = [fig_location(1) fig_location(2) fig_w fig_h];
fh=figure('position',fig_pos,...
          'Name',ftitle, ...
          'Numbertitle','off', ...
          'Tag','aedes_juigetfiles_main_fig', ...
          'Color',FigColor,...%[236 233 216]./255,...%GD.col.mainfig, ...
          'Toolbar','none', ...
          'Menubar','none', ...
          'DockControls','off',...
          'renderer','painters',...
          'closereq','uiresume(gcbf)',...'windowstyle','modal',...
          'handlevisibility','off',...
          'resize','off',...
          'keypressfcn',@l_KeyPressFcn);
H.fh=fh;
if ~GD.HG2graphics
	set(H.fh,'DoubleBuffer','on')
end
%pos = [5 60 200 420];
%pos=[5 60 320 420];
pos=[5 70 320 420];
%% Check if defaultdir exists
if ischar(defaultdir)
  if ~isdir(defaultdir)
    defaultdir=0;
  end
end

%% Try to find icon path
tmp = which('aedes_juigetfiles');
if ~isempty(tmp)
  [fp,fn,fe]=fileparts(tmp);
  iconpath = [fp,filesep,'icons',filesep];
else
  iconpath = [pwd,filesep,'icons',filesep];
end

%% Disable warnings about deprecated functions in Matlab R2008a->
version_number = aedes_getmatlabversion;
%if version_number>=7.06
%  warning('off','MATLAB:uitree:DeprecatedFunction');
%  warning('off','MATLAB:uitreenode:DeprecatedFunction');
%end

%% Create uitree for directory browsing --------------
rooticonpath = [iconpath,'MyComputer.gif'];
if version_number>=7.06
  rootnode = uitreenode('v0','My Computer','My Computer',rooticonpath,false);
  tree = uitree('v0',fh,'Root',rootnode,'position',pos);%,'ExpandFcn',@l_ExpFcn);
else
  rootnode = uitreenode('My Computer','My Computer',rooticonpath,false);
  tree = uitree(fh,'Root',rootnode,'position',pos);%,'ExpandFcn',@l_ExpFcn);
end
drawnow
drawnow

%% Create uitree for files ---------------------------
pos_filetree = [pos(1)+pos(3)+5 pos(2) 220 pos(4)];
if version_number>=7.06
  filetree = uitree('v0',fh,'position',pos_filetree);
else
  filetree = uitree(fh,'position',pos_filetree);
end
drawnow
drawnow

%'NodeSelectedCallback',{@l_FileTreeNodeSelected,filetree},...
if version_number>=7.06
  fileroot = uitreenode('v0','Loading','Loading...','',false);
else
  fileroot = uitreenode('Loading','Loading...','',false);
end
filetree.setRoot(fileroot);
filetree.reloadNode(fileroot);
drawnow
drawnow


%% Create selected files uitree
pos_selfiletree=[pos(1)+220+pos(3)+100 pos(2) 220 pos(4)];
if version_number>=7.06
  selfiletreeroot = uitreenode('v0','Selected Files','0 Files Selected','',false);
  selfiletree = uitree('v0',fh,'Root',selfiletreeroot,'position',pos_selfiletree);
else
  selfiletreeroot = uitreenode('Selected Files','0 Files Selected','',false);
  selfiletree = uitree(fh,'Root',selfiletreeroot,'position',pos_selfiletree);
end
drawnow
drawnow


%% UICONTROLS ----------------------------------
% Folders text
%tmp=get(tree,'position');
dir_tx = uicontrol('parent',fh,...
                   'position',[pos(1) pos(2)+pos(4)+5 150 13],...
                   'style','text',...
                   'string','Folders',...
                   'fontweight','bold',...
                   'horizontalalign','left',...
                   'backgroundcolor',FigColor);%[236 233 216]./255);

% Files text
%tmp=get(filetree,'position');
files_tx = uicontrol('parent',fh,...
                     'position',[pos(1)+pos(3)+5 pos(2)+pos(4)+5 150 13],...
                     'style','text',...
                     'string','Files',...
                     'fontweight','bold',...
                     'horizontalalign','left',...
                     'backgroundcolor',FigColor);%[236 233 216]./255);

% Selected Files text
%tmp=get(selfiletree,'position');
selfiles_tx = uicontrol('parent',fh,...
                        'position',[pos(1)+pos(3)+220+100+5 pos(2)+pos(4)+5 150 13],...
                        'style','text',...
                        'string','Selected Files',...
                        'fontweight','bold',...
                        'horizontalalign','left',...
                        'backgroundcolor',FigColor);%[236 233 216]./255);



% Add and Remove buttons
tmp=pos_filetree;
addallbnt = uicontrol('parent',fh,...
  'units','pixel',...
  'position',[tmp(1)+tmp(3)+5 (tmp(4))/1.5+tmp(2) 85 25],...
  'string','Add all >>',...
  'tooltip','Add all files from current folder to list',...
  'style','pushbutton',...
  'callback',{@l_AddRemFiles,'addall',filetree,selfiletree});
addbtn = uicontrol('parent',fh,...
  'units','pixel',...
  'position',[tmp(1)+tmp(3)+5 (tmp(4))/2+tmp(2) 85 25],...
  'string','Add >>',...
  'style','pushbutton',...
  'tooltip','Add selected files to list',...
  'callback',{@l_AddRemFiles,'add',filetree,selfiletree});
tmp=get(addbtn,'position');
rembtn = uicontrol('parent',fh,...
  'units','pixel',...
  'position',[tmp(1) tmp(2)-25-5 85 25],...
  'string','<< Remove',...
  'tooltip','Remove selected files from list',...
  'style','pushbutton',...
  'callback',{@l_AddRemFiles,'remove',filetree,selfiletree});
tmp=pos_filetree;
rembtn = uicontrol('parent',fh,...
  'units','pixel',...
  'position',[tmp(1)+tmp(3)+5 (tmp(4))/4+tmp(2) 85 25],...
  'string','<< Remove all',...
  'tooltip','Remove all files from list',...
  'style','pushbutton',...
  'callback',{@l_AddRemFiles,'removeall',filetree,selfiletree});


% Show path -checkbox
tmp=pos_selfiletree;
showpathcb = uicontrol('parent',fh,...
                       'units','pixel',...
                       'position',[tmp(1) tmp(2)-17 150 15],...
                       'style','checkbox',...
                       'value',0,...
                       'backgroundcolor',FigColor,...%[236 233 216]./255,...
                       'string','Show full path',...
                       'callback',{@l_ShowHidePathNames,selfiletree});

% Back/Forward buttons
tmp=pos;
back_btn = uicontrol('parent',fh,...
                     'units','pixel',...
                     'position',[tmp(1) tmp(2)+tmp(4)+20 80 25],...
                     'string','<< Back',...
                     'style','pushbutton',...
                     'callback',{@l_BackForward,'back',tree},...
                     'enable','off',...
                     'busyaction','cancel',...
                     'Interruptible','off');
forward_btn = uicontrol('parent',fh,...
                        'units','pixel',...
                        'position',[tmp(1)+85 tmp(2)+tmp(4)+20 80 25],...
                        'string','Forward >>',...
                        'style','pushbutton',...
                        'callback',{@l_BackForward,'forward',tree},...
                        'enable','off',...
                        'busyaction','cancel',...
                        'Interruptible','off');

% Path edit
tmp=get(forward_btn,'position');
pathedit = uicontrol('parent',fh,...
                     'units','pixel',...
                     'position',[tmp(1)+tmp(3)+5 tmp(2) 690 25],...
                     'string','',...
                     'style','edit',...
                     'horizontalalign','left',...
                     'backgroundcolor','w',...
                     'callback',{@l_OpenDirectory,tree,[]});


% Open and Cancel buttons
tmp=pos_selfiletree;
openbtn = uicontrol('parent',fh,...
                    'units','pixel',...
                    'position',[tmp(1)+tmp(3)-80 10 80 25],...
                    'string','Open',...
                    'style','pushbutton',...
                    'userdata',0,...
                    'callback','set(gcbo,''userdata'',1),uiresume(gcbf)',...
                    'enable','off');
tmp=get(openbtn,'pos');
cancelbtn = uicontrol('parent',fh,...
                      'units','pixel',...
                      'position',[tmp(1)-85 10 80 25],...
                      'string','Cancel',...
                      'style','pushbutton',...
                      'userdata',1,...
                      'callback','uiresume(gcbf)');

% File Filter popup
tmp=pos;
filefilter_popup = uicontrol('parent',fh,...
                             'units','pixel',...
                             'position',[tmp(1)+105 tmp(2)-35 440 25],...
                             'string',{' '},...
                             'value',1,...
                             'style','popup',...
                             'backgroundcolor','w',...
                             'callback',{@l_DirNodeSelected,tree,filetree});
ff_tx = uicontrol('parent',fh,...
                  'units','pixel',...
                  'position',[tmp(1) tmp(2)-35 100 20],...
                  'string','Show Files of Type:',...
                  'horizontalalign','left',...
                  'style','text',...
                  'backgroundcolor',FigColor);%[236 233 216]./255);
if isunix
  set(ff_tx,'fontsize',8)
end
ffstr = {filefilter{:,2}};

% Custom Filter Edit
tmp=get(ff_tx,'position');
custom_filter_tx = uicontrol('parent',fh,...
  'units','pixel',...
  'position',[tmp(1) tmp(2)-30 100 20],...
  'string','Custom Filter:',...
  'horizontalalign','left',...
  'style','text',...
  'backgroundcolor',FigColor);
if isunix
  set(custom_filter_tx,'fontsize',8)
end
tmp=get(filefilter_popup,'position');
customfilter_popup = uicontrol('parent',fh,...
  'units','pixel',...
  'position',[tmp(1) tmp(2)-30 tmp(3)/2 25],...
  'string',' ',...
  'tooltip','Filter using regular expressions',...
  'style','edit',...
  'backgroundcolor','w',...
  'callback',{@l_DirNodeSelected,tree,filetree});


% Construct find string for regexp
%regexp_str = strrep(strrep(strrep({filefilter{:,1}},';','|'),'*.','\.'),'\.*','.*');
regexp_str={filefilter{:,1}};
for ii=1:numel(regexp_str)
  str=regexp_str{ii};
  if any(str==';')
    [delim,tmp_str]=regexp(str,';','match','split');
    tmp_str=regexprep(regexprep(regexprep(regexprep(tmp_str,'\*\.','\\\.'),'^([^\\\.\*])','\^$1'),'([^\*])$','$1\$'),'^\\\.\*$','\.\*');
    str = sprintf('(%s)|',tmp_str{:});
    str = str(1:end-1);
  else
    str=regexprep(regexprep(regexprep(regexprep(str,'\*\.','\\\.'),'^([^\\\.\*])','\^$1'),'([^\*])$','$1\$'),'^\\\.\*$','\.\*');
  end
  regexp_str{ii} = str;
end
set(filefilter_popup,'string',ffstr,'userdata',regexp_str)

H.tree = tree;
H.filetree = filetree;
H.selfiletree = selfiletree;
H.openbtn = openbtn;
H.cancelbtn = cancelbtn;
H.filefilter_popup = filefilter_popup;
H.pathedit = pathedit;
H.back_btn = back_btn;
H.forward_btn = forward_btn;
H.showpathcb = showpathcb;
H.addbtn = addbtn;
H.rembtn = rembtn;
H.history = {};
H.historyind = 1;
H.updateHistory = true;
H.iconpath = iconpath;
H.customfilter_popup = customfilter_popup;
setappdata(H.fh,'H',H)

%% Set tree callbacks and options
set(tree,'NodeExpandedCallback',{@l_ExpFcnAlt,tree})
set(filetree,'MultipleSelectionEnabled',true,...
             'NodeCollapsedCallback',{@l_FileTreeNodeCollapsed,filetree},...
             'NodeExpandedCallback',{@l_OpenDirectory,tree,1});
set(tree,'NodeSelectedCallback',{@l_DirNodeSelected,tree,filetree});
set(selfiletree,'DndEnabled',true);
set(selfiletree,'MultipleSelectionEnabled',true);
set(selfiletree,'NodeDroppedCallback',{@l_SelFileTreeNodeDropped,selfiletree});

% Expand My Computer
tree.expand(rootnode);

% Try to open default directory
if ~ischar(defaultdir)
  if isunix
	defaultdir = getenv('HOME');
  else
	defaultdir = getenv('USERPROFILE');
  end
  %defaultdir=pwd;
end
l_OpenDirectory([],[],tree,defaultdir)
drawnow

% Wait for exit
uiwait(H.fh)
if get(H.openbtn,'userdata')==0
  % Action canceled
  filename=0;
  filepath=0;
  filterindex=0;
  delete(H.fh)
  return
end

% Get selected files
filename = {};
filepath = {};
filterindex = get(H.filefilter_popup,'value');
selfileroot = handle(selfiletree.getRoot);
childcount=selfileroot.ChildCount;
for ii=1:childcount
  node = selfileroot.getChildAt(ii-1);
  value = node.getValue;
  [fpath,fname,fext]=fileparts(value);
  filename{ii}=[fname,fext];
  filepath{ii}=[fpath,filesep];
end

fig_h = H.fh;
clear H tree filetree selfiletree

% Close window
delete(fig_h)
clear aedes_juigetfiles
return


% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ % Expand function
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ function nodes=l_ExpFcn(tree,value)
% $$$ %disp('Calling ExpFcn...')
% $$$ set(gcf,'pointer','watch')
% $$$ drawnow
% $$$ if strcmpi(value,'My Computer')
% $$$   
% $$$   %% Use Java to determine available drives
% $$$   drives = {};
% $$$   r=java.io.File.listRoots;
% $$$   for n=1:length(r)
% $$$     drive s{n} = char(r(n).toString);
% $$$   end
% $$$   
% $$$   % Add drives to the tree
% $$$   %nodes = [];
% $$$   for ii=1:length(drives)
% $$$     if any(strcmpi(drives{ii},{'A:\','B:\'}))
% $$$       iconpath = '.\icons\floppyicon.gif';
% $$$     else
% $$$       iconpath = '.\icons\driveicon.gif';
% $$$     end
% $$$     nodes(ii) = uitreenode(drives{ii},drives{ii},iconpath,0);
% $$$   end
% $$$ else
% $$$   count = 0;
% $$$   ch = dir(value);
% $$$   if isempty(ch)
% $$$     tmp=tree.SelectedNodes;
% $$$     node = handle(tmp(1));
% $$$     
% $$$     % Set loaded state for the node back to false
% $$$     tree.setLoaded(node,false)
% $$$     
% $$$     % Set leafing back to "off" position
% $$$     %node.Leaf(false)
% $$$     
% $$$     % Refresh node
% $$$     tree.reloadNode(node)
% $$$     
% $$$     nodes=[];
% $$$     h=warndlg(['Cannot access "' value '"'],'Access error');
% $$$     set(gcf,'pointer','arrow')
% $$$     return
% $$$   else
% $$$     directories = {ch([ch(:).isdir]).name};
% $$$   end
% $$$   %files = {ch(~[ch(:).isdir]).name};
% $$$   
% $$$   for ii=1:length(directories)
% $$$     if (any(strcmp(directories{ii}, {'.', '..', ''})) == 0)
% $$$       count = count + 1;
% $$$       iconpath = '.\icons\Folder.gif';
% $$$       
% $$$       leafing=l_CheckLeafing([value, directories{ii}, filesep]);
% $$$       nodes(count) = uitreenode([value, directories{ii}, filesep], ...
% $$$                                 directories{ii}, iconpath, leafing);
% $$$     end
% $$$   end
% $$$   
% $$$   if count==0
% $$$     nodes = [];
% $$$   end
% $$$ end
% $$$ set(gcf,'pointer','arrow')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alternate expand function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_ExpFcnAlt(h,evd,tree)

H=getappdata(findall(0,'tag','aedes_juigetfiles_main_fig'),'H');
set(H.fh,'pointer','watch')
drawnow

evdnode = evd.getCurrentNode;
value=evdnode.getValue;
version_number = aedes_getmatlabversion;

if ~tree.isLoaded(evdnode)


	if strcmpi(value,'My Computer')
    
    %% Use Java to determine available drives, as using DIR is veeryy
    %% slooow
    drives = {};
    r=java.io.File.listRoots;
    for n=1:length(r)
      drives{n} = char(r(n).toString);
    end
    
    % Add drives to the tree
	for ii=1:length(drives)
	  if any(strcmpi(drives{ii},{'A:\','B:\'}))
		iconpath = [H.iconpath,'floppyicon.gif'];
	  else
		iconpath = [H.iconpath,'driveicon.gif'];
    end
    if version_number >= 7.06
      nodes(ii) = uitreenode('v0',drives{ii},drives{ii},iconpath,0);
    else
      nodes(ii) = uitreenode(drives{ii},drives{ii},iconpath,0);
    end
	end
	else
		
		if ispc && length(value)>=4 && strcmpi(value(1:2),'\\') && ...
				length(find(value=='\'))==3
			
			% Handle Windows network shares
			[t,s]=dos(['net view ',value]);
			if t~=0
				directories = {'.','..'};
			else
				% Use regexp to find share names
				tmp=regexp(s,'([\w-\$]+)\s+Disk','tokens');
				directories=cat(1,tmp{:});
			end
			count = 0;
		else
			count = 0;
			ch = dir(value);
			if isempty(ch)
				directories = {};
			else
				directories = {ch([ch(:).isdir]).name};
			end
		end
		
    
		if isempty(directories)
      %tmp=tree.SelectedNodes;
      %node = handle(tmp(1));
      
      % Set loaded state for the node back to false
      tree.setLoaded(evdnode,false);
      
      % Refresh node
      tree.reloadNode(evdnode);
      
      h=warndlg(['Cannot access "' value '"'],'Access error','modal');
      set(H.fh,'pointer','arrow')
      return
    else
      %directories = {ch([ch(:).isdir]).name};
			if isunix
				% Hide hidden directories in Linux
				ind = regexp(directories,'^\.');
				directories = {directories{cellfun('isempty',ind)}};
			end
		end
    
    for ii=1:length(directories)
      if (any(strcmp(directories{ii}, {'.', '..', ''})) == 0)
        count = count + 1;
        iconpath = [H.iconpath,'Folder.gif'];
        
        %leafing=l_CheckLeafing([value, directories{ii}, filesep]);
        leafing=false;
        if version_number >= 7.06
          nodes(count) = uitreenode('v0',[value, directories{ii}, filesep], ...
            directories{ii}, iconpath, leafing);
        else
          nodes(count) = uitreenode([value, directories{ii}, filesep], ...
            directories{ii}, iconpath, leafing);
        end
      end
    end
    
    if count==0
      % Set loaded state for the node back to false
      tree.setLoaded(evdnode,false);
      
      % Refresh node
      tree.reloadNode(evdnode);
      
      set(H.fh,'pointer','arrow')
      return
    end
  end
  
  %% Add child nodes
  version_number = aedes_getmatlabversion;
  if length(nodes)==1 && version_number==7.01
    % Make nodes a JavaArray
    tmpnodes = nodes;
    nodes=javaArray('com.mathworks.hg.peer.UITreeNode',1);
    nodes(1) = java(tmpnodes);
  end
  
  tree.add(evdnode,nodes);
  tree.setLoaded(evdnode,true);
  tree.reloadNode(evdnode);
end
set(H.fh,'pointer','arrow')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check leafing for directories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function leafing=l_CheckLeafing(path)

s=dir(path);
if isempty(s)
  leafing = true;
  return
end

if length(s([s(:).isdir]))>2
  leafing = false;
else
  leafing = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute when directory node is selected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_DirNodeSelected(hCaller,evd,tree,filetree)

H=getappdata(findall(0,'tag','aedes_juigetfiles_main_fig'),'H');
set(H.fh,'pointer','watch')
drawnow
version_number = aedes_getmatlabversion;

% Clear file tree view
filetree.setRoot([]);
if version_number >=7.06
  fileroot = uitreenode('v0','Loading','Loading...','',false);
else
  fileroot = uitreenode('Loading','Loading...','',false);
end
filetree.setRoot(fileroot);
filetree.reloadNode(fileroot);


% Get Root
%fileroot = handle(filetree.getRoot);

% Remove all children
%filetree.removeAllChildren(fileroot);
%filetree.reloadNode(fileroot);


% Get selected node
currentnode=tree.getSelectedNodes;
if isempty(currentnode)
	set(H.fh,'pointer','arrow')
  return
else
	currentnode=handle(currentnode(1));
end

% Get path
path = currentnode.getValue;
if strcmp(path,'My Computer')
  %selfiletree.add(selfileroot,new_nodes);
  %selfiletree.setLoaded(selfileroot,true);
  %selfiletree.expand(selfileroot);
  %selfileroot.setName([num2str(length(new_nodes)) ' Files Selected'])
  fileroot.setValue('My Computer');
  fileroot.setName('My Computer');
  filetree.reloadNode(fileroot)
  set(H.fh,'pointer','arrow')
  return
end


% Try to get file information from the path
if strncmp(path,'\\',2) && length(find(path=='\'))==3
	[t,s]=dos(['net view ',path]);
	if t~=0
		%directories = {'.','..'};
		ch = [];
	else
		% Use regexp to find share names
		tmp=regexp(s,'([\w-\$]+)\s+Disk','tokens');
		directories=cat(1,tmp{:});
		ch = 1;
	end
	filenames = {};
else
	ch=dir(path);
	filenames={ch(~[ch(:).isdir]).name};
	directories={ch([ch(:).isdir]).name};
end

if isempty(ch)
  %rootnode = uitreenode(['Cannot access ' path],['Cannot access 'path],'',false);
  %filetree.setRoot(rootnode);
  %filetree.reloadNode(rootnode);
  fileroot.setName(['Cannot access ' path]);
  fileroot.setValue(['Cannot access ' path]);
  filetree.reloadNode(fileroot)
  set(H.pathedit,'string',path)
  h=warndlg(['Cannot access "' path '"'],'Access error','modal');
  set(H.fh,'pointer','arrow')
  return
end

if H.updateHistory
  % Add path to history
  if H.historyind<length(H.history)
    H.history(H.historyind+1:end)=[];
  end
  if isempty(H.history) || ~strcmp(H.history{end},path)
    H.history{end+1} = path;
  end
  H.historyind = length(H.history);
  set(H.forward_btn,'enable','off')
  if length(H.history)>1
    set(H.back_btn,'enable','on')
  end
  
  % Keep maximum of 50 folders in history
  if length(H.history)>50
    H.history(1)=[];
  end
  setappdata(H.fh,'H',H)
end

% Get Name
name = currentnode.getName;

% Set up file tree root
%rootnode = uitreenode(path,path,'',false);
%filetree.setRoot(rootnode);
fileroot.setName(path);
fileroot.setValue(path);


% Filter filenames

filter_val=get(H.filefilter_popup,'value');
if ~isempty(hCaller) && hCaller~=H.customfilter_popup
  regexp_str = get(H.filefilter_popup,'userdata');
  regexp_str = regexp_str{filter_val};
else
  regexp_str = get(H.customfilter_popup,'string');
end
if ~strcmp(regexp_str,'.*')
  [s_ind,e_ind]=regexpi(filenames,regexp_str);
  ind=~cellfun('isempty',s_ind);
  filenames = {filenames{ind}};
  
%   if ~isempty(s_ind)
%     ind = (double(char(s_ind{:}))-double(char(e_ind{:})))~=0;
%     filenames = {filenames{logical(ind)}};
%   else
%     filenames = {};
%   end
end


if isunix
  % Hide hidden directories in Linux
  ind = regexp(directories,'^\.');
  directories = {directories{cellfun('isempty',ind)}};
end

% Set regexp filter string to custom filter editbox
if ~isempty(hCaller) && hCaller~=H.customfilter_popup
  set(H.customfilter_popup,'string',regexp_str)
end

% Create child nodes
count=0;
for ii=1:length(directories)
  if ~any(strcmp(directories{ii},{'.','..'}))
    count=count+1;
    if version_number>=7.06
      nodes(count)=uitreenode('v0',[path,directories{ii},filesep],directories{ii},[H.iconpath,'Folder.gif'],false);
    else
      nodes(count)=uitreenode([path,directories{ii},filesep],directories{ii},[H.iconpath,'Folder.gif'],false);
    end
  end
end
for ii=1:length(filenames)
  count=count+1;
  if version_number>=7.06
    nodes(count)=uitreenode('v0',[path,filenames{ii}],filenames{ii},[H.iconpath,'new.gif'],true);
  else
    nodes(count)=uitreenode([path,filenames{ii}],filenames{ii},[H.iconpath,'new.gif'],true);
  end
end

set(H.pathedit,'string',path)
if count==0;
  filetree.reloadNode(fileroot);
  set(H.fh,'pointer','arrow')
  return
end

%% Add child nodes
version_number = aedes_getmatlabversion;
if length(nodes)==1 && ( version_number==7.01 | version_number==7.02)
  % Make nodes a JavaArray
  tmpnodes = nodes;
  nodes=javaArray('com.mathworks.hg.peer.UITreeNode',1);
  nodes(1) = java(tmpnodes);
end

filetree.add(fileroot,nodes);
filetree.setLoaded(fileroot,true);
filetree.expand(fileroot);
filetree.reloadNode(fileroot)
set(H.fh,'pointer','arrow')
%filetree.add(rootnode,nodes);
%filetree.setLoaded(rootnode,true);
%filetree.expand(rootnode);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Files selected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_FileTreeNodeSelected(h,evd,filetree)
%disp('l_FileTreeNodeSelected')

nodes=filetree.getSelectedNodes;
if isempty(nodes)
  return
end

% Don't allow selection of directories
count=0;
for ii=1:length(nodes)
  if ~nodes(ii).isRoot
    count=count+1;
    selnodes(count)=nodes(ii);
  end
end

if count==0
  filetree.setSelectedNodes([]);
else
  % Reselect only files
  filetree.setSelectedNodes(selnodes);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filetree collapsed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_FileTreeNodeCollapsed(h,evd,filetree)

% Don't let the filetree root node to collapse
filetree.expand(filetree.getRoot);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Node Dropped to selected files tree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_SelFileTreeNodeDropped(h,evd,selfiletree)

tmp=handle(evd);
SourceNode=tmp.getSourceNode;
TargetNode=tmp.getTargetNode;

target_value = TargetNode.getValue;
source_value = SourceNode.getValue;
if SourceNode.isRoot || TargetNode.isRoot
  return
end
if strcmp(target_value,source_value)
  return
end


% Reorder nodes
selfileroot = handle(selfiletree.getRoot);
childcount=selfileroot.ChildCount;

% Construct a new javaArray for nodes
nodes=javaArray('com.mathworks.hg.peer.UITreeNode',childcount);


count=0;
fromUp = false;
if childcount ~= 0
  for ii=0:childcount-1
    tmp=selfileroot.getChildAt(ii);
    if strcmp(tmp.getValue,target_value)
      if fromUp
        count=count+1;
        nodes(count) = TargetNode;
        count=count+1;
        nodes(count) = SourceNode;
      else
        count=count+1;
        nodes(count) = SourceNode;
        count=count+1;
        nodes(count) = TargetNode;
      end
    elseif ~strcmp(tmp.getValue,source_value)
      count=count+1;
      nodes(count) = tmp;
    else
      fromUp = true;
    end
  end
end

%selfiletree.setRoot([]);
%selfileroot = uitreenode('Selected Files','Selected Files','',false);
%selfiletree.setRoot(selfileroot);
selfiletree.removeAllChildren(selfileroot);

% Add new nodes to the tree
selfiletree.add(selfileroot,nodes);
selfiletree.expand(selfileroot);
selfiletree.reloadNode(selfileroot);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD/REMOVE Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_AddRemFiles(h,evd,opt,filetree,selfiletree)

% Get handles
H=getappdata(findall(0,'tag','aedes_juigetfiles_main_fig'),'H');

version_number = aedes_getmatlabversion;


switch opt
 case {'add','addall'}
  
  % Look for a quick return
  if isempty(filetree.getSelectedNodes) && strcmpi(opt,'add')
    return
  end
  
  % Get Rootnode
  selfileroot = handle(selfiletree.getRoot);
  
   
   if get(H.showpathcb,'value')
     ShowPaths = true;
   else
     ShowPaths = false;
   end
   
   if strcmpi(opt,'add')
	 % Get selected nodes
	 selnodes=filetree.getSelectedNodes;
	 filecount = length(selnodes);
   elseif strcmpi(opt,'addall')
	 selnodes=filetree.getRoot;
	 filecount=selnodes.getChildCount;
   end
   
   % Get existing nodes in the selfiletree
   childcount=selfileroot.ChildCount;
   count=0;
   if childcount ~= 0
     for ii=0:childcount-1
       count=count+1;
       selfilenodes(count) = selfileroot.getChildAt(ii);
     end
   end
   
   % Add existing nodes to new_nodes
   selfilenodevalues={};
   for ii=1:childcount
     if ShowPaths
       name_str = selfilenodes(ii).getValue;
     else
       [fpath,fname,fext] = fileparts(selfilenodes(ii).getValue);
       name_str = [fname,fext];
     end
     if version_number >= 7.06
       new_nodes(ii) = uitreenode('v0',selfilenodes(ii).getValue,...
         name_str,...
         [H.iconpath,'new.gif'],...
         true);
     else
       new_nodes(ii) = uitreenode(selfilenodes(ii).getValue,...
         name_str,...
         [H.iconpath,'new.gif'],...
         true);
     end
     selfilenodevalues{ii}=selfilenodes(ii).getValue;
   end
   
   % Add selected nodes to new_nodes
   count=childcount;
   count2=0;
   for ii=childcount+1:childcount+filecount
	 count2=count2+1;
	 if strcmpi(opt,'add')
	   % Add selected
	   nodeVal=selnodes(count2).getValue;
	   isNodeRoot = selnodes(count2).isRoot;
	 else
	   % Add all files
	   nodeVal=selnodes.getChildAt(count2-1).getValue;
	   isNodeRoot = selnodes.getChildAt(count2-1).isRoot;
	 end
	 if ~any(strcmp(selfilenodevalues,nodeVal)) && ...
		 ~strcmp(nodeVal(end),filesep) && ~isNodeRoot
	   if ShowPaths
		 name_str = nodeVal;
	   else
		 [fpath,fname,fext] = fileparts(nodeVal);
		 name_str = [fname,fext];
	   end
	   count=count+1;
     if version_number >= 7.06
       new_nodes(count) = uitreenode('v0',nodeVal,...
         name_str,...
         [H.iconpath,'new.gif'],...
         true);
     else
       new_nodes(count) = uitreenode(nodeVal,...
         name_str,...
         [H.iconpath,'new.gif'],...
         true);
     end
	 end
   end
	
   
   
   if count==0
     return
   end
   
   % Enable Open button
   set(H.openbtn,'enable','on')
   
   
   version_number = aedes_getmatlabversion;
   if length(new_nodes)==1 && version_number==7.01
     % Make nodes a JavaArray
     tmpnodes = new_nodes;
     new_nodes=javaArray('com.mathworks.hg.peer.UITreeNode',1);
     new_nodes(1) = java(tmpnodes);
   end
   
   % Remove old nodes from the tree
   %selfiletree.setRoot([]);
   %selfileroot = uitreenode('Selected Files','Selected Files','',false);
   %selfiletree.setRoot(selfileroot);
   selfileroot.removeAllChildren;
   selfiletree.reloadNode(selfileroot);
   %selfiletree.repaint;
   %selfiletree.removeAllChildren(selfileroot);
   
   % Add new nodes to the tree
   selfiletree.add(selfileroot,new_nodes);
   %selfiletree.setLoaded(selfileroot,true);
   selfiletree.expand(selfileroot);
   selfileroot.setName([num2str(length(new_nodes)) ' Files Selected'])
   selfiletree.reloadNode(selfileroot)

   
 case {'remove','removeall'}
  %keyboard
  selnodes = selfiletree.getSelectedNodes;
  
  % Look for a quick return
  if isempty(selnodes) && strcmpi(opt,'remove')
    return
  end
  
  % Get values to the selected nodes
  selnodevalues = cell(1,length(selnodes));
  for ii=1:length(selnodes)
    tmp=handle(selnodes(ii));
    selnodevalues{ii}=tmp.getValue;
  end
  
  % Get Root
  selfileroot=selfiletree.getRoot;
  
  if strcmpi(opt,'removeall')
	selfileroot.removeAllChildren;
	selfiletree.reloadNode(selfileroot);
	
	% Disable Open button
	set(H.openbtn,'enable','off')
	selfileroot.setName('0 Files Selected')
	selfiletree.reloadNode(selfileroot)
	selfiletree.expand(selfileroot);
	return
  end
  
  % Get number of files in the tree
  childcount = selfileroot.getChildCount;
  
  count = 0;
  for ii=0:childcount-1
    tmp=selfileroot.getChildAt(ii);
    if ~any(strcmp(tmp.getValue,selnodevalues))
      count=count+1;
      new_nodes(count) = tmp;
    end
  end
  
  %selfiletree.setRoot([]);
  %selfileroot = uitreenode('Selected Files','Selected Files','',false);
  %selfiletree.setRoot(selfileroot);
  selfileroot.removeAllChildren;
  selfiletree.reloadNode(selfileroot);
  %selfiletree.removeAllChildren(selfileroot);
   
  % Add new nodes to the tree
  if count~=0
    selfiletree.add(selfileroot,new_nodes);
    selfileroot.setName([num2str(length(new_nodes)) ' Files Selected'])
    selfiletree.reloadNode(selfileroot)
  else
    % Unenable Open button
    set(H.openbtn,'enable','off')
    
    selfileroot.setName('0 Files Selected')
    selfiletree.reloadNode(selfileroot)
  end
  %selfiletree.setLoaded(selfileroot,true);
  selfiletree.expand(selfileroot);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_OpenDirectory(h,evd,tree,directory)

% Get handles
H=getappdata(findall(0,'tag','aedes_juigetfiles_main_fig'),'H');
set(H.fh,'pointer','watch')
drawnow

if isempty(directory)
  directory = get(h,'string');
elseif isnumeric(directory)
  tmp=evd.getCurrentNode;
  if tmp.isRoot
    set(H.fh,'pointer','arrow')
    return
  end
  directory=tmp.getValue;
end

% Check that directory exists
if ~isdir(directory) && ( ~ispc | ~strncmp(directory,'\\',2) )
  h=errordlg({['Cannot open folder "' directory '".'],['Make sure the path is ' ...
                      'correct.']},'Cannot open folder','modal');
  set(H.fh,'pointer','arrow')
  return
elseif strncmp(directory,'\\',2) && ispc
	% Test if the Windows computer is accessible
	ind=find(directory=='\');
	if ind<3
		str = directory;
	else
		str = directory(1:ind(3)-1);
	end
	[t,s]=dos(['net view ',str]);
	if t~=0
		h=errordlg(s,'Cannot open folder','modal');
		set(H.fh,'pointer','arrow')
		return
	end
end
if ~strcmp(directory(end),filesep)
  directory = [directory,filesep];
end

% See if we are in the specified directory already
%keyboard
tmp=tree.getSelectedNodes;
if ~isempty(tmp)
  tmp=handle(tmp(1));
  if strcmp(tmp.getValue,directory)
    set(H.fh,'pointer','arrow')
    return
  end
end

if H.updateHistory
  % Add path to history
  if H.historyind<length(H.history)
    H.history(H.historyind+1:end)=[];
  end
  if isempty(H.history) || ~strcmp(H.history{end},directory)
    H.history{end+1} = directory;
  end
  H.historyind = length(H.history);
  set(H.forward_btn,'enable','off')
  if length(H.history)>1
    set(H.back_btn,'enable','on')
  end
  
  % Keep maximum of 50 folders in history
  if length(H.history)>50
    H.history(1)=[];
  end
  setappdata(H.fh,'H',H)
end

% Get root
rootnode = tree.getRoot;

% Get Matlab version
version_number = aedes_getmatlabversion;

% Split the directory using file separator
if ispc
	sDir = regexp(directory,'\\','split');
	if isempty(sDir{1})
		% Windows network drive
		sDir([1 2])=[];
		sDir{1} = ['\\',sDir{1},'\'];
		childcount=rootnode.getChildCount;
		nodeFound = false;
		for kk=1:childcount
			tmp=rootnode.getChildAt(kk-1);
			name=tmp.getName;
			if strcmp(name,sDir{1})
				nodeFound = true;
			end
		end
		if ~nodeFound
			if version_number >= 7.06
				nodes = uitreenode('v0',sDir{1},sDir{1},'',0);
			else
				nodes = uitreenode(sDir{1},sDir{1},'',0);
			end
			if length(nodes)==1 && version_number==7.01
				% Make nodes a JavaArray
				tmpnodes = nodes;
				nodes=javaArray('com.mathworks.hg.peer.UITreeNode',1);
				nodes(1) = java(tmpnodes);
			end
			tree.add(rootnode,nodes);
			%tree.setSelectedNode(nodes);
			tree.setLoaded(rootnode,true);
			tree.reloadNode(rootnode);
			tree.setSelectedNode(nodes);
		end
	else
		sDir{1} = [sDir{1},'\'];
	end
	if isempty(sDir{end})
		sDir(end)=[];
	end
else
	sDir = regexp(directory,'/','split');
	sDir{1}='/';
end

%ind=find(directory==filesep);
%if ind(1)==1 && ind(2)==2
%	% Windows network share
%	ind(1)=[];
%end



%keyboard
drawnow
drawnow
start_ind=0;
currentnode = rootnode;
for ii=1:length(sDir)
  %if ii==1
  %  currentdir=sDir{ii};
  %else
    currentdir=sDir{ii};
  %end
  childcount=currentnode.getChildCount;
  for kk=1:childcount
    tmp=currentnode.getChildAt(kk-1);
    name=tmp.getName;
    if strcmp(name,currentdir)
      %tree.expand(tmp);
      %tree.setLoaded(tmp,true)
      %tree.reloadNode(tmp);
      currentnode=tmp;
      if ii~=length(sDir)
        tree.expand(currentnode)
      end
      %tree.setLoaded(tmp,true)
      %tree.setSelectedNode(currentnode);
      %tree.reloadNode(currentnode);
      %pause(12)
      drawnow
      break;
    end
  end
  %start_ind = ind(ii);
end

%tree.expand(currentnode)
%tree.setLoaded(currentnode,true)
tree.setSelectedNode(currentnode);
drawnow, pause(0.01)
tree.reloadNode(currentnode);
set(H.fh,'pointer','arrow')
%tree.setSelectedNode(currentnode)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SHOW/HIDE path names from selected files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_ShowHidePathNames(h,evd,selfiletree)

if get(h,'value')==1
  ShowPaths = true;
else
  ShowPaths = false;
end

% Redraw all nodes in the tree
selfileroot = handle(selfiletree.getRoot);
childcount=selfileroot.ChildCount;
if childcount==0
  return
end

count=0;
if childcount ~= 0
  for ii=0:childcount-1
    count=count+1;
    childnode=selfileroot.getChildAt(ii);
    if ShowPaths
      childnode.setName(childnode.getValue);
    else
      [fpath,fname,fext] = fileparts(childnode.getValue);
      childnode.setName([fname,fext])
    end
  end
end

% Reload root node
selfiletree.reloadNode(selfileroot);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BACK/FORWARD functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_BackForward(h,evd,opt,tree)

% Get handles
H=getappdata(findall(0,'tag','aedes_juigetfiles_main_fig'),'H');
%set(H.back_btn,'enable','inactive')
%set(H.forward_btn,'enable','inactive')

% Put up a flag that prevents the history to be updated
H.updateHistory = false;
setappdata(H.fh,'H',H)

switch opt
 case 'back'
  %keyboard
  if length(H.history)>1
    H.historyind = H.historyind-1;
    drawnow % Java seems to be slow to update...
    l_OpenDirectory([],[],tree,H.history{H.historyind})
    drawnow
    set(H.forward_btn,'enable','on')
    if H.historyind==1
      set(H.back_btn,'enable','off')
    end
  end
  
  % Put history update flag down
  H.updateHistory = true;
  setappdata(H.fh,'H',H)
  %H.history
  %H.historyind
 case 'forward'
  if length(H.history)==H.historyind
    return
  end
  
  if H.historyind<length(H.history)
    H.historyind = H.historyind+1;
    drawnow
    l_OpenDirectory([],[],tree,H.history{H.historyind})
    drawnow
    set(H.back_btn,'enable','on')
    if H.historyind==length(H.history)
      set(H.forward_btn,'enable','off')
    end
  end
  
  % Put history update flag down
  H.updateHistory = true;
  setappdata(H.fh,'H',H)
  %H.history 
  %H.historyind
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resize Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_Resize(h,evd)

% To be written...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Key press callback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_KeyPressFcn(h,evd)

% Exit with Esc
if strcmp(evd.Key,'escape')
  uiresume(gcbf)
end
