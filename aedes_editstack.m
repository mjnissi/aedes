function [filelist,rem_ind,add_ind,sort_ind] = aedes_editstack(DATA)
% AEDES_EDITSTACK - GUI for editing image stack
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



% Public variables
H = [];
Out = [];
Dat = [];
cancel = true;


% Get file names and paths
for ii=1:length(DATA)
  Dat.fnames{ii} = DATA{ii}.HDR.fname;
  Dat.fpaths{ii} = DATA{ii}.HDR.fpath;
  Dat.fullfiles{ii} = [DATA{ii}.HDR.fpath,...
                      DATA{ii}.HDR.fname];
end
H = l_DrawGUI;

% Initial selection check
l_CheckSelection([],[])

% Wait for quit
waitfor(H.FH)
if cancel
  clear H Dat cancel DATA
  filelist = {};
  rem_ind = [];
  add_ind = [];
  sort_ind = [];
  return
end

filelist = Dat.liststr;

if ispc
  rem_ind = ~ismember(lower(Dat.fullfiles),lower(filelist));
  [tmp,sort_ind] = ismember(lower(filelist),lower(Dat.fullfiles));
  add_ind = ~ismember(lower(filelist),lower(Dat.fullfiles));
else
  rem_ind = ~ismember(Dat.fullfiles,filelist);
  [tmp,sort_ind] = ismember(filelist,Dat.fullfiles);
  add_ind = ~ismember(filelist,Dat.fullfiles);
end
clear H Dat cancel DATA

return

%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Draw Main GUI
%
%%%%%%%%%%%%%%%%%%%%%%%%%
function H=l_DrawGUI()
  

%% Load default font and colors
%FigColor=get(0,'DefaultUicontrolBackgroundcolor');
GD=aedes_gui_defaults;
%GD.col.mainfig = FigColor;
fig_h = 305;
fig_w = 550;
fig_location = aedes_dialoglocation([fig_w,fig_h]);
fig_pos = [fig_location(1) fig_location(2) fig_w fig_h];

%% The main figure
H.FH = figure('position',fig_pos,...
              'Units','Pixel', ...
              'Name','Edit Image Stack', ...
              'Numbertitle','off', ...
              'Tag','im_rotate_gui', ...
              'Color',GD.col.mainfig, ...
              'Toolbar','none', ...
              'Menubar','none', ...
              'DoubleBuffer','on', ...
              'DockControls','off',...
              'renderer','painters',...
              'KeyPressFcn','',...
              'resize','off',...
              'windowstyle','modal');

% Options uipanel
H.OPTUIPANEL = uipanel('parent',H.FH,...
                       'units','pixel',...
                       'position',[5 40 fig_w-10 fig_h-45],...
					   'backgroundcolor',GD.col.frame);


% Buttons uipanel
H.BTNPANEL = uipanel('parent',H.OPTUIPANEL,...
                     'units','pixel',...
                     'position',[10 10 150 105+85+15+30],...
                     'title','Options',...
					 'backgroundcolor',GD.col.frame);



% File listbox
tmp = get(H.BTNPANEL,'position');
H.FILELBOX = uicontrol('parent',H.OPTUIPANEL,...
                       'units','pixel',...
                       'style','listbox',...
                       'position',[tmp(1)+tmp(3)+10 tmp(2) fig_w-tmp(1)-tmp(3)-35 ...
                   tmp(4)-8],...
                       'backgroundcolor','w',...
                       'string',Dat.fullfiles,...
                       'Min',0,'Max',2,...
                       'value',1,...
                       'CallBack',@l_CheckSelection);
tmp = get(H.FILELBOX,'position');
files_tx = uicontrol('parent',H.OPTUIPANEL,...
                     'units','pixel',...
                     'style','text',...
                     'position',[tmp(1) tmp(2)+tmp(4) 150 12],...
                     'string','Files (slices)',...
                     'horizontalalign','left',...
                     'fontweig','bold',...
					 'backgroundcolor',GD.col.frame);

% Options buttons
tmp=get(H.BTNPANEL,'position');
H.ADDFILESBTN = uicontrol('parent',H.BTNPANEL,...
                          'position',[5 tmp(4)-45 tmp(3)-10 25],...
                          'string','Add Files',...
                          'Callback',@l_AddFiles,...
                          'Enable','on');
tmp=get(H.ADDFILESBTN,'position');
H.REMFILESBTN = uicontrol('parent',H.BTNPANEL,...
                          'position',[tmp(1) tmp(2)-tmp(4)-3 tmp(3) tmp(4)],...
                          'string','Remove Files',...
                          'Callback',@l_RemoveSelected,...
                          'Enable','on');
tmp=get(H.REMFILESBTN,'position');
H.MOVETOPBTN = uicontrol('parent',H.BTNPANEL,...
                         'position',[tmp(1) tmp(2)-tmp(4)-45 tmp(3) tmp(4)],...
                         'string','Move to Top',...
                         'Callback',{@l_MoveFcn,'top'},...
                         'Enable','off');
tmp=get(H.MOVETOPBTN,'position');
H.MOVEUPBTN = uicontrol('parent',H.BTNPANEL,...
                        'position',[tmp(1) tmp(2)-tmp(4)-3 tmp(3) tmp(4)],...
                        'string','Move Up',...
                        'Callback',{@l_MoveFcn,'up'},...
                        'Enable','off');
tmp=get(H.MOVEUPBTN,'position');
H.MOVEDOWNBTN = uicontrol('parent',H.BTNPANEL,...
                          'position',[tmp(1) tmp(2)-tmp(4)-3 tmp(3) tmp(4)],...
                          'string','Move Down',...
                          'Callback',{@l_MoveFcn,'down'},...
                          'Enable','off');
tmp=get(H.MOVEDOWNBTN,'position');
H.MOVEBOTTOMBTN = uicontrol('parent',H.BTNPANEL,...
                            'position',[tmp(1) tmp(2)-tmp(4)-3 tmp(3) tmp(4)],...
                            'string','Move Bottom',...
                            'Callback',{@l_MoveFcn,'bottom'},...
                            'Enable','off');

% Cancel button
tmp = get(H.OPTUIPANEL,'position');
H.CANCELBTN = uicontrol('parent',H.FH,...
                       'units','pixel',...
                       'position',[tmp(1)+tmp(3)-70 5 70 30],...
                       'string','Cancel',...
                       'callback','delete(gcbf)');

% OK Button
tmp = get(H.CANCELBTN,'position');
H.OKBTN = uicontrol('parent',H.FH,...
                    'units','pixel',...
                    'position',[tmp(1)-tmp(3)-5 tmp(2) tmp(3) tmp(4)],...
                    'string','OK',...
                    'callback',@l_OKCallBack);
tmp = get(H.OKBTN,'position');
H.RESETBTN = uicontrol('parent',H.FH,...
                       'units','pixel',...
                       'position',[5 tmp(2) tmp(3) tmp(4)],...
                       'string','Reset list',...
                       'callback',@l_ResetList);


end % function H=l_DrawGUI()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reset List
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_ResetList(h,evd)

% Ask for confirmatio
if ~isempty(h)
  resp = questdlg('This will reset file list to default. Are you sure?',...
                  'Reset File List?','Yes','No','No');
  if strcmpi(resp,'No')
    return
  end
end

% Update file list
set(H.FILELBOX,'string',Dat.fullfiles,...
               'value',[],...
               'listboxtop',1)

% Check if buttons need to be disabled
l_CheckSelection([],[])

end % function l_ResetList(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check listbox selection
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_CheckSelection(h,evd)

% Button handles
btn_h = [H.ADDFILESBTN,...
        H.REMFILESBTN,...
        H.MOVETOPBTN,...
        H.MOVEUPBTN,...
        H.MOVEDOWNBTN,...
        H.MOVEBOTTOMBTN];

% Get selected files
val = get(H.FILELBOX,'value');
str = get(H.FILELBOX,'string');
val_max = length(str);

if isempty(val)
  set(btn_h(2:end),'enable','off')
  return
end
if any(val==1) && any(val==val_max)
  set(btn_h(3:end),'enable','off')
elseif any(val==1)
  set(btn_h(3:4),'enable','off')
  set(setdiff(btn_h,btn_h(3:4)),'enable','on')
elseif any(val==val_max)
  set(btn_h(5:6),'enable','off')
  set(setdiff(btn_h,btn_h(5:6)),'enable','on')
else
  set(btn_h,'enable','on')
end

% Enable Remove button
set(btn_h(2),'enable','on')

end % function l_CheckSelectio(h,  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Move files in the list
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_MoveFcn(h,evd,opt)

% Get listbox string and value
liststr=get(H.FILELBOX,'string');
listval = get(H.FILELBOX,'value');
if isempty(listval)
  return
end

% Get selected listbox items
sel_items = {liststr{listval}};
unsel_items = {liststr{setdiff(1:length(liststr),listval)}};

switch lower(opt)
 case 'top'
  
  % Determine move length
  if min(listval)>1
    move_length=diff([1 min(listval)]);
  end
  
 case 'up'
  if min(listval)>1
    move_length=1;
  end
  
 case 'down'
  if max(listval)<length(liststr)
    move_length=-1;
  end
  
 case 'bottom'
  if max(listval)<length(liststr)
    move_length=diff([length(liststr) max(listval)]);
  end
  
 otherwise
  error('%s\n%s\n%s','Congratulations! You should not be able to produce this error!',...
         'In addition, I''m not all that sorry and there is no one to send an error message.',...
         'Please, stop bothering me and go away...')
  
end

% Compute new values
newsel_values = listval-move_length;
newunsel_values = setdiff(1:length(liststr),newsel_values);

% Rebuild list
newlist = cell(1,length(liststr));
newlist(newsel_values) = sel_items;
newlist(newunsel_values) = unsel_items;

% Refresh listbox
if any(strcmpi(opt,{'top','up'}))
  lboxtop = newsel_values(1);
else
  lboxtop = max(1,newsel_values(end)-10);
end

set(H.FILELBOX,'string',newlist,...
               'value',newsel_values,...
               'listboxtop',lboxtop)

% Check if buttons need to be disabled
l_CheckSelection([],[])

end % function l_MoveFcn(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Remove selected files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_RemoveSelected(h,evd)

% Get listbox string and value
liststr=get(H.FILELBOX,'string');
listval = get(H.FILELBOX,'value');
if isempty(listval)
  return
end

% Throw an error message if the user tries to remove all files
if length(listval)==length(liststr)
  h=errordlg({'The file list cannot be empty!',...
             'At least one file has to remain on the list.'},...
             'List cannot be empty!','modal');
  return
end

% Remove selected files
liststr(listval)=[];

% Refresh listbox
set(H.FILELBOX,'string',liststr,...
               'value',[],...
               'listboxtop',1)

% Check if buttons need to be disabled
l_CheckSelection([],[])


end % function l_RemoveSelected(h,


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Add Files to list
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_AddFiles(h,evd)
  
% Default filepath
try
  filepath=getpref('Aedes','GetDataFileDir');
catch
  filepath='';
end

% Ask for a file
[filefilt,dataformats] = aedes_getfilefilter;
[fn, fp, fi] = aedes_juigetfiles(filefilt,...
  'Select Data File(s)',filepath);
if isequal(fn,0) || isequal(fp,0)
  return
end

CurrentFiles = get(H.FILELBOX,'string');

% Construct full file cell
AddedFiles = cell(length(fn),1);
for ii=1:length(fn)
  AddedFiles{ii} = [fp{ii},fn{ii}];
end

% Check if some of the added files already exist in the list
tmp={CurrentFiles{ismember(CurrentFiles,AddedFiles)}};
rm_ind=find(ismember(AddedFiles,tmp));
rem_files = {AddedFiles{rm_ind}};
AddedFiles(rm_ind) = [];

% Add unique files to the list
if ~isempty(AddedFiles)
  NewFiles = {AddedFiles{:}, CurrentFiles{:}}';
  set(H.FILELBOX,'string',NewFiles,...
                 'value',1:length(AddedFiles),...
                 'listboxtop',1)
  
  % Check if buttons need to be disabled
  l_CheckSelection([],[])
  
  % Warn about dublicates that weren't added
  if ~isempty(rem_files)
    h=warndlg({['The following files were not added' ...
               ' because they were dublicates:'],'',...
              rem_files{:}},'Warning',...
              'modal');
  end
else
  h=warndlg({['All files were skipped' ...
              ' because they were dublicates!'],...
             'No files added.'},'Warning',...
            'modal');
end



end % function l_ResetFileList(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OK button is pressed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_OKCallBack(h,evd)

cancel = false;
liststr = get(H.FILELBOX,'string');
Dat.liststr = liststr';
delete(H.FH);

end % function l_OKCallBack(h,

end
