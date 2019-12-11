function aedes_resviewer(Res_in)
% AEDES_RESVIEWER - View and convert Aedes statistics results
%   
%
% Synopsis: 
%       aedes_resviewer(Res)
%
% Description:
%       The AEDES_RESVIEWER function opens a graphical user interface (GUI) for
%       handling Res-structures or *.res files saved from AEDES. If
%       the function is called with the optional input argument Res, the
%       GUI opens in a single-file mode. Multi-file mode is used if the
%       function is called without input arguments.
%
% Examples:
%       reviewer % Open the user interface
%
%       % or
%
%       tmp=load('resfile.res','mat'); % Read saved Res-file
%       Res=tmp.Res;                   
%       aedes_resviewer(Res)                 % Open the loaded Res-structure
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


% Public variables
H = [];   % Object handles
Dat = []; % GUI internal variables
Dat.SingleFileMode = false;
Dat.Res=[];

% If Res-structure is given as an input argument, go to "single-file" mode
if nargin == 1
  % Check that the Res structure is valid
  if isstruct(Res_in) && isfield(Res_in,'Stat') && ...
        isfield(Res_in,'FileInfo') && isfield(Res_in,'DateTime')
    Dat.SingleFileMode = true;
    Dat.Res = Res_in;
  else
    error('The Res-structure is not valid!')
	end
  
elseif nargin > 1
  error('Too many input arguments')
end


l_DrawGUI;
if Dat.SingleFileMode
  l_PreviewFile([],[])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Draw GUI objects
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_DrawGUI
  
%% Load default font and colors
  GD=aedes_gui_defaults;	

  fig_w = 1000;
  fig_h = 600;
	fig_location = aedes_dialoglocation([fig_w,fig_h]);
	fig_pos = [fig_location(1) fig_location(2) fig_w fig_h];
	
  %% Main Figure ----------------------------
  H.MAINFIG = figure('Units','Pixel', ...
                     'position',fig_pos,...
                     'Name','Results viewer', ...
                     'Numbertitle','off', ...
                     'Tag','aedes_resview_fig', ...
                     'Color',GD.col.mainfig, ...
                     'Toolbar','none', ...
                     'Menubar','none', ...
                     'DockControls','off',...
                     'renderer','painters',...
                     'KeyPressFcn','',...
                     'CloseRequestFcn',@l_Quit,...
                     'Handlevisibility','off');
  if ~GD.HG2graphics
		set(H.MAINFIG,'DoubleBuffer','on')
	end
  %% Options uipanel ----------------------------
  H.MAIN_UIPANEL = uipanel('parent',H.MAINFIG,...
                           'units','pixel',...
                           'position',...
                           [10 10 300 fig_h-2*10],...[10 fig_h-290 fig_w-2*10 280],...
						   'backgroundcolor',GD.col.frame);
  
  %% Output directory and browse button --------
  try
    out_str = getpref('Aedes','PutXLSFileDir');
  catch
    out_str = pwd;
  end
  if Dat.SingleFileMode
    edit_str = '';
    enabled = 'off';
  else
    edit_str = out_str;
    enabled = 'on';
  end
  tmp=get(H.MAIN_UIPANEL,'position');
  H.OUTDIR_EDIT=uicontrol('parent',H.MAIN_UIPANEL,...
                          'style','edit',...
                          'position',[10 tmp(4)-20-20 240 20],...
                          'backgroundcolor','w',...
                          'horizontalalign','left',...
                          'string',edit_str,...
                          'userdata',out_str,...
                          'callback',@l_CheckOutputDir,...
                          'enable',enabled);
  tmp = get(H.OUTDIR_EDIT,'position');
  H.OUTDIR_TX = uicontrol('parent',H.MAIN_UIPANEL,...
                          'style','text',...
                          'position',[tmp(1) tmp(2)+tmp(4)+2 100 12],...
                          'string','Output Folder',...
                          'Horizontalalign','left',...
                          'enable',enabled,...
						  'backgroundcolor',GD.col.frame);
  H.OUTPUT_BROWSE = uicontrol('parent',H.MAIN_UIPANEL,...
                              'style','pushbutton',...
                              'position',[tmp(1)+tmp(3)+5 tmp(2) 35 20],...
                              'string','...',...
                              'callback',@l_BrowseFolders,...
                              'enable',enabled);
                            
   %% Selected files -listbox ------------------------------------
  
  %% Disable multi-file selection in single-file mode
  if Dat.SingleFileMode
    isEnabled = 'off';
  else
    isEnabled = 'on';
  end
  btn_cdata = load('aedes_cdata.mat');
  tmp = get(H.OUTDIR_EDIT,'position');
  H.SELFILES_UIPANEL = uipanel('parent',H.MAIN_UIPANEL,...
    'units','pixel',...
    'position',[tmp(1) tmp(2)-10-185 280 ...
    185],...
    'title','Selected Files',...
    'backgroundcolor',GD.col.frame);
  tmp=get(H.SELFILES_UIPANEL,'position');
  H.SELFILES_LBOX = uicontrol('parent',H.SELFILES_UIPANEL,...
                              'units','pixel',...
                              'position',[60 25 tmp(3)-60-20 tmp(4)-45],...
                              'style','listbox',...
                              'backgroundcolor','w',...
                              'callback',@l_ListBoxCallBack,...
                              'min',0,...
                              'max',2,...
                              'Enable',isEnabled);
  H.SHOWPATH = uicontrol('parent',H.SELFILES_UIPANEL,...
                         'units','pixel',...
                         'position',[60 7 tmp(3)-60-20 15],...
                         'style','checkbox',...
                         'value',0,...
                         'string','Show full path',...
                         'callback',@l_ShowFullPath,...
						 'backgroundcolor',GD.col.frame,...
             'enable',isEnabled);
  
  H.ADDFILES = uicontrol('parent',H.SELFILES_UIPANEL,...
                         'units','pixel',...
                         'position',[10 tmp(4)-20-40 40 40],...
                         'style','pushbutton',...
                         'tooltip','Add files to list',...
                         'CData',btn_cdata.cdata.add,...
                         'callback',{@l_AddRemoveFiles,'add'},...
                         'Enable',isEnabled);
  
  H.REMFILES = uicontrol('parent',H.SELFILES_UIPANEL,...
                         'units','pixel',...
                         'position',[10 tmp(4)-20-2*40-17 40 40],...
                         'style','pushbutton',...
                         'tooltip','Remove selected files from list',...
                         'CData',btn_cdata.cdata.delete,...
                         'callback',{@l_AddRemoveFiles,'remove'},...
                         'Enable','off');
  
  H.REMALLFILES = uicontrol('parent',H.SELFILES_UIPANEL,...
                            'units','pixel',...
                            'position',[10 tmp(4)-20-3*40-2*17 40 40],...
                            'style','pushbutton',...
                            'tooltip','Remove all files from list',...
                            'CData',btn_cdata.cdata.deleteall,...
                            'callback',{@l_AddRemoveFiles,'remove_all'},...
                            'Enable','off');
  
                            
  
  %% Export file types uicontrols ----------------------------
  tmp = get(H.SELFILES_UIPANEL,'position');
  H.EXPORT_FILETYPES_GRP = uibuttongroup('parent',H.MAIN_UIPANEL,...
    'units','pixel',...
    'position',[10 tmp(2)-70-10 280 70],...
    'title','Export file types',...
    'backgroundcolor',GD.col.frame);
  
  if ispref('Aedes','ResViewerWriteXLS')
    val = getpref('Aedes','ResViewerWriteXLS');
  else
    val = 1;
    setpref('Aedes','ResViewerWriteXLS',val);
  end
  H.OUTPUTTYPE_XLS = uicontrol('parent',H.EXPORT_FILETYPES_GRP,...
    'style','checkbox',...
    'position',[10 30 150 15],...
    'value',val,...
    'string','Excel sheet (XLS)',...
    'backgroundcolor',GD.col.frame);
  if isunix
    % Excel export is only supported in Windows...
    set(H.OUTPUTTYPE_XLS,'value',0,'enable','off')
  end
  tmp = get(H.OUTPUTTYPE_XLS,'position');
  if ispref('Aedes','ResViewerWriteTXT')
    val = getpref('Aedes','ResViewerWriteTXT');
  else
    val = 1;
    setpref('Aedes','ResViewerWriteTXT',val);
  end
  H.OUTPUTTYPE_TXT = uicontrol('parent',H.EXPORT_FILETYPES_GRP,...
    'style','checkbox',...
    'position',[10 tmp(2)-15-5 150 15],...
    'value',val,...
    'string','Text file (TXT)',...
    'backgroundcolor',GD.col.frame);
  
  %% Export options ------------------------------------------
  tmp = get(H.EXPORT_FILETYPES_GRP,'position');
  H.EXPORT_OPTIONS_GRP = uibuttongroup('parent',H.MAIN_UIPANEL,...
    'units','pixel',...
    'position',[10 tmp(2)-10-220 280 220],...
    'title','Export options',...
    'backgroundcolor',GD.col.frame);
  tmp = get(H.EXPORT_OPTIONS_GRP,'position');
  H.SORTBY_TX = uicontrol('parent',H.EXPORT_OPTIONS_GRP,...
                          'units','pixel',...
                          'position',[10 tmp(4)-2*10-15 70 15],...
                          'style','text',...
                          'string','Sort by:',...
                          'horizontalalign','left',...
						  'backgroundcolor',GD.col.frame);
  tmp=get(H.SORTBY_TX,'position');
  H.SORTBY = uicontrol('parent',H.EXPORT_OPTIONS_GRP,...
                       'units','pixel',...
                       'position',[280-tmp(1)-10-100 tmp(2)+3 100 15],...
                       'style','popup',...
                       'string',{'ROI','File'},...
                       'value',1,...
                       'backgroundcolor','w',...
                       'callback',@l_PreviewFile);
  
  H.DECSEPARATOR_TX = uicontrol('parent',H.EXPORT_OPTIONS_GRP,...
                                'units','pixel',...
                                'position',[tmp(1) tmp(2)-tmp(4)-10 ...
                   130 15],...
                                'style','text',...
                                'string','Decimal separator:',...
                                'horizontalalign','left',...
								'backgroundcolor',GD.col.frame);
  tmp = get(H.DECSEPARATOR_TX,'position');
  try
	val=getpref('Aedes','ResViewerDecSep');
  catch
	val=1;
	setpref('Aedes','ResViewerDecSep',1);
  end
  H.DECSEPARATOR = uicontrol('parent',H.EXPORT_OPTIONS_GRP,...
                             'units','pixel',...
                             'position',[280-tmp(1)-10-100 tmp(2)+3 100 15],...
                             'style','popup',...
                             'string',{', (comma)','. (point)'},...
                             'value',val,...
                             'backgroundcolor','w',...
                             'callback',@l_PreviewFile);
  
  H.NUMDEC_TX =  uicontrol('parent',H.EXPORT_OPTIONS_GRP,...
                                'units','pixel',...
                                'position',[tmp(1) tmp(2)-tmp(4)-10 ...
                      150 15],...
                           'style','text',...
                           'string','Number of decimals:',...
                           'horizontalalign','left',...
						   'backgroundcolor',GD.col.frame);
  tmp = get(H.NUMDEC_TX,'position');
  if ispref('Aedes','ResViewerNumDec')
    val=getpref('Aedes','ResViewerNumDec');
    val = num2str(val);
  else
    val = '3';
  end
  H.NUMDEC = uicontrol('parent',H.EXPORT_OPTIONS_GRP,...
    'units','pixel',...
    'position',[280-tmp(1)-10-100 tmp(2)-1 100 19],...
    'style','edit',...
    'string',val,...
    'backgroundcolor','w',...
    'callback',@l_PreviewFile);
  
  % Get default print directions
  if ispref('Aedes','StatPrintDirs')
    Dat.dirs = getpref('Aedes','StatPrintDirs');
  else
    Dat.dirs = 'TXYZ';
  end
  if any(Dat.dirs=='T')
    total_val = 1;
  else
    total_val = 0;
  end
  if any(Dat.dirs=='X')
    x_val = 1;
  else
    x_val = 0;
  end
  if any(Dat.dirs=='Y')
    y_val = 1;
  else
    y_val = 0;
  end
  if any(Dat.dirs=='Z')
    z_val = 1;
  else
    z_val = 0;
  end
  if any(Dat.dirs=='V')
    v_val = 1;
  else
    v_val = 0;
  end
  H.PRINT_TOTAL = uicontrol('parent',H.EXPORT_OPTIONS_GRP,...
    'units','pixel',...
    'position',[tmp(1) tmp(2)-5-19 250 19],...
    'style','checkbox',...
    'string','Print total',...
    'backgroundcolor',GD.col.frame,...
    'value',total_val,...
    'callback',@l_PreviewFile);
  tmp = get(H.PRINT_TOTAL,'position');
  H.PRINT_XDIR = uicontrol('parent',H.EXPORT_OPTIONS_GRP,...
    'units','pixel',...
    'position',[tmp(1) tmp(2)-5-19 250 19],...
    'style','checkbox',...
    'string','Print X-Dir',...
    'backgroundcolor',GD.col.frame,...
    'value',x_val,...
    'callback',@l_PreviewFile);
  tmp = get(H.PRINT_XDIR,'position');
  H.PRINT_YDIR = uicontrol('parent',H.EXPORT_OPTIONS_GRP,...
    'units','pixel',...
    'position',[tmp(1) tmp(2)-5-19 250 19],...
    'style','checkbox',...
    'string','Print Y-Dir',...
    'backgroundcolor',GD.col.frame,...
    'value',y_val,...
    'callback',@l_PreviewFile);
  tmp = get(H.PRINT_YDIR,'position');
  H.PRINT_ZDIR = uicontrol('parent',H.EXPORT_OPTIONS_GRP,...
    'units','pixel',...
    'position',[tmp(1) tmp(2)-5-19 250 19],...
    'style','checkbox',...
    'string','Print Z-Dir',...
    'backgroundcolor',GD.col.frame,...
    'value',z_val,...
    'callback',@l_PreviewFile);
  tmp = get(H.PRINT_ZDIR,'position');
  H.PRINT_VDIR = uicontrol('parent',H.EXPORT_OPTIONS_GRP,...
    'units','pixel',...
    'position',[tmp(1) tmp(2)-5-19 250 19],...
    'style','checkbox',...
    'string','Print V-Dir',...
    'backgroundcolor',GD.col.frame,...
    'value',v_val,...
    'callback',@l_PreviewFile);
  
  %% RESTABLE -----------------------------------
  %ResTable=aedes_res2table('..\ROI_test\t1001001.res');
  tmp=get(H.MAIN_UIPANEL,'position');
  H.RESTABLE = uitable('Parent',H.MAINFIG,'position',[tmp(1)+tmp(3)+10 10 fig_w-2.5*10-tmp(3) fig_h-3*10]);
  H.RESFRAME = uicontrol('parent',H.MAINFIG,...
	'units','pixel',...
	'style','frame',...
	'position',[tmp(1)+tmp(3)+10 10 fig_w-2.5*10-tmp(3) fig_h-3*10],...
	'backgroundcolor',GD.col.frame);
  
  % Check Matlab version since uitable properties have changed in R2008a
  Dat.MatlabVersion = aedes_getmatlabversion;
  if Dat.MatlabVersion>=7.06
	set(H.RESTABLE,'Enable','on','visible','off');
  else
	set(H.RESTABLE,'Editable',false,'visible',false)
  end
  
  tmp=get(H.RESTABLE,'position');
  H.RESTABLE_TX = uicontrol('Parent',H.MAINFIG,...
    'Units','pixel',...
    'style','text',...
    'position',[tmp(1) tmp(2)+tmp(4) tmp(3) 15],...
    'string','Output preview ("N/A")',...
    'horizontalalign','left',...
    'fontweig','bold',...
    'backgroundcolor',GD.col.mainfig);
  
  % Display "Current results" text in single file mode
  if Dat.SingleFileMode
    set(H.RESTABLE_TX,'string','Output preview ("Current results")')
  end
  
  %% Export and Cancel buttons ----------------------
  tmp=get(H.SELFILES_UIPANEL,'position');
  H.EXPORT = uicontrol('parent',H.MAIN_UIPANEL,...
                       'units','pixel',...
                       'position',[tmp(1)+tmp(3)-80 5 80 25],...
                       'string','Export',...
                       'callback',@l_Export);
  tmp = get(H.EXPORT,'position');
  H.CANCEL = uicontrol('parent',H.MAIN_UIPANEL,...
                       'units','pixel',...
                       'position',[tmp(1)-80-5 5 80 25],...
                       'string','Close',...
                       'callback',@l_Quit);
    
  % Some internal variable defaults
  if isfield(Dat.Res,'Stat') && Dat.Res.Stat(1).isMixed==0
    if strcmpi(Dat.Res.FileInfo.DataFileName{1},'fid')
      [fp,fn,fe]=fileparts(Dat.Res.FileInfo.DataPathName{1}(1:end-1));
      Dat.ResFileName = fn;
    else
      Dat.ResFileName = Dat.Res.FileInfo.DataFileName{1};
    end
  else
    Dat.ResFileName = '';
  end
  Dat.OldFigSize = get(H.MAINFIG,'position');
  
  % Set resize function
  set(H.MAINFIG,'resizefcn',@l_Resize)
  
end % function l_DrawGUI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Resize Function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_Resize(h,evd)

  fig_pos = get(H.MAINFIG,'position');
  if fig_pos(4)<400 || fig_pos(3)<400
    set(H.MAINFIG,'position',Dat.OldFigSize)
    return
  end
  tmp=get(H.MAIN_UIPANEL,'position');
  set(H.MAIN_UIPANEL,'position',[tmp(1) fig_pos(4)-tmp(4)-10 tmp(3:4)]);
  %set(H.MAIN_UIPANEL,'position',...
  %                  [10 10 300 fig_pos(4)-2*10]);
  main_panel_pos = get(H.MAIN_UIPANEL,'position');
  set(H.RESTABLE,'position',[main_panel_pos(1)+main_panel_pos(3)+10 10 fig_pos(3)-2.5*10-main_panel_pos(3) fig_pos(4)-3*10]);
  set(H.RESFRAME,'position',[main_panel_pos(1)+main_panel_pos(3)+10 10 fig_pos(3)-2.5*10-main_panel_pos(3) fig_pos(4)-3*10]);
  tmp=get(H.RESFRAME,'position');
  set(H.RESTABLE_TX,'position',[tmp(1) tmp(2)+tmp(4) tmp(3) 15])
  
%   % Reposition uipanels
%   tmp=get(H.MAIN_UIPANEL,'position');
%   set(H.OUTDIR_EDIT,'position',[10 tmp(4)-20-20 240 20])
%   tmp = get(H.OUTDIR_EDIT,'position');
%   set(H.OUTDIR_TX,'position',[tmp(1) tmp(2)+tmp(4)+2 100 12])
%   set(H.OUTPUT_BROWSE,'position',[tmp(1)+tmp(3)+5 tmp(2) 35 20])
%   
%   tmp = get(H.OUTDIR_EDIT,'position');
%   set(H.SELFILES_UIPANEL,'position',[tmp(1) tmp(2)-10-185 280 185])
%   tmp = get(H.SELFILES_UIPANEL,'position');
%   set(H.EXPORT_FILETYPES_GRP,'position',[10 tmp(2)-70-10 280 70])
%   tmp = get(H.EXPORT_FILETYPES_GRP,'position');
%   set(H.EXPORT_OPTIONS_GRP,'position',[10 tmp(2)-10-220 280 220])
  
  Dat.OldFigSize = fig_pos;
    
end % function l_Resize(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ADD/REMOVE FILES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_AddRemoveFiles(h,evd,opt)
    
  %% Add files -----------------------
  if strcmpi(opt,'add')
    
    [fname,pname]=aedes_juigetfiles({'*.res;*.RES','Aedes Res-files (*.res)';...
                   '*.*','All Files (*.*)'},...
                              'Select Res-files',pwd);
    if ~iscell(fname)
      % Canceled
      return
    end
    
    %% Add selected files to the list
    FileList = get(H.SELFILES_LBOX,'string');
    FileListUD = get(H.SELFILES_LBOX,'userdata');
    if isempty(FileList)
      if get(H.SHOWPATH,'value')
        tmpstr = {};
        for ii=1:length(fname)
          tmpstr{ii} = [pname{ii},fname{ii}];
        end
        set(H.SELFILES_LBOX,'string',tmpstr);
      else
        set(H.SELFILES_LBOX,'string',fname);
      end
      FileListNewUD=cell(length(fname),2);
      FileListNewUD(:,1) = pname(:);
      FileListNewUD(:,2) = fname(:);
      set(H.SELFILES_LBOX,'userdata',FileListNewUD)
    else
      FileListNewUD=cell(length(fname),2);
      FileListNewUD(:,1) = pname(:);
      FileListNewUD(:,2) = fname(:);
      
      %% Remove dublicates
      tmp={FileListUD{all(ismember(FileListUD, ...
                                         FileListNewUD),2),:}};
      rm_ind=find(all(ismember(FileListNewUD,tmp),2));
      if ~isempty(rm_ind)
        isDublicate = true;
        DublicateList = FileListNewUD(rm_ind,:);
        FileListNewUD(rm_ind,:)=[];
      else
        isDublicate = false;
        DublicateList = {};
      end
      NewFileList = cell(size(FileListUD,1)+size(FileListNewUD,1),2);
      NewFileList(1:size(FileListUD,1),:)=FileListUD;
      NewFileList(size(FileListUD,1)+1:end,:)=FileListNewUD;
      if get(H.SHOWPATH,'value')
        FileListNewCell={}
        for ii=1:size(NewFileList,1)
          FileListNewCell{ii} = [NewFileList{:,1},NewFileList{:,2}];
        end
      else
        FileListNewCell = {NewFileList{:,2}};
      end
      set(H.SELFILES_LBOX,'string',FileListNewCell);
      set(H.SELFILES_LBOX,'userdata',NewFileList);
      set(H.SELFILES_LBOX,'value',[]);
      
      %% Warn about dublicates that were not added to the list
      if isDublicate
        warndlg({'The following dublicate files were ignored:',...
                 '',DublicateList{:,2}},...
                'Dublicate files ignored','modal');
      end

    end
    % Enable remove buttons
    set([H.REMFILES,H.REMALLFILES],'enable','on')
    
    %% Remove files ------------------
  elseif strcmpi(opt,'remove')
    
    % Ask confirmation
    resp=questdlg('Remove selected files from list?',...
                  'Remove selected files?','Yes','No','No');
    if strcmpi(resp,'No')
      return
    end
    
    val = get(H.SELFILES_LBOX,'value');
    str = get(H.SELFILES_LBOX,'string');
    % Return if files are not selected
    if isempty(val)
      return
    end
    
    NewStr = str;
    NewUD = get(H.SELFILES_LBOX,'userdata');
    NewStr(val)=[];
    NewUD(val,:)=[];
    
    set(H.SELFILES_LBOX,'value',[],'string',NewStr,...
                      'userdata',NewUD)
    if isempty(NewStr)
      % Disable remove buttons
      set([H.REMFILES,H.REMALLFILES],'enable','off')
    end

    % Check if preview needs to be cleared
    if ~isempty(Dat.ResFileName) && ~isempty(NewUD)
      [fp,fn,fe]=fileparts(Dat.ResFileName);
      fp = [fp,filesep];
      fn = [fn,fe];
      tmp1 = strcmp(fn,{NewUD{:,2}});
      tmp2 = strcmp(fp,{NewUD{:,1}});
      tmp1=tmp1(:);
      tmp2=tmp2(:);
      bool=prod(double([tmp2 tmp1]),2);
      if ~any(bool==1)
        l_ResetPreview([],[])
      end
    else
      l_ResetPreview([],[])
    end
    
    %% Remove ALL files --------------
  elseif strcmpi(opt,'remove_all')
    
    % Ask confirmation
    resp=questdlg('Remove ALL files from list?',...
                  'Remove ALL files?','Yes','No','No');
    if strcmpi(resp,'No')
      return
    end
    
    set(H.SELFILES_LBOX,'value',[],'string',{},'userdata',{})
    set([H.REMFILES,H.REMALLFILES],'enable','off')
    
    % Clear preview
    l_ResetPreview([],[])
  end
  
end % function l_AddRemoveFiles(h,


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LISTBOX CALLBACK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_ListBoxCallBack(h,evd)
  
  % Return immediately if not double-clicked
  SelType = get(H.MAINFIG,'selectiontype');
  if ~strcmpi(SelType,'open')
    return
  end
  
  % Set mouse pointer to hourglass
  set(H.MAINFIG,'pointer','watch')
  drawnow
  
  % Get file for preview
  val=get(H.SELFILES_LBOX,'value');
  if length(val)==0 || length(val)>1
    set(H.MAINFIG,'pointer','arrow')
    return
  end
  UD = get(H.SELFILES_LBOX,'userdata');
  fname = [UD{val,1},UD{val,2}];
  
  % Try to load Res-file
  try
    tmp = load(fname,'-mat');
  catch
    errordlg({'Could not open Res-file',fname,...
              '','Either the file is not a valid Res-file or it is corrupt.'},...
             'Could not open Res-file','modal');
    set(H.MAINFIG,'pointer','arrow')
    return
  end
  
  if ~isfield(tmp,'Res') || ~isfield(tmp.Res,'Stat') || ...
        ~isfield(tmp.Res,'FileInfo') || ~isfield(tmp.Res,'DateTime')
    errordlg({'Could not open Res-file',fname,...
              '','Required fields are missing from the structure!'},...
             'Could not open Res-file','modal');
    set(H.MAINFIG,'pointer','arrow')
    return
  end
  Dat.Res = tmp.Res;
  Dat.ResFileName = fname;
  
  % Preview file
  l_PreviewFile([],[]);
  set(H.MAINFIG,'pointer','arrow')
  drawnow
  
end % function l_ListBoxCallBack(h,


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PREVIEW FILE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_PreviewFile(h,evd)
  
  if isempty(Dat.Res)
    return
  end
  
  % Get export options
  sortbyroi = get(H.SORTBY,'value');
  if sortbyroi==1
    sortbyroi = true;
  else
    sortbyroi = false;
  end
  decsepval = get(H.DECSEPARATOR,'value');
  setpref('Aedes','ResViewerDecSep',decsepval)
  if decsepval==1
    decsep = ',';
  else
    decsep = '.';
  end
  numdec = get(H.NUMDEC,'string');
  numdec = str2num(numdec);
  numdec = floor(numdec);
  if isempty(numdec) || ~isreal(numdec) || numdec<0
    errordlg('The "Number of decimals" value must be a positive integer!',...
             'Error in "Number of decimals" field','modal')
    return
  end
  
  % Get print directions
  Dat.dirs = '';
  if get(H.PRINT_TOTAL,'value')
    Dat.dirs(end+1) = 'T';
  end
  if get(H.PRINT_XDIR,'value')
    Dat.dirs(end+1) = 'X';
  end
  if get(H.PRINT_YDIR,'value')
    Dat.dirs(end+1) = 'Y';
  end
  if get(H.PRINT_ZDIR,'value')
    Dat.dirs(end+1) = 'Z';
  end
  if get(H.PRINT_VDIR,'value')
    Dat.dirs(end+1) = 'V';
  end
  
  % Save number of decimals
  setpref('Aedes','ResViewerNumDec',numdec)
  
  % Save printed directions
  setpref('Aedes','StatPrintDirs',Dat.dirs)
  
  % Generate cell table from Res
  ResTable = aedes_res2table(Dat.Res,'sortbyroi',sortbyroi,...
                       'decsep',decsep,'numdec',numdec,...
                       'resfilename',Dat.ResFileName,...
                       'dirs',Dat.dirs);
  
  set(H.RESFRAME,'visible','off')
  drawnow
  drawnow
  
  % Clear restable
  if Dat.MatlabVersion>=7.06
	set(H.RESTABLE,'data',[])
	set(H.RESTABLE,'Data',ResTable,...
	  'visible','on',...
    'ColumnWidth',...
    {140,80,80,80,80,80,80},...
    'ColumnFormat',repmat({'char'},1,7),...
    'ColumnName',[])
  else
	set(H.RESTABLE,'NumRows',2,...
	  'NumColumns',2);
	drawnow
	drawnow
	tmp={'','';'',''};
	set(H.RESTABLE,'Data',tmp)
	drawnow
	drawnow
	
	% Update res table
	set(H.RESTABLE,'visible',true)
	drawnow
	drawnow
	set(H.RESTABLE,'NumRows',size(ResTable,1),...
	  'NumColumns',size(ResTable,2))
	drawnow
	drawnow
	set(H.RESTABLE,'Data',ResTable);
	drawnow
	drawnow
  end
  
  
 
  %drawnow
  if ~Dat.SingleFileMode
    set(H.RESTABLE_TX,'string',['Output preview ("' Dat.ResFileName '")'])
  end

end % function l_PreviewFiles(h,


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% RESET FILE PREVIEW
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_ResetPreview(h,evd)
  
  Dat.Res = [];
  Dat.ReFileName = '';
  if Dat.MatlabVersion>=7.06
	set(H.RESTABLE,'visible','off',...
	  'Data',[]);
  else
	set(H.RESTABLE,'visible',false,...
	  'Data',[]);
  end
  set(H.RESFRAME,'visible','on')
  set(H.RESTABLE_TX,'string','Output preview ("N/A")')
  drawnow
end % function l_ResetPreview(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EXPORT RESULTS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_Export(h,evd)
 
  % Check if there is something to export
  if ~Dat.SingleFileMode
    UD = get(H.SELFILES_LBOX,'userdata');
    if isempty(UD)
      h=errordlg('Nothing to export!','Error Exporting Files','modal');
      return
    end
  end
  
  % Get Export types
  WriteXLS = get(H.OUTPUTTYPE_XLS,'value');
  WriteTXT = get(H.OUTPUTTYPE_TXT,'value');
  if ~WriteXLS && ~WriteTXT
    h=errordlg('Export file type not selected!','Error Exporting Files','modal');
    return
  end
  
  % Save values to preferences
  setpref('Aedes','ResViewerWriteXLS',WriteXLS)
  setpref('Aedes','ResViewerWriteTXT',WriteTXT)
  
  % Get export parameters
  outdir = get(H.OUTDIR_EDIT,'string');
  if ~isempty(outdir) && outdir(end)~=filesep
    outdir = [outdir,filesep];
  end
  sortbyroi = get(H.SORTBY,'value');
  if sortbyroi==1
    sortbyroi = true;
  else
    sortbyroi = false;
  end
  decsepval = get(H.DECSEPARATOR,'value');
  if decsepval==1
    decsep = ',';
  else
    decsep = '.';
  end
  numdec = get(H.NUMDEC,'string');
  numdec = str2num(numdec);
  numdec = floor(numdec);
  if isempty(numdec) || ~isreal(numdec) || numdec<0
    errordlg('The "Number of decimals" value must be a positive integer!',...
             'Error in "Number of decimals" field','modal')
    return
  end
  
  %% Start Exporting files --------------------------
  if ~Dat.SingleFileMode
    wbar_h=aedes_wbar(0,{'Processing file','""'});
    set(get(findall(wbar_h,'type','axes'),'title'),'interpreter','none')
    nFiles = size(UD,1);
    success = true(1,nFiles);
    for ii=1:nFiles
      resfilename = [UD{ii,1},UD{ii,2}];
      [fp,fn,fe]=fileparts(UD{ii,2});
      outfilename = [outdir,fn];
      aedes_wbar(ii/nFiles,wbar_h,{'Processing file',['"',resfilename,'"']});
      
      % Open file
      try
        tmp = load([UD{ii,1},UD{ii,2}],'-mat');
      catch
        success(ii)=false;
        fprintf(1,'Could not open file "%s". Skipping...\n',resfilename)
        continue;
      end
      
      if ~isfield(tmp,'Res') || ~isfield(tmp.Res,'Stat') || ...
            ~isfield(tmp.Res,'FileInfo') || ~isfield(tmp.Res,'DateTime')
        success(ii)=false;
        fprintf(1,'Required field(s) missing! Skipping file "%s".\n',...
                resfilename)
        continue;
      end
      
      % Generate cell table
      ResTable = aedes_res2table(tmp.Res,'sortbyroi',sortbyroi,...
                           'decsep',decsep,'numdec',numdec,...
                           'resfilename',resfilename,...
                           'dirs',Dat.dirs);
      
      % Write XLS
      if WriteXLS
        bool = xlswrite([outfilename,'.xls'],ResTable);
        if ~bool
          success(ii)=false;
          fprintf(1,'Could not write "%s". Skipping file...\n',[outfilename,'.xls'])
        end
      end
      
      % Write TXT
      if WriteTXT
        bool = aedes_cellwrite(ResTable,[outfilename,'.txt'],...
                         'delimitter','tab');
        if ~bool
          success(ii)=false;
          fprintf(1,'Could not write "%s". Skipping file...\n',[outfilename,'.txt'])
        end
      end
      
    end
    close(wbar_h)
    
    if all(success)
      h=helpdlg(['Successfully exported ',num2str(nFiles),' files!'],...
                'Export Success');
    else
      SkippedFiles={};
      count=1;
      for ii=1:nFiles
        if success(ii)==0
          SkippedFiles{count} = [UD{ii,1},UD{ii,2}];
          count=count+1;
        end
      end
      h=warndlg({['Errors occurred while exporting.',...
                  ' The following files could not be exported:'],'',...
                 SkippedFiles{:}},'Errors occurred while exporting');
    end
    
    % Save export directory to preferences
    setpref('Aedes','PutXLSFileDir',outdir);
    
  else % Export in single file mode
    
    success = true;
    try
      outdir = getpref('Aedes','PutXLSFileDir');
    catch
      outdir = pwd;
    end
    
    if ~strcmp(outdir(end),filesep)
      outdir = [outdir,filesep];
    end
    
    %[fp,fn,fe]=fileparts(Dat.Res.FileInfo.DataFileName{1});
    
    % Ask for a file
    [fname,fpath,findex]=uiputfile({'*.xls;*.XLS;*.txt;*.TXT',...
                   'Excel and Text Files (*.xls,*.txt)';...
                   '*.*','All Files (*.*)'},'Export Results As',...
                                   [outdir,Dat.ResFileName]);
    if isequal(fname,0) || isequal(fpath,0)
      % Action cancelled
      return
    end
    
    % Check filename
    CheckXLS = true;
    CheckTXT = true;
    [fp,fn,fe]=fileparts(fname);
    if strcmpi(fe,'.xls')
      CheckXLS = false;
      xlsfname = [fn,fe];
      txtfname = [fn,'.txt'];
    elseif strcmpi(fe,'.txt')
      CheckTXT = false;
      xlsfname = [fn,'.xls'];
      txtfname = [fn,fe];
    else
      xlsfname = [fn,fe,'.xls'];
      txtfname = [fn,fe,'.txt'];
    end
    
    % Check if files exist
    if WriteXLS && CheckXLS
      if exist(xlsfname,'file')==2
        resp = questdlg({'File',['"',xlsfname,'"'],'already exists. Overwrite?'},...
          'Overwrite file?','Yes','No','No');
        if isempty(resp) || strcmpi(resp,'No')
          return
        end
      end
    end
    if WriteTXT && CheckTXT
      if exist(txtfname,'file')==2
        resp = questdlg({'File',['"',txtfname,'"'],'already exists. Overwrite?'},...
          'Overwrite file?','Yes','No','No');
        if isempty(resp) || strcmpi(resp,'No')
          return
        end
      end
    end
    
    % Generate cell table
    ResTable = aedes_res2table(Dat.Res,'sortbyroi',sortbyroi,...
                         'decsep',decsep,'numdec',numdec,...
                         'resfilename',Dat.ResFileName,...
                         'dirs',Dat.dirs);
    
    % Show aedes_calc_wait
    [h,txh]=aedes_calc_wait('');
    
    % Write XLS
    if WriteXLS
      set(txh,'string',{'Writing Excel XLS file',[fpath,xlsfname]})
      drawnow
      bool = xlswrite([fpath,xlsfname],ResTable);
      if ~bool
        success=false;
        fprintf(1,'Could not write "%s". Skipping file...\n',[fpath,xlsfname])
      end
    end
      
    % Write TXT
    if WriteTXT
      set(txh,'string',{'Writing TXT file',[fpath,txtfname]})
      drawnow
      bool = aedes_cellwrite(ResTable,[fpath,txtfname],...
                       'delimitter','tab');
      if ~bool
        success=false;
        fprintf(1,'Could not write "%s". Skipping file...\n',[fpath,txtfname])
      end
    end
    delete(h);
    
    if ~success
      h=warndlg('Errors occurred while exporting. Export failed.','Export failed!');
    end
    
    % Save export directory to preferences
    setpref('Aedes','PutXLSFileDir',fpath);
    
  end
  
  
  
end % function l_Export(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BROWSE FOLDERS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_BrowseFolders(h,evd)
  
  % Get currently selected folder from the edit box
  OldDir=get(H.OUTDIR_EDIT,'string');
  if ~isdir(OldDir)
    OldDir = pwd;
  end
  
  % Show folder selection dialog
  NewDir = uigetdir(OldDir,'Select output folder');
  if all(NewDir==0)
    % Action canceled
    return
  end
  
  % Update output directory edit box
  set(H.OUTDIR_EDIT,'string',NewDir,...
                    'userdata',NewDir)
  
end % function l_BrowseFolders(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CHECK OUTPUT DIRECTORY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_CheckOutputDir(h,evd) % ---------------------------
  
  output_dir = get(H.OUTDIR_EDIT,'string');
  if ~isdir(output_dir)
    h=warndlg({'The selected folder',['"' output_dir '"'],...
               'doesn''t exist.'},'Folder doesn''t exist','modal');
    set(H.OUTDIR_EDIT,'string',get(H.OUTDIR_EDIT,'userdata'))
  else
    set(H.OUTDIR_EDIT,'userdata',output_dir)
  end
  
end % function l_CheckOutputDir(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% QUIT AEDES_RESVIEWER
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_Quit(h,evd) % --------------------------------------
  
  % Clear public variables
  fig_h = H.MAINFIG;
  clear H Dat
  
  % Close aedes_resviewer window
  delete(fig_h)
  
end % function l_Quit(h,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Show/hide full path
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_ShowFullPath(h,evd)
  
  UD = get(H.SELFILES_LBOX,'userdata');
  if isempty(UD)
    return
  end
  
  % Construct new string to the file listbox
  NewStr = {};
  if get(H.SHOWPATH,'value')
    for ii=1:size(UD,1)
      NewStr{ii} = [UD{ii,1},UD{ii,2}];
    end
  else
    NewStr = {UD{:,2}};
  end
  set(H.SELFILES_LBOX,'string',NewStr)
  
end % function l_ShowFullPath(h,

end % function aedes_resviewer(Res_in)

