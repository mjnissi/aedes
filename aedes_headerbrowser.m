function aedes_headerbrowser(inStruct)
% AEDES_HEADERBROWSER - Browse and search file header information or any
%                       Matlab structure
%   
%
% Synopsis:
%       aedes_headerbrowser(DATA) % Open the file format specific header
%
%       aedes_headerbrowser(STRUCT) % Browse a generic Matlab structure
%
% Description:
%       
%
% Examples:
%
% See also:
%       AEDES, AEDES_READPROCPAR, AEDES_READ_NIFTI

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
H = [];    % Handles structure
Dat = [];  % Internal data structure

% Do some error checking
if nargin < 1
  error('Invalid number of input arguments!')
end

% Check that the input argument is structure
if ~isstruct(inStruct)
  error('Input argument must be structure!')
end

% Get the saved Favourite list
try
  Dat.FavList = getpref('Aedes','HeaderBrowserFavList');
catch
  Dat.FavList.nifti = {};
  Dat.FavList.vnmr = {};
  Dat.FavList.dcm = {};
	Dat.FavList.bruker_reco = {};
	Dat.FavList.bruker_raw = {};
  Dat.FavList.genheader = {};
  Dat.FavList.generic = {};
end

% Check if input is a Aedes data structure
if isfield(inStruct,'HDR') && isfield(inStruct,'FTDATA')
  if isfield(inStruct,'DataFormat')
    if any(strcmpi(inStruct.DataFormat,{'NIfTI(1)','NIfTI(2)','Analyze75'}))
      % NIfTI and Analyze75 files
      Dat.HeaderFileName = [inStruct.HDR.fpath,inStruct.HDR.fname];
      Dat.FileFormatName = 'NIfTI/Analyze';
      try
        Dat.FavouriteList = Dat.FavList.nifti;
      catch
        Dat.FavList.nifti = {};
        Dat.FavouriteList = Dat.FavList.nifti;
      end
      inStruct=inStruct.HDR.FileHeader;
      
    elseif strcmpi(inStruct.DataFormat,'VNMR')
      % Varian VNMR files
      Dat.HeaderFileName = [inStruct.HDR.fpath,'procpar'];
      Dat.FileFormatName = 'VNMR Procpar';
      try
        Dat.FavouriteList = Dat.FavList.vnmr;
      catch
        Dat.FavList.vnmr = {};
        Dat.FavouriteList = Dat.FavList.vnmr;
      end
      inStruct=inStruct.PROCPAR;
      
    elseif strcmpi(inStruct.DataFormat,'DCM')
      % DICOM files
      Dat.HeaderFileName = [inStruct.HDR.fpath,inStruct.HDR.fname];
      Dat.FileFormatName = 'DICOM';
      try
        Dat.FavouriteList = Dat.FavList.dcm;
      catch
        Dat.FavList.dcm = {};
        Dat.FavouriteList = Dat.FavList.dcm;
      end
      inStruct=inStruct.HDR.FileHeader;
      
		elseif strcmpi(inStruct.DataFormat,'bruker_reco')
			% Bruker 2DSEQ files
			Dat.HeaderFileName = [inStruct.HDR.fpath,inStruct.HDR.fname];
			Dat.FileFormatName = 'Bruker 2DSEQ';
			try
				Dat.FavouriteList = Dat.FavList.bruker_reco;
			catch
				Dat.FavList.bruker_reco = {};
				Dat.FavouriteList = Dat.FavList.bruker_reco;
			end
			inStruct=inStruct.HDR.FileHeader;
			
		elseif strcmpi(inStruct.DataFormat,'bruker_raw')
			% Bruker Raw FID files
			Dat.HeaderFileName = [inStruct.HDR.fpath,inStruct.HDR.fname];
			Dat.FileFormatName = 'Bruker FID';
			try
				Dat.FavouriteList = Dat.FavList.bruker_raw;
			catch
				Dat.FavList.bruker_raw = {};
				Dat.FavouriteList = Dat.FavList.bruker_raw;
			end
			inStruct=inStruct.HDR.FileHeader;
			
    else
      % Generic header
      Dat.HeaderFileName = [inStruct.HDR.fpath,inStruct.HDR.fname];
      Dat.FileFormatName = 'generic header';
      try
        Dat.FavouriteList = Dat.FavList.genheader;
      catch
        Dat.FavList.genheader = {};
        Dat.FavouriteList = Dat.FavList.genheader;
      end
      inStruct=inStruct.HDR.FileHeader;
      
    end
  else
    % Generic header
    inStruct=inStruct.HDR.FileHeader;
    Dat.FileFormatName = 'generic header';
    try
      Dat.FavouriteList = Dat.FavList.genheader;
    catch
      Dat.FavList.genheader = {};
      Dat.FavouriteList = Dat.FavList.genheader;
    end
    Dat.HeaderFileName = [inStruct.HDR.fpath,inStruct.HDR.fname];
  end
else
  % Generic Matlab structure
  Dat.HeaderFileName = '';
  Dat.FileFormatName = 'generic structure';
  try
    Dat.FavouriteList = Dat.FavList.generic;
  catch
    Dat.FavList.generic = {};
    Dat.FavouriteList = Dat.FavList.generic;
  end
end

% Initialize the cell array for uitable
Dat.FullTable = {};
Dat.StructMaxDepth = 10; % Limit recursion to something reasonable...
Dat.FieldSeparator = '\';
Dat.TableFavSeparator = '<html><b>=============</b></html>';
Dat.SelectedCells = [];

% Draw main GUI
[H,Dat]=l_DrawGui(Dat);
l_ParseHeaderStructure(inStruct,'',0);
Dat.DataTable = Dat.FullTable;
l_SearchParams([],[]);

% Update Data Table
l_UpdateDataTable([],[]);

% Get java handle to the table for Matlab R2007a and older
if Dat.MatlabVersion<=7.05
  drawnow, pause(0.01) % Force refresh
  
  % Set Table to non-editable
  set(H.DataTable,'Editable',false)
  
  % Get java handle
  H.DataTableJavaHandle = handle(H.DataTable.getTable,'callbackproperties');
  
  % Set KeyPressedCallback
  set(H.DataTableJavaHandle.KeyPressedCallback,@l_KeyPressFcn)
end

% Set initial focus to search editbox
drawnow
uicontrol(H.SearchEdit);

% For debuging
%assignin('base','H',H) % debug

%% ------------------------------------------
  function [H,Dat]=l_DrawGui(Dat)

    % Load default font and colors
    Dat.GD=aedes_gui_defaults;
    
    % Get Matlab version
    Dat.MatlabVersion = aedes_getmatlabversion;

    fig_w = 600;
    fig_h = 700;
		fig_location = aedes_dialoglocation([fig_w,fig_h]);
		fig_pos = [fig_location(1) fig_location(2) fig_w fig_h];
		
    % Main Figure ----------------------------
    H.Fig = figure('Units','Pixel', ...
      'position',fig_pos,...
      'Name',['File Header Browser - ',Dat.FileFormatName], ...
      'Numbertitle','off', ...
      'Tag','header_browser_fig', ...
      'Color',Dat.GD.col.mainfig, ...
      'Toolbar','none', ...
      'Menubar','none', ...
      'DockControls','off',...
      'renderer','painters',...
      'resize','off',...
      'Handlevisibility','off',...
      'KeyPressFcn',@l_KeyPressFcn,...
      'CloseReq',@l_CloseFig);
    if ~Dat.GD.HG2graphics
			set(H.Fig,'DoubleBuffer','on')
		end
    % File name
    fname_txh = 0;
    if ~isempty(Dat.HeaderFileName)
      fname_txh = 70;
      H.FileNameUipanel = uipanel('parent',H.Fig,...
        'units','pixel',...
        'position',[10 fig_h-fname_txh fig_w-20 fname_txh-10],...
        'title','File Name',...
        'backgroundcolor',Dat.GD.col.mainfig);
      
      tmp_pos=get(H.FileNameUipanel,'position');
      H.FileNameText = uicontrol('parent',H.FileNameUipanel,...
        'style','text',...
        'string',Dat.HeaderFileName,...
        'units','pixel',...
        'position',[10 5 tmp_pos(3)-20 fname_txh-35],...
        'horizontalalign','left',...
        'backgroundcolor',Dat.GD.col.mainfig);
    end
    
    % Search uipanel ----------------------
    up_height = 80;
    H.SearchUipanel = uipanel('parent',H.Fig,...
      'units','pixel',...
      'position',[10 fig_h-up_height-10-fname_txh fig_w-20 up_height],...
      'title','Search Options',...
      'backgroundcolor',Dat.GD.col.mainfig);

    % Search text
    pos=get(H.SearchUipanel,'position');
    %H.SearchText = uicontrol('parent',H.SearchUipanel,...
    %  'style','text',...
    %  'string','Parameter:',...
    %  'units','pixel',...
    %  'position',[10 pos(4)-45 100 19],...
    %  'horizontalalign','left',...
    %  'backgroundcolor',Dat.GD.col.mainfig);
    try
      val=getpref('Aedes','HeaderBrowserSearchWhere');
    catch
      val=1;
    end
    H.SearchWherePopup = uicontrol('parent',H.SearchUipanel,...
      'style','popup',...
      'string',{'Parameter',...
      'Value'},...
      'units','pixel',...
      'position',[10 pos(4)-45 100 22],...
      'value',val,...
      'tag','',...
      'backgroundcolor',Dat.GD.col.popup,...
      'callback',@l_SearchParams,...
      'KeyPressFcn',@l_KeyPressFcn);
    
    % Search type popupmenu
    try
      val = getpref('Aedes','HeaderBrowserSearchHow');
    catch
      val = 1;
    end
    tmp_pos = get(H.SearchWherePopup,'position');
    %tmp_pos = get(H.SearchText,'position');
    H.SearchHowPopup = uicontrol('parent',H.SearchUipanel,...
      'style','popup',...
      'string',{'contains...',...
      'is exactly...','begins with...',...
      'ends with...',...
      'matches regexp...'},...
      'units','pixel',...
      'position',[tmp_pos(1)+tmp_pos(3)+5 tmp_pos(2) ...
      150 22],...
      'value',val,...
      'tag','',...
      'backgroundcolor',Dat.GD.col.popup,...
      'callback',@l_SearchParams,...
      'KeyPressFcn',@l_KeyPressFcn);

    % Search editbox
    tmp_pos = get(H.SearchHowPopup,'position');
    H.SearchEdit = uicontrol('parent',H.SearchUipanel,...
      'style','edit',...
      'units','pixel',...
      'position',[tmp_pos(1)+tmp_pos(3)+5 tmp_pos(2) pos(3)-tmp_pos(1)-tmp_pos(3)-15 tmp_pos(4)],...
      'Callback',@l_SearchParams,...
      'backgroundcolor',Dat.GD.col.edit,...
      'KeyPressFcn',@l_SearchParams);
    
    % Case sensitive
    try
      val = getpref('Aedes','HeaderBrowserCaseSensitivity');
    catch
      val = 1;
    end
    tmp_pos = get(H.SearchEdit,'position');
    H.SearchCaseSens = uicontrol('parent',H.SearchUipanel,...
      'style','checkbox',...
      'units','pixel',...
      'string','Case sensitive',...
      'position',[tmp_pos(1) tmp_pos(2)-25 tmp_pos(3) tmp_pos(4)],...
      'Callback',@l_SearchParams,...
      'value',val,...
      'backgroundcolor',Dat.GD.col.mainfig,...
      'KeyPressFcn',@l_KeyPressFcn);
    
    % --------------------------------------
    % Favourite params list
    favup_height = 150;
    tmp_pos = get(H.SearchUipanel,'position');
    H.FavUipanel = uipanel('parent',H.Fig,...
      'units','pixel',...
      'position',[10 fig_h-up_height-favup_height-fname_txh-15 fig_w-20 favup_height],...
      'title',['Favourite parameters for ',Dat.FileFormatName],...
      'backgroundcolor',Dat.GD.col.mainfig);
    
    tmp_pos = get(H.FavUipanel,'position');
    InfoTx = {['The parameter names entered here will be searched ',...
      'automatically. The matching parameters are placed on top of the parameter list ',...
      'in boldface. This automatic search is case sensitive and only returns exact matches.']};
    H.FavInfoText = uicontrol('parent',H.FavUipanel,...
      'style','text',...
      'units','pixel',...
      'string','',...
      'position',[10 10 250 120],...
      'horizontalalign','left',...
      'fontsize',10,...
      'backgroundcolor',[0.85 0.85 0.85]);%Dat.GD.col.mainfig);
    InfoTx=textwrap(H.FavInfoText,InfoTx);
    set(H.FavInfoText,'string',InfoTx)
    
    tmp_pos = get(H.FavInfoText,'position');
    H.FavAdd = uicontrol('parent',H.FavUipanel,...
      'style','pushbutton',...
      'units','pixels',...
      'position',[tmp_pos(1)+tmp_pos(3)+5 ...
      tmp_pos(2)+tmp_pos(4)-25 70 25],...
      'String','Add',...
      'callback',{@l_AddFavourite,[]},...
      'KeyPressFcn',@l_KeyPressFcn);
    
    H.FavRem = uicontrol('parent',H.FavUipanel,...
      'style','pushbutton',...
      'units','pixels',...
      'position',[tmp_pos(1)+tmp_pos(3)+5 ...
      tmp_pos(2)+tmp_pos(4)-2*25-5 70 25],...
      'String','Remove',...
      'callback',@l_RemFavourite,...
      'KeyPressFcn',@l_KeyPressFcn);
    
    tmp_pos = get(H.FavInfoText,'position');
    H.FavListBox = uicontrol('parent',H.FavUipanel,...
      'style','listbox',...
      'units','pixels',...
      'min',1,...
      'max',3,...
      'position',[tmp_pos(1)+tmp_pos(3)+80 ...
      tmp_pos(2) fig_w-30-tmp_pos(3)-90 tmp_pos(4)],...
      'String',Dat.FavouriteList,...
      'value',[],...
      'backgroundcolor','w',...
      'KeyPressFcn',@l_KeyPressFcn);
    
    % Close & Export to workspace btns -------
    tmp_pos = get(H.FavUipanel,'position');
    btn_w = (fig_w-20-2*5)/3;
    %btn_w = 200;
    btn_h = 25;
    H.ExportBtn = uicontrol('parent',H.Fig,...
      'style','pushbutton',...
      'units','pixels',...
      'position',[tmp_pos(1) ...
      tmp_pos(2)-btn_h-5 btn_w btn_h],...
      'string','Export to workspace',...
      'callback',@l_ExportToWorkspace,...
      'KeyPressFcn',@l_KeyPressFcn);
    
    %btn_w = 200;
    btn_h = 25;
    H.AddParamsToFavBtn = uicontrol('parent',H.Fig,...
      'style','pushbutton',...
      'units','pixels',...
      'position',[tmp_pos(1)+btn_w+5 tmp_pos(2)-btn_h-5 btn_w btn_h],...
      'string','Add to favourites',...
      'enable','off',...
      'callback',@l_AddToFavBtnCallback,...
      'KeyPressFcn',@l_KeyPressFcn);
    if Dat.MatlabVersion<=7.05
      % With Matlab R2007b and older this has to be handled differently
      set(H.AddParamsToFavBtn,'enable','on');
    end
    
    %btn_w = 200;
    btn_h = 25;
    H.CloseBtn = uicontrol('parent',H.Fig,...
      'style','pushbutton',...
      'units','pixels',...
      'position',[tmp_pos(1)+tmp_pos(3)-btn_w ...
      tmp_pos(2)-btn_h-5 btn_w btn_h],...
      'string','Close',...
      'callback',@l_CloseFig,...
      'KeyPressFcn',@l_KeyPressFcn);
    
    % Data Table -----------------------------
    pos=get(H.FavUipanel,'position');
    if Dat.MatlabVersion>=7.06
      H.DataTable = uitable('Parent',H.Fig,...
        'Data',{'loading...','loading...'},...
        'units','pixel',...
        'position',[10 10 fig_w-20 pos(2)-20-btn_h-5],...
        'KeyPressFcn',@l_KeyPressFcn,...
        'CellSelectionCallback',@l_DataTableCellSelectionFcn,...
        'ColumnFormat',{'char','char'},...
        'ColumnWidth',{280,280},...
        'RowName',[],...
        'ColumnName',{'Parameter','Value(s)'});
    else
      H.DataTable = uitable('Parent',H.Fig,...
        'position',[10 10 fig_w-20 pos(2)-20-btn_h-5],...
        'Data',{'loading...','loading...'},...
        'ColumnWidth',265,...
        'ColumnNames',{'Parameter','Value(s)'});
    end
    
    % If other header browser windows are currently open, the figure name
    % has to be unique for finding the correct javaFrame.
    figs=findall(0,'tag','header_browser_fig');
    if length(figs)>1
      set(H.Fig,'Name',['File Header Browser - ',Dat.FileFormatName,' (',num2str(H.Fig),')']);
      drawnow,pause(0.02);
    end
    H.SearchEditJavaHandle=l_GetEditboxJavaHandle(H.Fig);
    if isempty(H.SearchEditJavaHandle)
      % Fall back to basic matlab code -> no search as you type... :-(
      Dat.SearchEditUseJava = false;
      set(H.SearchEdit,'KeyPressFcn','')
    else
      Dat.SearchEditUseJava = true;
    end
    set(H.Fig,'Name',['File Header Browser - ',Dat.FileFormatName]); % Set Figure name back to original...
    
  end

  function l_ParseHeaderStructure(inStruct,parent,level)
    % This subfunction goes throught the input structure recursively and
    % constructs the n-by-2 cell array for the uitable.
    
    % Get fieldnames
    fldnames = fieldnames(inStruct);
    fldnames = sort(fldnames);
    
    % Loop over the field names
    for ii=1:length(fldnames)
      cVal = inStruct.(fldnames{ii});
      
      % Structures
      if isstruct(cVal)
        if level>=Dat.StructMaxDepth
          % Don't go beyond the StructMaxDepth level
          Dat.FullTable{end,2} = regexprep(regexprep(mat2str(size(cVal)),'\s','x'),'\]$',' struct]');
        else
          % Recursively go through the next level
          if isempty(parent) % First level
            l_ParseHeaderStructure(cVal,fldnames{ii},level+1);
          else
            l_ParseHeaderStructure(cVal,[parent,Dat.FieldSeparator,fldnames{ii}],level+1);
          end
        end
      else
        if isempty(parent)
          Dat.FullTable{end+1,1} = fldnames{ii};
        else
          Dat.FullTable{end+1,1} = [parent,Dat.FieldSeparator,fldnames{ii}];
        end
        
        % Numeric values
        if isnumeric(cVal)
          if ndims(cVal)<=2 %|| any(size(cVal)==1)
            Dat.FullTable{end,2} = mat2str(cVal);
          else
            Dat.FullTable{end,2} = regexprep(regexprep(mat2str(size(cVal)),'\s','x'),'\]$',' matrix]');
          end
          
          % Logical arrays
        elseif islogical(cVal)
          if ndims(cVal)<=2 || any(size(cVal)==1)
            Dat.FullTable{end,2} = num2str(cVal);
          else
            Dat.FullTable{end,2} = regexprep(regexprep(mat2str(size(cVal)),'\s','x'),'\]$',' logical]');
          end
          
          % Character strings
        elseif ischar(cVal)
          Dat.FullTable{end,2} = ['''',cVal,''''];
          
          % Cell arrays
        elseif iscell(cVal)
          str = '{';
          for kk=1:length(cVal)
            if kk==1
              sep='';
            else
              sep=', ';
            end
            if ischar(cVal{kk})
              str = [str,sep,'''',cVal{kk},''''];
            elseif isnumeric(cVal{kk}) || islogical(cVal{kk})
              if ~all(size(cVal{kk})==1)
                tmp_str = regexprep(regexprep(mat2str(size(cVal{kk})),'\s','x'),'\]$',' matrix]');
                str = [str,sep,tmp_str];
              else
                str = [str,sep,mat2str(cVal{kk})];
              end
            elseif iscell(cVal{kk})
              % Don't handle cell arrays inside cell arrays
              tmp_str = regexprep(regexprep(mat2str(size(cVal{kk})),'\s','x'),'\]$',' cell array]');
              str = [str,sep,tmp_str];
            elseif isstruct(cVal{kk})
              % Don't handle structures inside cell arrays
              tmp_str = regexprep(regexprep(mat2str(size(cVal{kk})),'\s','x'),'\]$',' struct]');
              str = [str,sep,tmp_str];
            else
              tmp_str = regexprep(regexprep(mat2str(size(cVal{kk})),'\s','x'),'\]$',[' ',class(cVal{kk}),']']);
              str = [str,sep,tmp_str];
            end
          end
          str(end+1) = '}';
          Dat.FullTable{end,2} = str;
          
          % Unsupported class/object
        else
          tmp_str = regexprep(regexprep(mat2str(size(cVal)),'\s','x'),'\]$',[' ',class(cVal),']']);
          Dat.FullTable{end,2} = tmp_str;
        end
      end
    end
    
    % Put favourites on top
    %Dat.DataTable = l_PutFavouritesOnTop(Dat.FullTable);
    
  end

%% ---------------------------------------------
  function l_SearchParams(fh,evd)

    if ~isempty(evd) && isfield(evd,'Key') && ...
        strcmpi(evd.Key,'escape')
      l_KeyPressFcn(fh,evd);
      return
    end
    
    % Get the search string
    if Dat.SearchEditUseJava
      drawnow
      SearchStr = get(H.SearchEditJavaHandle,'Text');
    else
      SearchStr = get(H.SearchEdit,'string');
    end
    
    % Get case sensitivity
    CaseSens = get(H.SearchCaseSens,'value');
    if CaseSens==0
      RegExpOpt = 'ignorecase';
    else
      RegExpOpt = 'matchcase';
    end
    
    % Search from parameters or values...?
    col_ind=get(H.SearchWherePopup,'value');
    
    if isempty(SearchStr)
      % If the search string is empty, show full table
      Dat.DataTable = l_PutFavouritesOnTop(Dat.FullTable);
    else
      OrigSearchStr = SearchStr;
      
      % Escape the regexp special characters
      SearchStr=regexptranslate('escape',SearchStr);
      
      if col_ind==1
        SearchStr=strrep(SearchStr,' ','|');
      end
      switch get(H.SearchHowPopup,'value')
        case 1
          % Search contains
          RegExpStr = SearchStr;
        case 2
          % Search is exactly
          SearchStr=regexprep(SearchStr,'\|','\$\|\^');
          RegExpStr = ['^',SearchStr,'$'];
        case 3
          % Search begins with
          SearchStr=regexprep(SearchStr,'\|','\|\^');
          RegExpStr = ['^',SearchStr];
        case 4
          % Search ends with
          SearchStr=regexprep(SearchStr,'\|','\$\|');
          RegExpStr = [SearchStr,'$'];
        case 5
          % Search matches regexp
          RegExpStr = OrigSearchStr;
          SearchStr = OrigSearchStr;
      end
      
      % Regexp the full table parameters
      ind=find(~cellfun('isempty',regexp(Dat.FullTable(:,col_ind),RegExpStr,RegExpOpt)));
      
      % Update DataTable
      if isempty(ind)
        if Dat.MatlabVersion>=7.06
          Dat.DataTable = {'<html><b>No matches!</b></html>',''};
          %set(H.SearchEdit,'string',['<html><font color="red">',SearchStr,'</font></html>']);
        else
          Dat.DataTable = {'<html><b>No matches!</b></html>','<html><b>No matches!</b></html>';...
            '=============','============='};
        end
      else
        Dat.DataTable = Dat.FullTable(ind,:);
        Dat.DataTable(:,col_ind)=regexprep(Dat.DataTable(:,col_ind),['(',SearchStr,')'],['<b>$1</b>'],RegExpOpt);
        Dat.DataTable(:,col_ind)=regexprep(Dat.DataTable(:,col_ind),'(^.*$)','<html>$1</html>',RegExpOpt);
      end
    end

    % Update Data Table
    l_UpdateDataTable([],[])
    
    
    
  end

%% ---------------------------------------------
  function l_UpdateDataTable(fh,evd)
    
    if Dat.MatlabVersion>=7.06
      set(H.DataTable,'Data',Dat.DataTable)
    else
      if size(Dat.DataTable,1)==1
        Dat.DataTable(end+1,:)={Dat.TableFavSeparator,Dat.TableFavSeparator};
      end
      set(H.DataTable,'Data',Dat.DataTable)
    end
    
  end

%% ---------------------------------------------
  function l_AddFavourite(h,evd,param)
    
    if isempty(param)
      % Prompt for a favourite search string
      resp = aedes_inputdlg('Search text to be added:');
      if isempty(resp) || isempty(resp{1})
        return
      end
      param = resp;
    end
    
    params_added = 0;
    for ii=1:length(param)
      
      % Check if parameter is separator
      if strcmp(param{ii},Dat.TableFavSeparator)
        continue
      end
      
      % Check if the parameter is already in favourites
      if any(strcmp(param{ii},Dat.FavouriteList))
        % Parameter already in the favourite list
        %if ~isempty(h) && h==H.FavAdd
          % Put the favourite on top of the list, if already exists
          ind = find(strcmp(param{ii},Dat.FavouriteList));
          inds = [ind setdiff(1:length(Dat.FavouriteList),ind)];
          Dat.FavouriteList = Dat.FavouriteList(inds);
          params_added = params_added + 1;
        %else
          % Add search string to list
          %Dat.FavouriteList{end+1} = param{ii};
          %params_added = params_added + 1;
        %  continue
        %end
      else
        % Add search string to list
        Dat.FavouriteList{end+1} = param{ii};
        params_added = params_added + 1;
      end
    end
    
    if params_added>0
      set(H.FavListBox,'string',Dat.FavouriteList,...
        'value',[]);
    
      % Run search
      l_SearchParams([],[]);
    end
    
  end

%% ---------------------------------------------
  function l_AddToFavBtnCallback(h,evd)
    
    % Matlab R2008a onwards...
    if Dat.MatlabVersion>=7.06
      if isempty(Dat.SelectedCells)
        % This should not happen, but return immediately if it does...
        return
      end
      
      ind = find(Dat.SelectedCells(:,2)==1);
      if isempty(ind)
        return
      end
      
      % Get index to favourites separator
      ind=Dat.SelectedCells(ind,1);
      sep_ind=find(strcmp(Dat.TableFavSeparator,Dat.DataTable(:,1)));
      
      % Remove separator from the list
      if ~isempty(sep_ind)
        ind=setdiff(ind,sep_ind);
        if isempty(ind)
          return
        end
      end
      %if ~isempty(sep_ind)
      %  ind(ind<=sep_ind)=[];
      %  if isempty(ind)
      %    return
      %  end
      %end
    else
      % Matlab R2007a and older
      selCols = H.DataTableJavaHandle.SelectedColumns;
      selRows = H.DataTableJavaHandle.SelectedRows;
      
      if any(selCols==0)
        ind = selRows+1;
      else
        return
      end
      
      % Get index to favourites separator
      sep_ind=find(strcmp(Dat.TableFavSeparator,Dat.DataTable(:,1)));
      
      % Remove separator from the list
      if ~isempty(sep_ind)
        ind=setdiff(ind,sep_ind);
        if isempty(ind)
          return
        end
      end
    end
    
    % Get the parameters to be added
    params=Dat.DataTable(ind,1);
    
    % Remove html tags from the strings
    params=regexprep(params,'<html>|</html>|<b>|</b>','');
    
    % Add to favourites
    l_AddFavourite([],[],params);
    
  end

%% ---------------------------------------------
  function l_RemFavourite(fh,evd)
    
    % Get selected indices
    ind = get(H.FavListBox,'value');
    if isempty(ind)
      return
    end
    
    % Ask for confirmation
    resp = questdlg('Remove selected favourites?',...
      'Remove favourites','Yes','No','No');
    if strcmpi(resp,'No')
      return
    end
    Dat.FavouriteList(ind)=[];
    set(H.FavListBox,'string',Dat.FavouriteList,...
      'value',[]);
    
    % Run search to update changes
    l_SearchParams([],[]);
    
  end


%% ---------------------------------------------
  function TableOut = l_PutFavouritesOnTop(TableIn)
    
    FavouriteList = Dat.FavouriteList;%get(H.FavListBox,'string');
    
    % Put favourite list matches on top of the search list and bold face
    ind = [];
    for ii=1:length(FavouriteList)
      SearchStr = regexptranslate('escape',FavouriteList{ii});
      tmp_ind=find(~cellfun('isempty',regexp(TableIn(:,1),['^',SearchStr,'$'])));
      if ~isempty(tmp_ind)
        ind = [ind tmp_ind];
      end
      if isempty(tmp_ind)
        % If favourite is not found, color it red in the listbox
        ListboxStr=get(H.FavListBox,'string');
        ListboxStr = regexprep(ListboxStr,['(^',SearchStr,'$)'],...
          '<html><font color="red"><b>$1 (not found)</b></font></html>');
        set(H.FavListBox,'string',ListboxStr)
      end
    end
    if ~isempty(ind)
      ind2 = setdiff(1:size(TableIn,1),ind);
      tmpFavList = TableIn(ind,:);
      tmpFavList(:,1)=regexprep(tmpFavList(:,1),'(.*)',...
        '<html><b>$1</b></html>');
      TableOut = cat(1,...
        tmpFavList,...
        {Dat.TableFavSeparator,Dat.TableFavSeparator},...
        TableIn(ind2,:));
    else
      TableOut = TableIn;
    end
  end

%% --------------------------------------
  function l_ExportToWorkspace(fh,evd)
    
    varname = 'header_struct';
    assignin('base',varname,inStruct);
    msgstr = sprintf('Structure exported to variable "%s".',varname);
    
    % Print notification to workspace
    fprintf(1,['AEDES_HEADERBROWSER: ',msgstr,'\n']);
    
    % Nag with a dialog...
    helpdlg(msgstr,'Structure exported!')
    
  end

%% --------------------------------------------
  function l_KeyPressFcn(fh,evd)
    
    % Call l_CloseFig if ESC is pressed
    if ~isempty(evd) && isfield(evd,'Key') && ...
        strcmpi(evd.Key,'escape')
      l_CloseFig([],[])
      
    elseif ~isempty(evd) && ~isempty(evd.Modifier)
      if ( strcmpi(evd.Modifier{1},'control') && strcmpi(evd.Key,'l') ) || ...
          ( strcmpi(evd.Modifier{1},'alt') && strcmpi(evd.Key,'d') )
        % Change focus to the search editbox
        % using CTRL-L or ALT-D
        uicontrol(H.SearchEdit);
      end
    end
    
  end

%% ---------------------------------------------
  function l_DataTableCellSelectionFcn(h,evd)

    if isfield(evd,'Indices')
      Dat.SelectedCells = evd.Indices;
      if isempty(Dat.SelectedCells)
        set(H.AddParamsToFavBtn,'enable','off');
      else
        set(H.AddParamsToFavBtn,'enable','on');
      end
    else
      Dat.SelectedCells = [];
      set(H.AddParamsToFavBtn,'enable','off');
    end
    
  end
%% --------------------------------------------
  function l_CloseFig(fh,evd)
    
    % Try to save preferences to disk before exiting
    try
      CaseSens = get(H.SearchCaseSens,'value');
      setpref('Aedes','HeaderBrowserCaseSensitivity',CaseSens);
      
      valhow = get(H.SearchHowPopup,'value');
      setpref('Aedes','HeaderBrowserSearchHow',valhow);
      
      valwhere = get(H.SearchWherePopup,'value');
      setpref('Aedes','HeaderBrowserSearchWhere',valwhere);
      
      % Update the corresponding FavList
      if strcmp(Dat.FileFormatName,'NIfTI/Analyze')
        Dat.FavList.nifti = Dat.FavouriteList;
      elseif strcmp(Dat.FileFormatName,'VNMR Procpar')
        Dat.FavList.vnmr = Dat.FavouriteList;
      elseif strcmp(Dat.FileFormatName,'DICOM')
        Dat.FavList.dcm = Dat.FavouriteList;
      elseif strcmp(Dat.FileFormatName,'generic header')
        Dat.FavList.genheader = Dat.FavouriteList;
			elseif strcmp(Dat.FileFormatName,'generic structure')
        Dat.FavList.generic = Dat.FavouriteList;
      end
      setpref('Aedes','HeaderBrowserFavList',Dat.FavList);
    catch
      disp('Warning: Could not save AEDES_HEADERBROWSER preferences!');
    end
    delete(H.Fig);
    
  end

%% -----------------------------------------------
  function jhEditBox=l_GetEditboxJavaHandle(hFig)
    % Parts of this code are originally written by Yair Altman, whose
    % wizardry in Java programming in Matlab has been of great help in many
    % occasions. Please see http://undocumentedmatlab.com and the
    % comp.soft-sys.matlab news group for further information. The
    % following license notice is for findjobj.m (by Yair Altman, available
    % in the Matlab File Exchange,
    % http://www.mathworks.com/matlabcentral/fileexchange/14317) from which
    % I have used code snippets in this subfunction.
    %
    %
    % Copyright (c) 2009, Yair Altman All rights reserved.
    %
    % Redistribution and use in source and binary forms, with or without
    % modification, are permitted provided that the following conditions
    % are met:
    %
    %   * Redistributions of source code must retain the above copyright
    %     notice, this list of conditions and the following disclaimer.
    %   * Redistributions in binary form must reproduce the above copyright
    %     notice, this list of conditions and the following disclaimer in
    %     the documentation and/or other materials provided with the
    %     distribution
    %
    % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    % "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    % LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    % A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    % OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    % SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    % LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    % DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    % THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    % (INCLUDING NEGLIGENCE OR OTHERWISE) RISING IN ANY WAY OUT OF THE USE
    % OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    
    % Ensure that objects have been drawn 
    drawnow;
    pause(0.01);
    
    % Get figure's root pane
     try
       figName = get(hFig,'Name');
       mde = com.mathworks.mde.desk.MLDesktop.getInstance;
       jFigPanel = mde.getClient(figName);
       jRootPane = jFigPanel;
       jRootPane = jFigPanel.getRootPane;
     catch
       try
         warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');  % R2008b compatibility
         jFrame = get(hFig,'JavaFrame');
         jFigPanel = get(jFrame,'FigurePanelContainer');
         jRootPane = jFigPanel;
         jRootPane = jFigPanel.getComponent(0).getRootPane;
       catch
         % Never mind
       end
     end
     try
      % If invalid RootPane - try another method...
      warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');  % R2008b compatibility
      jFrame = get(hFig,'JavaFrame');
      jAxisComponent = get(jFrame,'AxisComponent');
      jRootPane = jAxisComponent.getParent.getParent.getRootPane;
     catch
      % Never mind
     end
     %jRootPane=handle(jRootPane,'callbackproperties');
     
     % Find java handle to Search editbox recursively
     Dat.javaHandleFound = false;
     jhEditBox=l_FindJavaHandle(jRootPane);
  end

% Recursion function for finding java handle for search editbox
  function jhout=l_FindJavaHandle(jhin)
    jhout = handle([]);
    jhin=handle(jhin,'callbackproperties');
    
    if ~isempty(regexp(char(jhin.toString),'EditTextPeer\$hgTextField'))
      jhout = jhin;
      Dat.javaHandleFound=true;
      return;
    else
      try
        compCount = jhin.getComponentCount;
        for ii=0:compCount-1
          childComp = jhin.getComponent(ii);
          if ~Dat.javaHandleFound
            jhout=l_FindJavaHandle(childComp);
          else
            break;
          end
        end
      catch
        % Nothing to do here...
      end
    end
    
  end
     

end

  
