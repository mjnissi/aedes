function aedes_readfidprefs()
% AEDES_READFIDPREFS - GUI for editing preferences for AEDES and
% AEDES_DATA_READ for VNMR format files
%   
%
% Synopsis: 
%        aedes_readfidprefs;
%
% Description:
%
% Examples:
%
% See also:
%        AEDES, AEDES_DATA_READ

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


done=false;
msg='';


Dat = [];

%% Try to read defaults from preferences
%% Get defaults for Return
if ispref('Aedes','ReadfidReturn')
  if getpref('Aedes','ReadfidReturn')==1
    Dat.ReturnKSpace  = false;
    Dat.ReturnFTData  = true;
  elseif getpref('Aedes','ReadfidReturn')==2
      Dat.ReturnKSpace  = true;
      Dat.ReturnFTData  = false;
  elseif getpref('Aedes','ReadfidReturn')==3
    Dat.ReturnKSpace  = true;
    Dat.ReturnFTData  = true;
  else
    Dat.ReturnKSpace  = false;
    Dat.ReturnFTData  = true;
  end
else
  Dat.ReturnKSpace  = false;
  Dat.ReturnFTData  = true;
  setpref('Aedes','ReadfidReturn',1)
end

%% Get defaults for DC correction
if ispref('Aedes','ReadfidDCcorrection')
  if getpref('Aedes','ReadfidDCcorrection')
    Dat.DCcorrection  = true;
  else
    Dat.DCcorrection  = false;
  end
else
  Dat.DCcorrection  = false;
  setpref('Aedes','ReadfidDCcorrection',Dat.DCcorrection)
end

%% Get defaults for Zeropadding
if ispref('Aedes','ReadfidZeropadding')
  if getpref('Aedes','ReadfidZeropadding')==0
    Dat.ZeroPadding = 0;
  elseif getpref('Aedes','ReadfidZeropadding')==1
    Dat.ZeroPadding = 1;
  elseif getpref('Aedes','ReadfidZeropadding')==2
    Dat.ZeroPadding = 2;
  else
    Dat.ZeroPadding = 2;
  end
else
  Dat.ZeroPadding = 2;
  setpref('Aedes','ReadfidZeropadding',Dat.ZeroPadding)
end

%% Get defaults for Sorting
if ispref('Aedes','ReadfidSorting')
  if getpref('Aedes','ReadfidSorting')
    Dat.Sorting = true;
  else
    Dat.Sorting = false;
  end
else
  Dat.Sorting = true;
  setpref('Aedes','ReadfidSorting',Dat.Sorting)
end

%% Get defaults for FastRead
if ispref('Aedes','ReadfidFastRead')
  if getpref('Aedes','ReadfidFastRead')
    Dat.FastDataRead = true;
  else
    Dat.FastDataRead = false;
  end
else
  Dat.FastDataRead = true;
  setpref('Aedes','ReadfidFastRead',Dat.FastDataRead)
end

%% Get defaults for Precision
if ispref('Aedes','ReadfidPrecision')
  if strcmpi(getpref('Aedes','ReadfidPrecision'),'single')
	Dat.Precision = 'single';
  else
	Dat.Precision = 'double';
  end
else
  Dat.Precision = 'single';
  setpref('Aedes','ReadfidPrecision',Dat.Precision)
end

% %% Get defaults for Reorienting EPI images
% if ispref('Aedes','ReadfidReorientEPI')
%   Dat.ReorientEPI = getpref('Aedes','ReadfidReorientEPI');
% else
%   Dat.ReorientEPI = 'off';
%   setpref('Aedes','ReadfidReorientEPI','off')
% end

%% Get defaults for Reorienting images according to procpar.orient
if ispref('Aedes','ReadfidOrientImages')
  Dat.OrientImages = getpref('Aedes','ReadfidOrientImages');
else
  Dat.OrientImages = 'on';
  setpref('Aedes','ReadfidOrientImages','on')
end

%% Get defaults for Removing phase image from EPI data
if ispref('Aedes','ReadfidRemoveEPIphaseIm')
  Dat.RemoveEPIphaseIm = getpref('Aedes','ReadfidRemoveEPIphaseIm');
else
  Dat.RemoveEPIphaseIm = 'off';
  setpref('Aedes','ReadfidRemoveEPIphaseIm','off')
end

%% Get defaults for read function
if ispref('Aedes','VnmrUseOldReadFcn')
  if getpref('Aedes','VnmrUseOldReadFcn')
    Dat.VnmrUseOldReadFcn = true;
  else
    Dat.VnmrUseOldReadFcn = false;
  end
else
  Dat.VnmrUseOldReadFcn = false;
  setpref('Aedes','VnmrUseOldReadFcn',Dat.VnmrUseOldReadFcn)
end

%% Load default font and colors
%FigColor=get(0,'DefaultUicontrolBackgroundcolor');
GD=aedes_gui_defaults;
%GD.col.mainfig = FigColor;
fig_h = 300;
fig_w = 270;
fig_location = aedes_dialoglocation([fig_w,fig_h]);
fig_pos = [fig_location(1) fig_location(2) fig_w fig_h];

%% The main figure
fh = figure('position',fig_pos,...
            'Units','Pixel', ...
            'Name','Edit VNMR preferences', ...
            'Numbertitle','off', ...
            'Tag','readfid_prefs', ...
            'Color',GD.col.mainfig, ...
            'Toolbar','none', ...
            'Menubar','none', ... 
            'DockControls','off',...
            'renderer','painters',...
            'KeyPressFcn','',...
            'resize','off',...
            'windowstyle','normal');
if ~GD.HG2graphics
	set(fh,'DoubleBuffer','on')
end

%% Main Uipanel
uipanel_h = uipanel('parent',fh,...
                    'Units','pixel',...
                    'position',[5 40 fig_w-10 fig_h-45],...
                    'title','Preferences for reading VNMR files',...
                    'fontweig','bold',...
					'backgroundcolor',GD.col.frame);

%% OK button
ok_btn = uicontrol('parent',fh,...
                   'units','pixel',...
                   'position',[fig_w-170 5 80 30],...
                   'string','OK',...
                   'style','pushbutton',...
                   'callback',@l_ChangePrefs);

%% Cancel button
tmp=get(ok_btn,'position');
cancel_btn = uicontrol('parent',fh,...
                       'units','pixel',...
                       'position',[tmp(1)+tmp(3)+5 tmp(2:4)],...
                       'string','Cancel',...
                       'style','pushbutton',...
                       'callback','delete(gcbf)');

%% Returned data popup
if Dat.ReturnKSpace==true & Dat.ReturnFTData==false
  val=2;
elseif Dat.ReturnKSpace==false & Dat.ReturnFTData==true
  val=1;
else
  val=3;
end
tmp=get(uipanel_h,'position');
return_tx = uicontrol('parent',uipanel_h,...
                      'units','pixel',...
                      'position',[10 tmp(4)-50 100 20],...
                      'style','text',...
                      'horizontalalign','left',...
                      'string','Returned data',...
					  'backgroundcolor',GD.col.frame,...
					  'fontsize',GD.text_fs);
tmp=get(return_tx,'position');
return_popup = uicontrol('parent',uipanel_h,...
                         'units','pixel',...
                         'position',[tmp(1)+tmp(3)+1 tmp(2)+3 140 20],...
                         'style','popup',...
                         'backgroundcolor','w',...
                         'string',{'FT-Data','K-Space','FT-Data and K-Space'},...
                         'value',val);

%% DC correction
if Dat.DCcorrection
  val=1;
else
  val=2;
end
tmp=get(return_tx,'position');
dc_tx = uicontrol('parent',uipanel_h,...
                  'units','pixel',...
                  'position',[10 tmp(2)-tmp(4)-5 tmp(3) 20],...
                  'style','text',...
                  'horizontalalign','left',...
                  'string','DC correction',...
				  'backgroundcolor',GD.col.frame);
tmp=get(dc_tx,'position');
dc_popup = uicontrol('parent',uipanel_h,...
                     'units','pixel',...
                     'position',[tmp(1)+tmp(3)+1 tmp(2)+3 140 20],...
                     'style','popup',...
                     'backgroundcolor','w',...
                     'string',{'On','Off'},...
                     'value',val);

%% Zeropadding
if Dat.ZeroPadding==1
  val=1;
elseif Dat.ZeroPadding==0
  val=2;
elseif Dat.ZeroPadding==2
  val=3;
end
tmp=get(dc_tx,'position');
zeropadding_tx = uicontrol('parent',uipanel_h,...
                           'units','pixel',...
                           'position',[10 tmp(2)-tmp(4)-5 tmp(3) 20],...
                           'style','text',...
                           'horizontalalign','left',...
                           'string','Zeropadding',...
						   'backgroundcolor',GD.col.frame);
tmp=get(zeropadding_tx,'position');
zeropadding_popup = uicontrol('parent',uipanel_h,...
                              'units','pixel',...
                              'position',[tmp(1)+tmp(3)+1 tmp(2)+3 140 20],...
                              'style','popup',...
                              'backgroundcolor','w',...
                              'string',{'On','Off','Auto'},...
                              'value',val);

%% Sorting
if Dat.Sorting
  val=1;
else
  val=2;
end
tmp=get(zeropadding_tx,'position');
sorting_tx = uicontrol('parent',uipanel_h,...
                       'units','pixel',...
                       'position',[10 tmp(2)-tmp(4)-5 tmp(3) 20],...
                       'style','text',...
                       'horizontalalign','left',...
                       'string','Sorting',...
					   'backgroundcolor',GD.col.frame);
tmp=get(sorting_tx,'position');
sorting_popup = uicontrol('parent',uipanel_h,...
                          'units','pixel',...
                          'position',[tmp(1)+tmp(3)+1 tmp(2)+3 140 20],...
                          'style','popup',...
                          'backgroundcolor','w',...
                          'string',{'On','Off'},...
                          'value',val);
%% Fast Read
if Dat.FastDataRead
  val=1;
else
  val=2;
end
tmp=get(sorting_tx,'position');
fastread_tx = uicontrol('parent',uipanel_h,...
  'units','pixel',...
  'position',[10 tmp(2)-tmp(4)-5 tmp(3) 20],...
  'style','text',...
  'horizontalalign','left',...
  'string','FastRead',...
  'backgroundcolor',GD.col.frame);
tmp=get(fastread_tx,'position');
fastread_popup = uicontrol('parent',uipanel_h,...
  'units','pixel',...
  'position',[tmp(1)+tmp(3)+1 tmp(2)+3 140 20],...
  'style','popup',...
  'backgroundcolor','w',...
  'string',{'On','Off'},...
  'value',val);

%% Precision
if strcmpi(Dat.Precision,'single')
  val=1;
else
  val=2;
end
tmp=get(fastread_tx,'position');
precision_tx = uicontrol('parent',uipanel_h,...
  'units','pixel',...
  'position',[10 tmp(2)-tmp(4)-5 tmp(3) 20],...
  'style','text',...
  'horizontalalign','left',...
  'string','Precision',...
  'backgroundcolor',GD.col.frame);
tmp=get(precision_tx,'position');
precision_popup = uicontrol('parent',uipanel_h,...
  'units','pixel',...
  'position',[tmp(1)+tmp(3)+1 tmp(2)+3 140 20],...
  'style','popup',...
  'backgroundcolor','w',...
  'string',{'single','double'},...
  'value',val);

%% Default read function
if Dat.VnmrUseOldReadFcn
  val=1;
else
  val=2;
end
tmp=get(precision_tx,'position');
readfcn_tx = uicontrol('parent',uipanel_h,...
  'units','pixel',...
  'position',[10 tmp(2)-tmp(4)-5 tmp(3) 20],...
  'style','text',...
  'horizontalalign','left',...
  'string','Read Fcn',...
  'backgroundcolor',GD.col.frame);
tmp=get(readfcn_tx,'position');
readfcn_popup = uicontrol('parent',uipanel_h,...
  'units','pixel',...
  'position',[tmp(1)+tmp(3)+1 tmp(2)+3 140 20],...
  'style','popup',...
  'backgroundcolor','w',...
  'string',{'readfid (old)','readvnmr'},...
  'value',val);


%% Reorient images according to procpar.orient
if strcmpi(Dat.OrientImages,'on')
  val=1;
else
  val=0;
end
tmp=get(readfcn_tx,'position');
orient_images_chbox = uicontrol('parent',uipanel_h,...
  'units','pixel',...
  'position',[10 tmp(2)-tmp(4)-5 240 20],...
  'style','checkbox',...
  'horizontalalign','left',...
  'backgroundcolor',GD.col.frame,...
  'string','Orient images using PROCPAR.orient',...
  'value',val);

%% Remove phase image frim EPI data
if strcmpi(Dat.RemoveEPIphaseIm,'on')
  val=1;
else
  val=0;
end
tmp=get(orient_images_chbox,'position');
removeepiphaseim_chbox = uicontrol('parent',uipanel_h,...
  'units','pixel',...
  'position',[10 tmp(2)-tmp(4)-5 240 20],...
  'style','checkbox',...
  'horizontalalign','left',...
  'backgroundcolor',GD.col.frame,...
  'string','Remove phase image from EPI data',...
  'value',val);

% Store handles to a structure
H.fh = fh;
H.return_h = return_popup;
H.dccorr_h = dc_popup;
H.zeropadding_h = zeropadding_popup;
H.sorting_h = sorting_popup;
H.fastread_h = fastread_popup;
H.precision_h = precision_popup;
H.readfcn_h = readfcn_popup;
H.orient_images_h = orient_images_chbox;
H.removeepiphaseim_h = removeepiphaseim_chbox;
set(ok_btn,'userdata',H)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OK button press callback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_ChangePrefs(h,evd)

H = get(h,'userdata');

try
% Get popup values
return_val = get(H.return_h,'value');
dccorr_val = get(H.dccorr_h,'value');
zeropadding_val = get(H.zeropadding_h,'value');
sorting_val = get(H.sorting_h,'value');
fastread_val = get(H.fastread_h,'value');
precision_val = get(H.precision_h,'value');
readfcn_val = get(H.readfcn_h,'value');
orient_images_val = get(H.orient_images_h,'value');
removeepiphaseim_val = get(H.removeepiphaseim_h,'value');

% Set preferences
%% -------------------------------------------------------------------

%% Return value
setpref('Aedes','ReadfidReturn',return_val)

%% DC correction
if dccorr_val==1
  setpref('Aedes','ReadfidDCcorrection',true)
else
  setpref('Aedes','ReadfidDCcorrection',false)
end

%% Zeropadding
if zeropadding_val==1
  setpref('Aedes','ReadfidZeropadding',1)
elseif zeropadding_val==2
  setpref('Aedes','ReadfidZeropadding',0)
else
  setpref('Aedes','ReadfidZeropadding',2)
end

%% Sorting
if sorting_val==1
  setpref('Aedes','ReadfidSorting',true)
else
  setpref('Aedes','ReadfidSorting',false)
end

%% FastRead
if fastread_val==1
  setpref('Aedes','ReadfidFastRead',true)
else
  setpref('Aedes','ReadfidFastRead',false)
end

%% Precision
if precision_val==1
  setpref('Aedes','ReadfidPrecision','single')
else
  setpref('Aedes','ReadfidPrecision','double')
end

%% Read fcn
if readfcn_val==1
  setpref('Aedes','VnmrUseOldReadFcn',true)
else
  setpref('Aedes','VnmrUseOldReadFcn',false)
end

%% Orient images
if orient_images_val==1
  setpref('Aedes','ReadfidOrientImages','on')
else
  setpref('Aedes','ReadfidOrientImages','off')
end

%% Remove phase image from EPI
if removeepiphaseim_val==1
  setpref('Aedes','ReadfidRemoveEPIphaseIm','on')
else
  setpref('Aedes','ReadfidRemoveEPIphaseIm','off')
end

catch
  hh=errordlg('An error occurred! Could not save preferences.',...
              'Could not save preferences','modal');
end

delete(H.fh)
