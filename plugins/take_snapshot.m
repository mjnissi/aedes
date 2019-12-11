function take_snapshot(DATA,ROI,AddInfo)
% TAKE_SNAPSHOT - A snapshot plugin for Aedes
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


% Make sure that Aedes window is on top
figure(AddInfo.hFigure)

% Make sure that the window has time to refresh
pause(0.1),drawnow

if AddInfo.isDataMixed || AddInfo.AxView~=0
  % Take the snapshot from X-Dir axis
  SnapShot = getframe(AddInfo.hAxes(AddInfo.AxView));
else
  % Take snapshots of all 3 axis
  tmp1=getframe(AddInfo.hAxes(1));
  tmp2=getframe(AddInfo.hAxes(2));
  tmp3=getframe(AddInfo.hAxes(3));
  %tmp4 = zeros([size(tmp1.cdata,2)+size(tmp2.cdata,2)-size(tmp3.cdata,2) ...
  %              size(tmp3.cdata,1) 3],'uint8');
  tmp4 = zeros([size(tmp3.cdata,1) size(tmp1.cdata,2)+size(tmp2.cdata,2)-size(tmp3.cdata,2) ...
                 3],'uint8');
  %keyboard
  SnapShot.cdata = [tmp1.cdata,tmp2.cdata;tmp3.cdata,tmp4];
end

% Get default filepath
try
  filepath=getpref('Aedes','PutSnapshotPath');
catch
  filepath = [pwd,filesep];
end

% Ask where to save the image
[fname,fpath,findex]=uiputfile({'*.png;*.PNG','PNG Image Files (*.png)';...
                    '*.tif;*.tiff;*.TIF;*.TIFF','TIFF Image Files (*.tif)';...
                    '*jpg;*.JPG;*.jpeg;*.JPEG','JPEG Image Files (*.jpg)';...
                    '*.bmp;*.BMP','BMP Windows Bitmap Files (*.bmp)'},...
                               'Save Snapshot As',...
                               [filepath,'Snapshot']);
if isequal(fname,0) || isequal(fpath,0)
  return
end

% Save selected folder to preferences
setpref('Aedes','PutSnapshotPath',fpath)

% Divide the filepath into parts
[fp,fn,fe]=fileparts([fpath,fname]);

warnOverwrite = false;


% Check file format
if findex==1
  fformat = 'PNG';
  fext = '.png';
  if ~strcmpi(fe,'.png')
    filename = [fpath,fname,'.png'];
    warnOverwrite = true;
  else
    filename = [fpath,fname];
  end
elseif findex==2
  fformat = 'TIFF';
  fext = '.tif';
  if ~any(strcmpi(fe,{'.tiff','.tif'}))
    filename = [fpath,fname,'.tif'];
    warnOverwrite = true;
  else
    filename = [fpath,fname];
  end
elseif findex==3
  fformat = 'JPEG';
  fext = '.jpg';
  if ~any(strcmpi(fe,{'.jpg','.jpeg'}))
    filename = [fpath,fname,'.jpg'];
    warnOverwrite = true;
  else
    filename = [fpath,fname];
  end
elseif findex==4
  fformat = 'BMP';
  fext = '.bmp';
  if ~strcmpi(fe,'.bmp')
    filename = [fpath,fname,'.bmp'];
    warnOverwrite = true;
  else
    filename = [fpath,fname];
  end
end

% Warn possible overwrite
d=dir(fpath);
files_in_dir = {d(~[d(:).isdir]).name};
if any(strcmpi([fname,fext],files_in_dir)) && warnOverwrite
  resp = questdlg({[filename,'already exists.'],'Do you want to replace it?'},...
                  'Overwrite Existing File?','Yes','No','No');
  if strcmpi(resp,'No')
    return
  end
end

% Save the snapshot into a image file
try
  switch fformat
    case 'PNG'
      imwrite(SnapShot.cdata,filename,fformat)
    case 'TIFF'
      imwrite(SnapShot.cdata,filename,fformat)
    case 'JPEG'
      imwrite(SnapShot.cdata,filename,fformat,'Quality',85)
    case 'BMP'
      imwrite(SnapShot.cdata,filename,fformat)
  end
catch
  % If error occurs, show an error dialog
  h=errordlg({'Could not write image file.','',lasterr},...
             'Error','modal');
end
