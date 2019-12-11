function [done,msg] = aedes_saveres(DATA,ROI,filename,varargin)
% AEDES_SAVERES - Save Aedes results (statistics and ROIs)
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


done=false;
msg='';

%% Default values for properties
SaveType = 'all';      % Save ROIs and statistics by default
FileName = 'Untitled'; % Filename(s) Untitled.res and Untitled.roi
FilePath = pwd;        % Use working directory by default
ShowWbar = true;       % Show waitbar by default
RotateFlip = [];
ConfirmOverwrite = true;


%% Parse input arguments
if nargin<3
  msg = 'Too few input arguments';
  return
end

%% Parse varargin
for ii=1:2:length(varargin)
  switch lower(varargin{ii})
   case {'wbar','waitbar'}
    if ischar(varargin{ii+1})
      if strcmpi(varargin{ii+1},'on')
        ShowWbar = true;
      else
        ShowWbar = false;
      end
    else
      if varargin{ii+1}==1
        ShowWbar = true;
      else
        ShowWbar = false;
      end
    end
   case 'savetype'
    SaveType = varargin{ii+1};
    
   case 'rotateflip'
    RotateFlip = varargin{ii+1};
    
    case 'confirmoverwrite'
    ConfirmOverwrite = varargin{ii+1};
    
   otherwise
    done=false;
    msg=sprintf('Unknown parameter "%s"!',varargin{ii});
    return
  end
end

[fp,fn,fe]=fileparts(filename);
if any(strcmpi(fe,{'.roi','.res'}))
  FileName = fn;
else
  FileName = [fn,fe];
end
if ~isempty(fp)
  FilePath = [fp,filesep];
else
  FilePath = [pwd,filesep];
end

if isstruct(DATA)
  DATA = {DATA};
elseif ~iscell(DATA)
  msg='Input argument "DATA" is not valid';
  return
end

if isempty(ROI) || ~isstruct(ROI)
  msg='Input argument "ROI" is not valid';
  return
end

%% Save results and/or ROIs
switch SaveType
 case 'all' % Save both ROI and RES --------------------
  
  %% Check if files to be written already exist
  tmp=dir(FilePath);
  fnames = {tmp(~[tmp(:).isdir]).name};
  if any(strcmpi(fnames,[FileName,'.roi']))
    RoiFileExists = true;
  else
    RoiFileExists = false;
  end
  if any(strcmpi(fnames,[FileName,'.res']))
    ResFileExists = true;
  else
    ResFileExists = false;
  end
  
  if ResFileExists && ConfirmOverwrite
    resp=questdlg({['"',FilePath,FileName,'.res" already exists.'],...
                  'Overwrite?'},'Overwrite File?',...
                  'Yes','No','No');
    if strcmpi(resp,'No')
      msg='Overwrite cancel';
      return
    end
  end
  if RoiFileExists && ConfirmOverwrite
    resp=questdlg({['"',FilePath,FileName,'.roi" already exists.'],...
                  'Overwrite?'},'Overwrite File?',...
                  'Yes','No','No');
    if strcmpi(resp,'No')
      msg='Overwrite cancel';
      return
    end
  end
  
  if ShowWbar
    [h,txh]=aedes_calc_wait({'Saving ROI(s)...',...
                      ['(',FilePath,FileName,'.roi)']});
  end
  
  [ok,msg]=l_SaveRoi(DATA,ROI,FileName,FilePath,RotateFlip);
  if ~ok
    if ShowWbar
      delete(h)
    end
    return
  end
  set(txh,'string',{'Saving ROI Statistics...',...
                   ['(',FilePath,FileName,'.res)']})
  [ok,msg]=l_SaveRes(DATA,ROI,FileName,FilePath);
  if ~ok
    if ShowWbar
      delete(h)
    end
    return
  end
  if ShowWbar
    delete(h)
  end
  
 case 'roi' % Save only ROI ----------------------------
  
  %% Check if files to be written already exist
  tmp=dir(FilePath);
  fnames = {tmp(~[tmp(:).isdir]).name};
  if any(strcmpi(fnames,[FileName,'.roi']))
    RoiFileExists = true;
  else
    RoiFileExists = false;
  end
  
  if RoiFileExists && ConfirmOverwrite
    resp=questdlg({['"',FilePath,FileName,'.roi" already exists.'],...
                  'Overwrite?'},'Overwrite File?',...
                  'Yes','No','No');
    if strcmpi(resp,'No')
      msg='Overwrite cancel';
      return
    end
  end
  
  if ShowWbar
    [h,txh]=aedes_calc_wait({'Saving ROI(s)...',...
                      ['(',FilePath,FileName,'.roi)']});
  end
  
  [ok,msg]=l_SaveRoi(DATA,ROI,FileName,FilePath,RotateFlip);
  if ~ok
    if ShowWbar
      delete(h)
    end
    return
  end
  if ShowWbar
    delete(h)
  end
  
 case 'res' % Save only RES ----------------------------
  
  %% Check if files to be written already exist
  tmp=dir(FilePath);
  fnames = {tmp(~[tmp(:).isdir]).name};
  if any(strcmpi(fnames,[FileName,'.res']))
    ResFileExists = true;
  else
    ResFileExists = false;
  end
  
  if ResFileExists && ConfirmOverwrite
    resp=questdlg({['"',FilePath,FileName,'.res" already exists.'],...
                  'Overwrite?'},'Overwrite File?',...
                  'Yes','No','No');
    if strcmpi(resp,'No')
      msg='Overwrite cancel';
      return
    end
  end
  
  if ShowWbar
    [h,txh]=aedes_calc_wait({'Saving ROI Statistics...',...
                       ['(',FilePath,FileName,'.res)']});
  end  
  
  [ok,msg]=l_SaveRes(DATA,ROI,FileName,FilePath);
  if ~ok
    if ShowWbar
      delete(h)
    end
    return
  end
  if ShowWbar
    delete(h)
  end
  
 otherwise
  msg='Unknown save type';
  return
end

% No problems if this point is reached...
done=true;
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save ROIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ok,msg]=l_SaveRoi(DATA,ROI,FileName,FilePath,rotflip)

ok = false;
msg='';

% Generate timestamp
DateTime = datestr(now);

% Append file information to the save ROI file
for ii=1:length(DATA)
  FileInfo.DataFileName{ii} = DATA{ii}.HDR.fname;
  FileInfo.DataPathName{ii} = DATA{ii}.HDR.fpath;
end

% Append rotate and flip information
if iscell(rotflip)
  RotateFlip = rotflip;
else
  RotateFlip.Rotate = zeros(1,length(DATA));
  RotateFlip.Flip = zeros(1,length(DATA));
  if ~isempty(rotflip)
	RotateFlip.Rotate = rotflip(1,:);
	RotateFlip.Flip = rotflip(2,:);
  end
end

% Save ROI(s)
try
  save([FilePath,FileName,'.roi'],'ROI','DateTime','FileInfo','RotateFlip','-mat')
catch
  msg={'Could not save ROI(s). Following error was returned:',lasterr};
  return
end

% All went well
ok=true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save RES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ok,msg]=l_SaveRes(DATA,ROI,FileName,FilePath)

ok = false;
msg = '';


try
  Res = [];
  % Calculate results
  Res = aedes_roi_stats(DATA,ROI);
  if isempty(Res)
    msg = 'Error while calculating statistics. Could not save results.';
    return
  end
  
  % Generate timestamp
  %DateTime = datestr(now);
  
  % Generate file info
  %for ii=1:length(DATA)
  %  FileInfo.DataFileName{ii} = DATA{ii}.HDR.fname;
  %  FileInfo.DataPathName{ii} = DATA{ii}.HDR.fpath;
  %end
  %Res.DateTime = DateTime;
  %Res.FileInfo = FileInfo;
  %Res.Stat = Stat;
  
  % Save results
  save([FilePath,FileName,'.res'],'Res','-mat')
catch
  msg = {'Could not save results. Following error was returned:', ...
         lasterr};
  return
end

% All went well
ok = true;
