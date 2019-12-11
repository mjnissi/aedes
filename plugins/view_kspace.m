function view_kspace(DATA,ROI,AddInfo)
% VIEW_KSPACE - Aedes plugin for viewing k-space of vnmr files
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



%% Check that the current file is vnmr data
isVnmr = false;
if AddInfo.isDataMixed
  currentSlice = AddInfo.CurrentSlice;
else
  currentSlice = 1;
end
if isfield(DATA{currentSlice},'DataFormat') && ...
    strcmpi(DATA{currentSlice}.DataFormat,'vnmr')
  isVnmr = true;
elseif isfield(DATA{currentSlice},'HDR') && ...
    strcmpi(DATA{currentSlice}.HDR.fname,'fid')
  isVnmr = true;
end
if ~isVnmr
  hh=errordlg('The current file does not contain Varian VNMR data!',...
              'Error','modal');
  return
end

%% Check if the DATA structure already contains k-space. Otherwise the
%% k-space data has to be read.
if ~isfield(DATA{currentSlice},'KSPACE') || ...
    isempty(DATA{currentSlice}.KSPACE)
  
  % Get file name and path
  if isfield(DATA{currentSlice},'HDR') && ...
      isfield(DATA{currentSlice}.HDR,'fname') && ...
      isfield(DATA{currentSlice}.HDR,'fpath')
    
    fname = DATA{currentSlice}.HDR.fname;
    fpath = DATA{currentSlice}.HDR.fpath;
  else
    hh=errordlg('The current file does not contain Varian VNMR data!',...
                'Error','modal');
    return
  end
  
  resp=questdlg({'The k-space data has to be read from the file',...
                sprintf('"%s"',[fpath,fname]),'',...
                'Do you want to continue?'},...
                'Read k-space from file?','OK','Abort','OK');
  if strcmpi(resp,'Abort')
    return
  end
  
  %% Check if file exists
  if exist([fpath,fname],'file')~=2
    hh=errordlg({'Cannot read file',...
                ['"',fpath,fname,'"'],'',...
                'The file does not exist!'},...
                'File does not exist!','modal');
    return
  end
  
  %% Read K-space
  [tmp,msg]=aedes_readfid([fpath,fname],'return',2,...
    'FastRead','on','OrientImages','on');
  if isempty(tmp)
    if iscell(msg)
      hh==errordlg({'Error while reading file.',...
                    '',msg{:}},'Error reading file',...
                   'modal');
      return
    else
      hh==errordlg({'Error while reading file.',...
                    '',msg},'Error reading file',...
                   'modal');
      return
    end
  end
  DATA_NEW = tmp;
  DATA_NEW.FTDATA = abs(DATA_NEW.KSPACE);
  DATA_NEW.FTDATA(:,:,:,2) = real(DATA_NEW.KSPACE);
  DATA_NEW.FTDATA(:,:,:,3) = imag(DATA_NEW.KSPACE);
  DATA_NEW.KSPACE = [];
else
  DATA_NEW = DATA{currentSlice};
  DATA_NEW.FTDATA = abs(DATA_NEW.KSPACE);
  DATA_NEW.FTDATA(:,:,:,2) = real(DATA_NEW.KSPACE);
  DATA_NEW.FTDATA(:,:,:,3) = imag(DATA_NEW.KSPACE);
  DATA_NEW.KSPACE=[];
end

%% Open the new DATA-structure in a new Aedes window
aedes(DATA_NEW)

