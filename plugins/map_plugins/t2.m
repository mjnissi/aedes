function t2(DATA,ROI,AddInfo)
% This Aedes plugin calculates T2 map

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

fit_vals = [];
if AddInfo.isDataMixed
  nSlices = length(DATA);
  fit_vals = [];
else
  nSlices = size(DATA{1}.FTDATA,3);
  
  % Try to get fit values from PROCPAR
  if isfield(DATA{1},'PROCPAR') && ~isempty(DATA{1}.PROCPAR) && ...
        isfield(DATA{1}.PROCPAR,'ne')
    
    if DATA{1}.PROCPAR.ne==8
        fit_vals=cumsum(ones(1,6)*DATA{1}.PROCPAR.te).*1000;
        nSlices = 6;
        DATA{1}.FTDATA = DATA{1}.FTDATA(:,:,1:6);
    else
        fit_vals=DATA{1}.PROCPAR.te.*1000;
    end 
  end
end


resp = aedes_inputdlg('Type TE values','Input dialog',...
  {mat2str(fit_vals)});
if isempty(resp)
  return
else
  resp=resp{1};
  fit_vals = str2num(resp);
end

% Prompt for file name
[fname,fpath,findex]=uiputfile({'*.t2;*.T2;*.s2;*.S2',...
                    'T2-Files (*.t2, *.s2)';...
                    '*.*','All Files (*.*)'},...
                               'Save T2-file',[DATA{1}.HDR.fpath, ...
                    't2_map']);
if isequal(fname,0) || isequal(fpath,0)
  return
end

% Calculate the map
[fp,fn,fe]=fileparts([fpath,fname]);
try
  aedes_fitmaps(DATA,'t2',fit_vals,'FileName',[fp,filesep,fn]);
catch
  errordlg({'Could not calculate T2 maps. The following error was returned',...
           '',lasterr},'modal')
end

