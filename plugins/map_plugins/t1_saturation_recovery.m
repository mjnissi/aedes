function t1_saturation_recovery(DATA,ROI,AddInfo)
% This Aedes plugin calculates T1 saturation recovery maps

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
else
  nSlices = size(DATA{1}.FTDATA,3);
end


resp = aedes_inputdlg('Type fit values','Input dialog',...
  {mat2str(fit_vals)});
if isempty(resp)
  return
else
  resp=resp{1};
  fit_vals = str2num(resp);
end

% Prompt for file name
[fname,fpath,findex]=uiputfile({'*.t1;*.T1;*.s1;*.S1',...
                    'T1-Files (*.t1, *.s1)';...
                    '*.*','All Files (*.*)'},...
                               'Save T1-file',[DATA{1}.HDR.fpath, ...
                    't1_sr_map']);
if isequal(fname,0) || isequal(fpath,0)
  return
end

% Calculate the map
[fp,fn,fe]=fileparts([fpath,fname]);
try
  aedes_fitmaps(DATA,'T1_SR',fit_vals,'FileName',[fp,filesep,fn]);
catch
  errordlg({'Could not calculate T1 SR maps. The following error was returned',...
           '',lasterr},'modal')
end
