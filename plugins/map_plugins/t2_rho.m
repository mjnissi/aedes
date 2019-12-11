function t2_rho(DATA,ROI,AddInfo)
% This Aedes plugin calculates T1 rho map

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
  
  % Try to get fit values from PROCPAR.pi
%   if isfield(DATA{1},'PROCPAR') && ~isempty(DATA{1}.PROCPAR) && ...
%         isfield(DATA{1}.PROCPAR,'pi') && length(DATA{1}.PROCPAR.pi)>1
%     fit_vals = DATA{1}.PROCPAR.pi./1000;
%   end
end

resp = aedes_inputdlg('Type spin-lock values');
if isempty(resp)
  return
else
  resp=resp{1};
  fit_vals = str2num(resp);
end


% Prompt for file name
[fname,fpath,findex]=uiputfile({'*.t2r;*.T2R;*.s2r;*.S2R',...
                    'T2R-Files (*.t2r, *.s2r)';...
                    '*.*','All Files (*.*)'},...
                               'Save T2R-file',[DATA{1}.HDR.fpath, ...
                    't2rho_map']);
if isequal(fname,0) || isequal(fpath,0)
  return
end

% Calculate the map
[fp,fn,fe]=fileparts([fpath,fname]);
try
  aedes_fitmaps(DATA,'t2r',fit_vals,'FileName',[fp,filesep,fn]);
catch
  errordlg({'Could not calculate T2 rho maps. The following error was returned',...
           '',lasterr},'modal')
end

