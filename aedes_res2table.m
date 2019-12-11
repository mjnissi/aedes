function ResTable = aedes_res2table(Res,varargin)
% AEDES_RES2TABLE - 
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


% Defaults
SortByRoi = false; % Sort by ROIs or FileNames.
DecSep = ',';
NumDec = 3;
ResTable = {};
ResFileName = '';
dirs = [0 1 2 3 4]; % By default show all directions (1,2,3,4) and total (0)

% Check Res-structure
if ~isstruct(Res) && ~isfield(Res,'Stat') && ...
      ~isfield(Res,'FileInfo') && ~isfield(Res,'DateTime')
  error('Invalid Res-structure')
end


% Parse varargin
for ii=1:2:length(varargin)
  switch lower(varargin{ii})
   case 'sortbyroi'
    SortByRoi = varargin{ii+1};
   case 'decsep'
    DecSep = varargin{ii+1};
   case 'numdec'
    NumDec = varargin{ii+1};
   case 'resfilename'
    ResFileName = varargin{ii+1};
   case {'directions','dirs'}
    tmp = varargin{ii+1};
    dirs = find(ismember('TXYZV',tmp))-1;
   otherwise
    error('Unknown parameter %s',varargin{ii})
  end
end

FileInfo = Res.FileInfo;
DateTime = Res.DateTime;

nRois = length(Res.Stat);

% Ensure backward compatibility
if ~isfield(Res.Stat(1),'Sum')
  for ii=1:length(Res.Stat)
    Res.Stat(ii).Sum = NaN(1,length(Res.Stat(ii).Mean));
    if isfield(Res.Stat(1),'XD')
      Res.Stat(ii).XD.Sum = NaN(1,length(Res.Stat(ii).XD.Mean));
      Res.Stat(ii).YD.Sum = NaN(1,length(Res.Stat(ii).YD.Mean));
      Res.Stat(ii).ZD.Sum = NaN(1,length(Res.Stat(ii).ZD.Mean));
      Res.Stat(ii).VD.Sum = NaN(1,length(Res.Stat(ii).VD.Mean));
    end
  end
end

% Allocate space for ResTable
if Res.Stat(1).isMixed
  if SortByRoi
    ResTable=cell(1,7);
  else
    ResTable=cell(1,8);
  end
else
  %ResTable=cell(7+nRois+5*3+nRois*3*length(Res.Stat(1).Mean),6);
  ResTable=cell(1,7);
end

%size(ResTable)
% Write date, time and Res filename
ResTable{1,1} = ResFileName; % Res-file name
ResTable{2,1} = DateTime; % 


% concatenate file and path names
fnames={};
for ii=1:length(Res.FileInfo.DataFileName)
  fnames{ii} = [Res.FileInfo.DataPathName{ii},Res.FileInfo.DataFileName{ii}];
end
% $$$ fnames=cellstr([char(Res.FileInfo.DataPathName) ...
% $$$                char(Res.FileInfo.DataFileName)]);

%% Parse Res structure
start_ind = 3;

% Parse mixed type results
if Res.Stat(1).isMixed
  for ii=1:nRois
    if SortByRoi
      start_ind = start_ind+2;
      mtrx=[Res.Stat(ii).Mean(:) Res.Stat(ii).Std(:) Res.Stat(ii).Sum(:) ...
            Res.Stat(ii).Min(:) Res.Stat(ii).Max(:) ...
            Res.Stat(ii).PixelCount(:)];
      ResTable{start_ind,1}=['ROI: ',Res.Stat(ii).Label];
      ResTable(start_ind+1:start_ind+2,:) = {'Filename','Mean','STD','Sum','Min','Max',...
                          'Pixel count';'------','------','------','------','------',...
                          '------','------'};
      ResTable(start_ind+3:start_ind+3+length(Res.Stat(ii).Mean)-1,1) = fnames;
      ResTable(start_ind+3:start_ind+3+length(Res.Stat(ii).Mean)-1,2:end-1) = ...
          strrep(aedes_cellsprintf(['%.' num2str(NumDec) 'f'],mtrx(:,1:end-1)),'.', ...
                 DecSep);
      ResTable(start_ind+3:start_ind+3+length(Res.Stat(ii).Mean)-1,end) = ...
          strrep(aedes_cellsprintf(['%.0f'],mtrx(:,end)),'.', ...
                 DecSep);
      start_ind = start_ind+3+size(mtrx,1);
    else
      mtrx=[Res.Stat(ii).Mean(:) Res.Stat(ii).Std(:) Res.Stat(ii).Sum(:) Res.Stat(ii).Min(:) Res.Stat(ii).Max(:) ...
            Res.Stat(ii).PixelCount(:)];
      ResTable(start_ind+1:start_ind+2,:) = {'Filename','ROI','Mean','STD','Sum','Min','Max',...
                          'Pixel count';'------','------','------','------','------',...
                          '------','------','------'};
      ResTable(start_ind+2+ii:nRois:start_ind+3+nRois*length(Res.Stat(ii).Mean)-1,1) ...
          = fnames;
      ResTable(start_ind+2+ii:nRois:start_ind+3+nRois*length(Res.Stat(ii).Mean)-1,2) ...
          = {Res.Stat(ii).Label};
      ResTable(start_ind+2+ii:nRois:start_ind+3+nRois*length(Res.Stat(ii).Mean)-1,3:end-1) ...
          = strrep(aedes_cellsprintf(['%.' num2str(NumDec) 'f'],mtrx(:,1:end-1)),'.', ...
                   DecSep);
      ResTable(start_ind+2+ii:nRois:start_ind+3+nRois*length(Res.Stat(ii).Mean)-1,end) = ...
          strrep(aedes_cellsprintf(['%.0f'],mtrx(:,end)),'.', ...
                 DecSep);
    end
  end
else % Parse normal type results

    % Construct totals
    if any(dirs==0)
      ResTable{5,1} = 'TOTAL';
      ResTable(6:7,:) = {'ROI','Mean','STD','Sum','Min','Max','Pixel count';...
        '------','------','------','------',...
        '------','------','------'};
      RowInd = 8;
      for kk=1:nRois
        totals = [Res.Stat(kk).Mean,Res.Stat(kk).Std,Res.Stat(kk).Sum,Res.Stat(kk).Min,...
          Res.Stat(kk).Max];
        ResTable{RowInd,1} = Res.Stat(kk).Label;
        ResTable(RowInd,2:6) = strrep(aedes_cellsprintf(['%.' num2str(NumDec) 'f'],...
          totals),'.',DecSep);
        ResTable(RowInd,7) =  strrep(aedes_cellsprintf('%.0f',Res.Stat(kk).PixelCount),'.',...
          DecSep);
        RowInd=RowInd+1;
      end
    else
      RowInd = 3;
    end
    
    dir_str = {'X-Direction','XD';...
               'Y-Direction','YD';...
               'Z-Direction','ZD';...
               'V-Direction','VD'};
    
    
    for kk=dirs
      if kk==0
        continue
      end
      
      % Construct X, Y, and Z directions
      RowInd = RowInd+2;
      D = dir_str{kk,2};
      ResTable{RowInd,1} = dir_str{kk,1};
      ResTable(RowInd+1:RowInd+2,:) = {'ROI','Mean','STD','Sum','Min','Max','Pixel count';...
                          '------','------','------','------',...
                          '------','------','------'};
      RowInd=RowInd+3;
      for ii=1:nRois
        totals = [Res.Stat(ii).(D).Mean(:),Res.Stat(ii).(D).Std(:),Res.Stat(ii).(D).Sum(:),...
                  Res.Stat(ii).(D).Min(:),...
                  Res.Stat(ii).(D).Max(:)];
        ResTable(RowInd:RowInd+length(totals(:,1))-1,1) = ...
            cellstr(repmat(Res.Stat(ii).Label,length(totals(:,1)),1));
        ResTable(RowInd:RowInd+length(totals(:,1))-1,2:end-1) = ...
            strrep(aedes_cellsprintf(['%.' num2str(NumDec) 'f'],...
                               totals),'.',DecSep);
        ResTable(RowInd:RowInd+length(totals(:,1))-1,end) = ...
            strrep(aedes_cellsprintf('%.0f',Res.Stat(ii).(D).PixelCount),'.',...
                   DecSep);
        RowInd = RowInd+length(totals(:,1));
      end
    end

end

vers = aedes_getmatlabversion;
if vers>=7.06
  % Convert empty matrices to empty strings
  ResTable(cellfun('isempty',ResTable))={''};
end
