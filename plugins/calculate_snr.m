function calculate_snr(DATA,ROI,AddInfo)
% calculate_snr - Calculate Signal-to-noise ratio (SNR) from two ROIs
%   (this is an Aedes plugin)
%
% This plugin looks for two ROIs labeled "signal" and "noise", which it
% will then use to calculate SNR. SNR value is printed in the workspace and
% show in a messagebox.

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


% -------------------------------------------------------------------------
% First check if two ROIs labeled "signal" and "noise" are defined. If not,
% display an error and return.
% -------------------------------------------------------------------------

% Check if ROI-structure is empty, i.e. no ROIs have been defined. Display
% an error ROIs have not been defined.
if isempty(ROI)
  errordlg(['No ROIs defined. Please draw ROIs for "signal" and "noise" ',...
    'and label them correspondingly.'],'SNR Plugin Error');
  return
end

% Check that there are at least two ROIs defined with the labels "signal"
% and "noise"
if length(ROI)<2
    errordlg(['Too few ROIs found. Please draw ROIs for "signal" and "noise" ',...
    'and label them correspondingly.'],'SNR Plugin Error');
  return
end

% Get indices for "signal" and "noise" ROIs
SignalRoiInd = find(strcmpi('signal',{ROI(:).label}));
NoiseRoiInd = find(strcmpi('noise',{ROI(:).label}));

% If either signal or noise ROI is not found, display an error and return
if isempty(SignalRoiInd) || isempty(NoiseRoiInd)
  errordlg(['Please draw ROIs for "signal" and "noise" ',...
    'and label them correspondingly.'],'SNR Plugin Error');
  return
end

% Display general info
fprintf(1,'\n') % empty line
disp('==========================')
disp('= Signal-to-noise ratios =')
disp('==========================')

% We almost always have to handle stacked and non-stacked data differently.
% Stacked or mixed data means that multiple slices have been opened in
% Aedes using the "Open multiple files" dialog. Non-stacked or non-mixed
% data set is e.g. a regular fid-file, which contains one or more slices.
if AddInfo.isDataMixed
  fprintf(1,'------------------------------\n')
  fprintf(1,'Slice nbr\tSNR\tFile Name\n')
  fprintf(1,'------------------------------\n')
  
  % Calculate SNR for every slice
  snr = zeros(length(DATA),1);
  for ii=1:length(DATA)
    signal = DATA{ii}.FTDATA(ROI(SignalRoiInd).voxels{ii});
    noise = DATA{ii}.FTDATA(ROI(NoiseRoiInd).voxels{ii});
    if isempty(signal) || isempty(noise)
      snr(ii) = NaN;
    else
      mean_signal = mean(signal);
      std_noise = std(noise);
      snr(ii) = mean_signal/std_noise;
    end
    fprintf(1,'Slice %d:\t%f\t%s%s\n',ii,snr(ii),...
      DATA{ii}.HDR.fpath,DATA{ii}.HDR.fname)
  end
else
  % Calculate SNR for every slice in the slice (X) direction. If there are
  % more than one volumes, use the current volume.
  fprintf(1,'Filename:\n%s%s\n',DATA{1}.HDR.fpath,DATA{1}.HDR.fname)
  fprintf(1,'------------------------------\n')
  fprintf(1,'Slice nbr\tSNR\n')
  fprintf(1,'------------------------------\n')
  
  snr = zeros(size(DATA{1}.FTDATA,3),1);
  for ii=1:size(DATA{1}.FTDATA,3)
    vol=AddInfo.CurrentVol;
    data = DATA{1}.FTDATA(:,:,ii,vol);
    signal=data(ROI(SignalRoiInd).voxels{1}(:,:,ii,vol));
    noise=data(ROI(NoiseRoiInd).voxels{1}(:,:,ii,vol));
    if isempty(signal) || isempty(noise)
      snr(ii) = NaN;
    else
      mean_signal = mean(signal);
      std_noise = std(noise);
      snr(ii) = mean_signal/std_noise;
    end
    fprintf(1,'Slice %d:\t%f\n',ii,snr(ii))
  end
end



