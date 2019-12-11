function fat_analysis(DATA,ROI,AddInfo)
% FAT_ANALYSIS - Calculate fat percentage (Aedes plugin)
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


%% Check if the data to be analyzed is valid
if AddInfo.isDataMixed
  errordlg('Sorry. This plugin does not work with mixed data.','Error',...
    'modal');
  return
end

%% Check if the size of the data is correct
data_sz = size(DATA{1}.FTDATA);
if length(data_sz)<4 || data_sz(4)~=2
  errordlg(['The data should be a 4D-matrix containing fat and water ',...
    'volumes in the 4th dimesion.'],'Error',...
    'modal');
  return
end


%% Warn if data does not appear to be correct but give a choice to continue
%% anyway...
if ~isfield(DATA{1},'PROCPAR') || ...
    ~isfield(DATA{1}.PROCPAR,'seqfil') || ...
    ~strcmp(DATA{1}.PROCPAR.seqfil,'ge3d_csi2')
    
  resp=questdlg(['Data does not seem to be measured with a valid sequence. ',...
    'Do you want to continue anyway?'],...
    'Invalid sequence detected',...
    'Yes','No','No');
  if isempty(resp) || strcmpi(resp,'no')
    % Cancelled
    return
  end
end

%% Prompt for the slices that correspond to right anatomical area and have
%% tolerable level of noise/artefacts
canceled=false;
invalid_input = true;
while invalid_input
  resp = aedes_inputdlg(['Input the appropriate slice range'],...
    'Input slice range','1:end');
  if isempty(resp)
    % Cancelled
    canceled=true;
    invalid_input = false;
    continue
  end
  resp = strrep(resp{1},'end',num2str(data_sz(3)));
  z_range = str2num(resp);
  if ~isempty(z_range) && isnumeric(z_range) && isreal(z_range)
    invalid_input = false;
  else
    h=errordlg('Invalid range input','Invalid range input','modal');
    uiwait(h);
  end
end
if canceled
  return
end

%% Get FOV
if ~isfield(DATA{1},'PROCPAR') || ~isfield(DATA{1}.PROCPAR,'lro')
  %% Prompt for FOV
  resp = aedes_inputdlg('Input FOV in millimeters',...
    'Input FOV','[40 40 80]');
  if isempty(resp)
    % Canceled
    return
  end
  FOV = str2num(resp{1});
  if isempty(FOV) || length(FOV)~=3
    errordlg('Invalid FOV! Aborting...','Invalid FOV','modal');
    return
  end
else
  % FOV in millimeters
  FOV = [DATA{1}.PROCPAR.lro,DATA{1}.PROCPAR.lpe,...
    DATA{1}.PROCPAR.lpe2]*10;
end

% A brute force -method for correct DC artifacts, may not always work
data_water = DATA{1}.FTDATA(:,:,:,2);
data_fat = DATA{1}.FTDATA(:,:,:,1);
[mx,DC_Ind_fat]=max(max(max(data_fat,[],2),[],3));
[mx,DC_Ind_water]=max(max(max(data_water,[],2),[],3));

% Check that the found noise peak is near the center
if DC_Ind_fat < data_sz(1)/2-5 || DC_Ind_fat > data_sz(1)/2+5
  warning('Largest noise peak not near the center line in fat data! Check images and results.');
  DC_Ind_fat=data_sz(1)/2+1;
end
if DC_Ind_water < data_sz(1)/2-5 || DC_Ind_water > data_sz(1)/2+5
  warning('Largest noise peak not near the center line in water data! Check images and results.');
  DC_Ind_water=data_sz(1)/2+1;
end
data_fat(DC_Ind_fat,:,:) = mean(data_fat([DC_Ind_fat-1 DC_Ind_fat+1],:,:));
data_water(DC_Ind_water,:,:) = mean(data_water([DC_Ind_water-1 DC_Ind_water+1],:,:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FAT analysis -----------------------------------------


%% Set the default noise ROI
% noise_fat = 2*squeeze(mean(mean(DATA{1}.FTDATA(data_sz(1)-10:data_sz(1),1:data_sz(2),:,1),1),2));
% noise_water = 2*squeeze(mean(mean(DATA{1}.FTDATA(data_sz(1)-10:data_sz(1),1:data_sz(2),:,2),1),2));
noise_fat = 2*squeeze(mean(mean(data_fat(1:30,1:40,:),1),2));
noise_water = 2*squeeze(mean(mean(data_water(1:30,1:40,:),1),2));

mean_fat = squeeze(mean(mean(data_fat(:,:,:),1),2));
mean_water = squeeze(mean(mean(data_water(:,:,:),1),2));

fat_th = noise_fat+mean_fat;
water_th = noise_water+mean_water;
fat_th = fat_th(:);
water_th = water_th(:);

% Construct matrices for threshold comparison
fat_ind = repmat(permute(fat_th(z_range),[3 2 1]),...
  [data_sz(1) data_sz(2) 1]);
water_ind = repmat(permute(water_th(z_range),[3 2 1]),...
  [data_sz(1) data_sz(2) 1]);

% Calculate number of water and fat voxels
nVoxFat = numel(find(data_fat(:,:,z_range)>fat_ind));
nVoxWater = numel(find(data_water(:,:,z_range)>water_ind));

% Calculate fat and water concentrations
density_fat = 0.9; %g/ml;
voxel_weight_fat = (FOV(1)/data_sz(1)*FOV(2)/data_sz(2)*FOV(3)/data_sz(3))*density_fat/1000;
fat_weight = nVoxFat*voxel_weight_fat;

density_water = 1; %g/ml;
voxel_weight_water = (FOV(1)/data_sz(1)*FOV(2)/data_sz(2)*FOV(3)/data_sz(3))*density_water/1000;
water_weight = nVoxWater*voxel_weight_water;

% Fat-%
fat_percent = fat_weight/(fat_weight+water_weight)*100;

% Print results to the command window
fprintf(1,'\n####################################################\n');
fprintf(1,'FAT ANALYSIS RESULTS\n\n');
fprintf(1,'File: %s\n',[DATA{1}.HDR.fpath,DATA{1}.HDR.fname]);
fprintf(1,'FOV: %dx%dx%d mm\n',FOV);
fprintf(1,'Number of fat voxels:\t%d\n',nVoxFat);
fprintf(1,'Number of H2O voxels:\t%d\n\n',nVoxWater);
fprintf(1,'==========================================\n');
fprintf(1,'%s\n',fliplr(sprintf('%15s%15s%15s',...
  fliplr('Fat prc (%)'),fliplr('H2O (mg/ml)'),fliplr('Fat (mg/ml)'))));
fprintf(1,'------------------------------------------\n')
fprintf(1,'%s\n',fliplr(sprintf('%15s%15s%15s',...
  fliplr(sprintf('%.4f',fat_percent)),...
  fliplr(sprintf('%.4f',water_weight)),fliplr(sprintf('%.4f',fat_weight)))));
fprintf(1,'==========================================\n');
fprintf(1,'####################################################\n');


% Show results in uitable
fig_h = 150;
fig_w = 350;
fig_location = aedes_dialoglocation([fig_w,fig_h]);
fig_pos = [fig_location(1) fig_location(2) fig_w fig_h];
fh = figure('units','pixels',...
  'position',fig_pos,...
  'Name','Fat Analysis Results',...
  'IntegerHandle','off',...
  'Numbertitle','off',...
  'MenuBar','none');

matlab_ver = version;
if str2num(matlab_ver(1:3))<=7.5
  h_uitable = uitable('parent',fh,...
    'position',[5 5 340 55],...
    'Data',[fat_weight water_weight fat_percent],...
    'ColumnNames',{'Fat (mg/ml)','H2O (mg/ml)','Fat prc (%)'},...
    'ColumnWidth',100);
else
  h_uitable = uitable('parent',fh,...
    'position',[5 5 340 55],...
    'Data',[fat_weight water_weight fat_percent],...
    'ColumnName',{'Fat (mg/ml)','H2O (mg/ml)','Fat prc (%)'},...
    'ColumnWidth',{100 100 100});
end

pos = get(h_uitable,'position');
h_text = uicontrol('parent',fh,...
  'style','text',...
  'position',[pos(1) pos(2)+pos(4)+2 pos(3) fig_h-pos(2)-pos(4)-4],...
  'string',...
  {['File: ',DATA{1}.HDR.fpath,DATA{1}.HDR.fname],...
  sprintf('FOV: %dx%dx%d mm',FOV),...
  sprintf('Number of fat voxels: %d',nVoxFat),...
  sprintf('Number of H2O voxels: %d',nVoxWater)},...
  'horizontalalign','left');

if ispc
  set(h_text,'fontsize',8);
else
  set(h_text,'fontsize',8);
end


% Prompt for saving image
resp = questdlg(['Do you also want to save images of the fat/water distributions results?'],...
  'Save images?','Yes','No','Yes');
if isempty(resp) || strcmpi(resp,'No')
  % Cancelled
  return
end

% Prompt for output directory and file name
if ~isfield(DATA{1}.HDR,'fpath') || isempty(DATA{1}.HDR.fpath)
  default_name = ['fat_analysis_',datestr(now,30)];
else
  if strcmpi(DATA{1}.HDR.fpath(end),filesep)
    [fp,fn,fext]=fileparts(DATA{1}.HDR.fpath(1:end-1));
  else
    [fp,fn,fext]=fileparts(DATA{1}.HDR.fpath);
  end
  default_name = [fn,'_fat_analysis'];
end

try
  default_path = getpref('FatAnalysis','ImagePath');
catch
  default_path = [pwd,filesep];
end

[fname,fpath,findex] = uiputfile(...
  {'*.png','PNG Image Files (*.png)';...
  '*.tif','TIF Image Files (*.tif)'},...
  'Save Images to...',...
  [default_path,default_name]);
if isequal(fname,0) || isequal(fpath,0)
  % Cancelled
  return
end
setpref('FatAnalysis','ImagePath',fpath);

% Get sw from procpar or prompt for it if it is missing
if ~isfield(DATA{1},'PROCPAR') || ~isfield(DATA{1}.PROCPAR,'sw')
  resp = aedes_inputdlg('Please give sw value',...
    'Input sw value','');
  if isempty(resp)
    % Cancelled
    return
  end
  sw = str2num(resp{1});
  if isempty(sw) || ~isreal(sw)
    errordlg('Invalid value for sw! Aborting...','Invalid sw value','modal');
    return
  end
else
  sw = DATA{1}.PROCPAR.sw;
end


% A brute force -method for correct DC artifacts, may not always work
% data_water = DATA{1}.FTDATA(:,:,:,2);
% data_fat_shifted = DATA{1}.FTDATA(:,:,:,1);
% [mx,DC_Ind_fat]=max(max(DATA{1}.FTDATA(:,:,size(DATA{1}.FTDATA,3)/2,1),[],2));
% [mx,DC_Ind_water]=max(max(DATA{1}.FTDATA(:,:,size(DATA{1}.FTDATA,3)/2,2),[],2));
% data_fat_shifted(DC_Ind_fat,:,:) = mean(data_fat_shifted([DC_Ind_fat-1 DC_Ind_fat+1],:,:));
% data_water(DC_Ind_water,:,:) = mean(data_water([DC_Ind_water-1 DC_Ind_water+1],:,:));

%% Calculate Fat shift ------------------

% Show waitbar
[h,txt]=aedes_calc_wait({'Correcting fat shift.','This may take some time...'});
fat_shift = 680/sw;

data_fat = interp1(1:size(data_fat,1),...
  data_fat,...
  (1:size(data_fat,1))+fat_shift,'spline',0);

% If there is a huge DC artifact in the data the spile interpolation tends
% to go crazy near the artifacts and produce huge negative values...
data_fat(data_fat<0) = 0;

% Close waitbar
close(h)
% ---------------------------------------

% Try to set the intensity levels to reasonable values
%tmp_water=data_water;
tmp_water = data_water((size(data_water,1)/2-30):(size(data_water,1)/2+30),...
  (size(data_water,2)/2-30):(size(data_water,2)/2+30),:);
tmp_fat = data_fat((size(data_fat,1)/2-30):(size(data_fat,1)/2+30),...
  (size(data_fat,2)/2-30):(size(data_fat,2)/2+30),:);
data_fat = data_fat./max(tmp_fat(:));
data_water = data_water./max(tmp_water(:));
data_fat(data_fat>1) = 1;
data_water(data_water>1) = 1;

% Construct an RGB image of water and fat content
RGB_im = cat(4,data_fat,...
  zeros(size(data_fat)),...
  data_water);

% Write png or tif images
nImages = size(RGB_im,3);
[fp,fn,fext]=fileparts(fname);
if ~isempty(findex) && findex==2
  imformat = 'TIF';
  f_ext = '.tif';
else
  imformat = 'PNG';
  f_ext = '.png';
end
  
h = aedes_wbar(0,'');
for ii=1:nImages
  filename = sprintf([fn,'_%03d',f_ext],ii);
  aedes_wbar(ii/nImages,h,{sprintf('Writing image %d/%d to',ii,nImages),...
    [fpath,filename]});
  imwrite(squeeze(RGB_im(:,:,ii,:)),[fpath,filename],imformat);
end
close(h)
