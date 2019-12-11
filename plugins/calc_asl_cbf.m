function calc_asl_cbf(DATA,ROI,AddInfo)
% CALC_ASL_CBF - Calculate CBF maps from ASL data using T1-map data if
% possible (Aedes plugin)
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
lambda = 0.9;
t1map = 1.2;

if AddInfo.isDataMixed
  % Check that the number of images is OK
  if rem(length(DATA),2)~=0 || length(DATA)==1
	h=errordlg('The number of input images has to be divisible by 2!',...
	  'Error','modal');
	return
  end
  controlIm = zeros([size(DATA{1}.FTDATA) length(DATA)/2]);
  labelIm = zeros([size(DATA{1}.FTDATA) length(DATA)/2]);
  for ii=1:2:size(controlIm,3)
	controlIm(:,:,ii) = DATA{ii}.FTDATA;
	labelIm(:,:,ii) = DATA{ii+1}.FTDATA;
  end
  for ii=1:length(DATA)
	ParentFileName{ii} = fullfile(DATA{ii}.HDR.fpath,DATA{ii}.HDR.fname);
  end
else
	
	if strcmpi(DATA{1}.DataFormat,'bruker_reco')
		% Assume that Bruker data is 4D
		if rem(size(DATA{1}.FTDATA,4),2)~=0 || size(DATA{1}.FTDATA,4)==1
			h=errordlg('The number of input images has to be divisible by 2!',...
				'Error','modal');
			return
		end
		
		controlIm = squeeze(DATA{1}.FTDATA(:,:,2,:));
		labelIm = squeeze(DATA{1}.FTDATA(:,:,1,:));
		ParentFileName = fullfile(DATA{1}.HDR.fpath,DATA{1}.HDR.fname);
	else
		% Check that the number of images is OK
		if rem(size(DATA{1}.FTDATA,3),2)~=0 || size(DATA{1}.FTDATA,3)==1
			h=errordlg('The number of input images has to be divisible by 2!',...
				'Error','modal');
			return
		end
		if isfield(DATA{1}.PROCPAR,'ltype') && ...
				length(DATA{1}.PROCPAR.ltype)==size(DATA{1}.FTDATA,3)
			control_ind = find(strcmpi(DATA{1}.PROCPAR.ltype,'c'));
			label_ind = find(strcmpi(DATA{1}.PROCPAR.ltype,'l'));
			controlIm = DATA{1}.FTDATA(:,:,control_ind);
			labelIm = DATA{1}.FTDATA(:,:,label_ind);
			ParentFileName = fullfile(DATA{1}.HDR.fpath,DATA{1}.HDR.fname);
		else
			controlIm = DATA{1}.FTDATA(:,:,1:2:size(DATA{1}.FTDATA,3));
			labelIm = DATA{1}.FTDATA(:,:,2:2:size(DATA{1}.FTDATA,3));
			ParentFileName = fullfile(DATA{1}.HDR.fpath,DATA{1}.HDR.fname);
		end
	end
end

% Ask for lambda
resp = aedes_inputdlg('Lambda coefficient?','Lambda coefficient?',...
  num2str(lambda));
if isempty(resp)
  % Canceled
  return
end
lambda = str2num(resp{1});

% Ask for T1-map
default_path = DATA{1}.HDR.fpath;
if isempty(default_path)
  default_path = [pwd,filesep];
end

resp = questdlg(['Select if you want to use a constant value for T1 ',...
  'or an existing T1-map.'],...
  'Use constant T1 value or T1-map?',...
  'Use T1-map','Use constant T1','Cancel',...
  'Use T1-map');
if isempty(resp) || strcmpi(resp,'Cancel')
  % Canceled
  return
elseif strcmpi(resp,'Use constant T1')
  % Ask for T1 value
  resp = aedes_inputdlg('Input T1 value','Input T1 value',...
	num2str(t1map));
  if isempty(resp)
	% Canceled
	return
  end
  fname = '';
  fpath = '';
  t1map = str2num(resp{1});
  t1val = t1map;
else
  [fname,fpath,findex] = uigetfile({'*.t1','T1-Files (*.t1)';...
	'*.*','All Files (*.*)'},'Select T1-map to open',...
	default_path);
  if isequal(fname,0)
	% Canceled
	return
  end
  try
	tmp=load(fullfile(fpath,fname),'-mat');
	t1map = tmp.Data;
  catch
	h=errordlg({'Could not open T1-map from file',...
	  fullfile(fpath,fname)},'Error','modal');
	return
  end
  

  % Check that the T1-map size matches with ASL. If it doesn't, resize it to
  % match...
  t1map = t1map(:,:,1);
  if ~isequal(size(t1map),[size(controlIm,1),size(controlIm,2)])
	t1map = imresize(t1map,[size(controlIm,1),size(controlIm,2)]);
  end
  t1val = [];
end

% Ask if the user wants to calculate a single asl map or time series
resp = questdlg(['Calculate an average CBF map (recommended) ',...
	'or a CBF time series (i.e. one CBF map / image pair)?'],...
  'Calculate an average CBF map?',...
  'Average CBF','CBF map / image pair','Cancel',...
  'Average CBF');
if isempty(resp) || strcmpi(resp,'Cancel')
	% Canceled
	return
end
if strcmpi(resp,'Average CBF')
	CBFts = false;
else
	CBFts = true;
end

% Convert milliseconds to seconds
%t1map=t1map./1000;

if CBFts
	% Calculate differences and CBF maps
	asl_map = zeros(size(controlIm));
	cbf_map = zeros(size(controlIm));
	for ii=1:size(controlIm,3)
		asl_map(:,:,ii) = controlIm(:,:,ii)-labelIm(:,:,ii);
		cbf_map(:,:,ii) = (100.*60.*lambda.*asl_map(:,:,ii))./(t1map.*2.*controlIm(:,:,ii));
	end

else
	% Calculate difference fotr ASL
	mean_control = mean(controlIm,3);
	mean_label = mean(labelIm,3);
	asl_map = mean_control-mean_label;
	
	% Calculate CBF and save it into a MAT-file
	cbf_map = (100.*60.*lambda.*asl_map)./(t1map.*2.*mean_control);

end
cbf_map(isinf(cbf_map))=0;
cbf_map(isnan(cbf_map))=0;

% Ask where to save the resulting CBF map
[fn,fp,fi] = uiputfile({'*.mat','Matlab MAT-Files (*.mat)';...
	'*.*','All Files (*.*)'},'Save CBF-map as...',...
	[default_path,'asl_cbf_map.mat']);
if isequal(fname,0)
	% Canceled
	return
end

% Construct the map structure
Data = cbf_map;
Param.ParentFileName = ParentFileName;
Param.T1FileName = fullfile(fpath,fname);
Param.T1value = t1val;
Param.Lambda = lambda;

% Save the map
try
  save(fullfile(fp,fn),'Data','Param','-mat');
catch
  h=errordlg({'Could not save CBF-map to',...
	fullfile(fp,fn)},'Error','modal');
  return
end

% Open in a new Aedes window
aedes(fullfile(fp,fn))










