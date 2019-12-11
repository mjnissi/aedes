function map_struct = aedes_fitmaps(DATA,maptype,fit_vals,varargin)
% AEDES_FITMAPS - Calculate various parameter maps (T1, T2, T1rho, B1,
%               Perfusion, etc.) 
%   
%
% Synopsis: 
%       MapStruct = aedes_fitmaps(data,maptype,fit_vals,...
%                                'property1',value1,'property2',value2,...)
%
% Description:
%       Function fits various parameter maps defined by MAPTYPE from the
%       slice data DATA using the fit values FIT_VALS. The DATA can be a 3D
%       matrix or an Aedes DATA-structure. The MAPTYPE variable is a string
%       variable that defines the type of the fitted map, i.e. the function
%       used in the fitting. Valid values for MAPTYPE are the following:
%
%       'T1_IR'     <-> T1 map (inversion recovery)
%       'T1_SR'     <-> T1 map (saturation recovery)
%       'T1_3P'     <-> T1 map (3-parameter fit)
%       'T2'        <-> T2 map
%       'R2'        <-> R2 map
%       'T1r'       <-> T1 rho map
%       'T2r'       <-> T2 rho map
%       'ADC'       <-> Apparent diffusion coefficient map
%       'perfusion' <-> Perfusion map
%       'B1'        <-> B1 map
%
%       If you want to omit some slices from the map calculations, just
%       set the corresponding value in FIT_VALS to NaN.
%
%       Property-value pairs can be used to control additional options in
%       the fitting of maps (The { } denotes default value for the
%       corresponding property): 
%
%       Property:           Value:              Description:
%       ********            ********            ************
%       'FileName'          String              % Save maps to files
%
%       'SaveSpinDensities' ['on' | {'off'}]    % Save also the
%                                               % S0-parameter in a   
%                                               % separate file.
%                                               % This property is ignored
%                                               % if the FileName -property
%                                               % is not defined
%
%       'Wbar'              [{'on'} | 'off' ]   % Show/hide waitbar
%
%       'Linear'            [{'on'} | 'off']    % Perform linear or
%                                               % non-linear fit
%
%       'InitVal'           vector              % Initial values for
%                                               % non-linear fit
%
%       'MaxIter'           scalar              % Maximum number of
%                                               % iterations in non-linear
%                                               % fits (default=50)
%
%       The 'FileName' property can be used to write the map into a file
%       with the corresponding file extension (.t1, .t2, .t1r, etc.). The
%       files are normal Matlab MAT-files with a different file extension
%       and all the parameters used to calculate the map are also stored in
%       the file. Aedes can read these files normally. You can also load
%       these files into Matlab workspace by typing
%       maps=load('/path/to/my/mapfile','-mat'). If the 'FileName' property
%       does not include full path, the map-files a written into the same
%       directory as the data by default if possible. Otherwise the maps
%       are written into the current directory.
%
%       The 'Linear' property is ignored in the cases where only linear or
%       non-linear fit is available.
%
%       The 'InitVal' property contains the global initial values that are
%       used for all individual fittings. If the 'InitVal' property is
%       omitted, some "optimal" initial values are estimated using the data
%       and the fit values.
%
%       The function outputs a structure containing the following fields:
%       
%       MapStruct
%               |-> Map             :(The calculated map(s))
%               |-> S0              :(The "spin density")
%               |-> FitValues       :(Values used in the fitting)
%               |-> Type            :(Used map type, e.g. 'T2')
%               |-> Linear          :(1=linear fitting, 0=nonlinear fitting)
%               |-> MaxIter         :(Maximum number of iterations, non-linear only)
%               |-> ParentFileName  :(File(s) from which the map was calculated)
%
%
% Examples:
%
% See also:
%       AEDES

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

map_struct=[];
if nargin<3
  error('Too few input arguments!')
end

Dat.filename = '';
Dat.linear = true;
Dat.wbar = true;
Dat.MaxIter = 75; % Maximum number of iterations in nonlinear fits
Dat.initVal = []; % Empty means that defaults are used for initial values
ParentFileName = '';
Dat.UseComplexData = false;
Dat.SaveSpinDensities = false;
Dat.Mask = [];

% Parse additional input arguments
for ii=1:2:length(varargin)
  switch lower(varargin{ii})
    case 'filename'
      Dat.filename = varargin{ii+1};
    case 'linear'
      if ischar(varargin{ii+1})
        if strcmpi(varargin{ii+1},'off')
          Dat.linear = false;
        end
      elseif islogical(varargin{ii+1})
        Dat.linear = varargin{ii+1};
      else
        error('Invalid value for the LINEAR property!');
      end
    case {'wbar','waitbar'}
      if ischar(varargin{ii+1})
        if strcmpi(varargin{ii+1},'off')
          Dat.wbar = false;
        end
      elseif islogical(varargin{ii+1})
        Dat.wbar = varargin{ii+1};
      else
        error('Invalid value for the WBAR property!');
      end
    case {'initval','initialvalues'}
      Dat.initVal = varargin{ii+1};
    case 'maxiter'
      Dat.MaxIter = varargin{ii+1};
    case 'usecomplexdata'
      Dat.UseComplexData = varargin{ii+1};
    case 'savespindensities'
      if ischar(varargin{ii+1})
        if strcmpi(varargin{ii+1},'on')
          Dat.SaveSpinDensities = true;
        end
      elseif islogical(varargin{ii+1})
        Dat.SaveSpinDensities = varargin{ii+1};
      else
        error('Invalid value for the WBAR property!');
      end
    case 'mask'
      Dat.Mask = varargin{ii+1};
    otherwise
      error('Invalid property "%s"!',varargin{ii})
  end
end

if iscell(DATA) && length(DATA)==1
  DATA=DATA{1};
end

%% Check map type
if ~ischar(maptype)
  error(['Map type has to be one of the following strings: ',...
	'%s, %s, %s, %s, %s, %s, %s.'],...
	'T1','T2','T1r','T2r','perfusion','ADC','B1')
end

%% Convert data to 3D-matrix
if iscell(DATA) && length(DATA)>1
  sz = size(DATA{1}.FTDATA);
  sz(3) = length(DATA);
  datamtx = zeros(sz);
  ParentFileName = {};
  for ii=1:sz(3)
    datamtx(:,:,ii)=double(DATA{ii}.FTDATA);
	ParentFileName{ii} = fullfile(DATA{ii}.HDR.fpath,DATA{ii}.HDR.fname);
  end
elseif isstruct(DATA)
  if ndims(DATA.FTDATA)==4
    DATA.FTDATA = squeeze(DATA.FTDATA);
  end
  sz = size(DATA.FTDATA);
  datamtx = double(DATA.FTDATA);
  ParentFileName = fullfile(DATA.HDR.fpath,DATA.HDR.fname);
else
  if ndims(DATA)==4
    DATA = squeeze(DATA);
  end
  sz = size(DATA);
  datamtx = double(DATA);
end

%% Check the number of maps to be calculated
if rem(sz(3),length(fit_vals))~=0
  error('Number of slices doesn''t match with the length of fit values.');
  return
elseif length(fit_vals)<2
  % Display error if the number of fit values is not acceptable
  error('Cannot calculate maps with less than 2 fit values!')
else
  Dat.nMaps = sz(3)/length(fit_vals);
  
  %% NaN:s in fit_vals means that the corresponding slices will be omitted
  %% from the calculations
  indNans = find(isnan(fit_vals));
  if not(isempty(indNans))
    l=length(fit_vals);
    ind = repmat(indNans,1,Dat.nMaps)+reshape(repmat(l.*[0:Dat.nMaps-1],length(indNans),1),[],1)';
    datamtx(:,:,ind) = [];
    fit_vals(indNans)=[];
    sz=size(datamtx);
  end
  
end

%% Check that mask is of correct size
if ~isempty(Dat.Mask)
  if ~( size(Dat.Mask,1)==sz(1) && ...
      size(Dat.Mask,2)==sz(2) && ...
      size(Dat.Mask,3)==Dat.nMaps)
    error('Mask size does not correspond with data size!')
  end
else
  Dat.Mask = true(sz(1),sz(2),Dat.nMaps);
end

switch lower(maptype)
  %% Calculate ADC-maps ----------------------------
  case {'adc','df','diffusion'}
    
    % Fit ADC map
    map_out=l_ADC_map(datamtx,fit_vals,'',Dat);
    map_out.Type = 'ADC';
    map_out.Linear = true;
    fext{1} = 'df';
    fext{2} = 'sf';
    
    %% Calculate T1-maps ----------------------------
  case {'t1_ir','t1_sr','t1_3p','t1ir','t1sr','t13p'}
    
	if any(strcmpi(maptype,{'t1_ir','t1ir'}))
	  maptype = 'T1_IR';
	elseif any(strcmpi(maptype,{'t1_sr','t1sr'}))
	  maptype = 'T1_SR';
	else
	  maptype = 'T1_3P';
	end
	
    % Fit T1 maps
    map_out = l_T1_map(datamtx,fit_vals,maptype,Dat);
    map_out.Type = maptype;
    map_out.Linear = false;%Dat.linear;
	map_out.MaxIter = Dat.MaxIter;
    fext{1} = 't1';
    fext{2} = 's1';
    
    %% Calculate T1rho-maps ----------------------------
  case {'t1r','t1rho'}
    
	maptype = 'T1r';
	
    % Fit T1r maps
    map_out = l_T2_map(datamtx,fit_vals,maptype,Dat);
    map_out.Type = 'T1r';
    map_out.Linear = Dat.linear;
    fext{1} = 't1r';
    fext{2} = 's1r';
    
	%% Calculate T2rho-maps ----------------------------
  case {'t2r','t2rho'}
    
	maptype = 'T2r';
	
	% Fit T2r maps
	map_out = l_T2_map(datamtx,fit_vals,maptype,Dat);
    map_out.Type = 'T2r';
    map_out.Linear = Dat.linear;
    fext{1} = 't2r';
    fext{2} = 's2r';
	
	
    %% Calculate T2 maps ----------------------------
  case {'t2','r2'}
    
	if strcmpi(maptype,'t2')
	  % Fit T2 maps
	  map_out = l_T2_map(datamtx,fit_vals,maptype,Dat);
	  map_out.Type = 'T2';
	  map_out.Linear = true;
	  fext{1} = 't2';
	  fext{2} = 's2';
	elseif strcmpi(maptype,'r2')
	  % Fit R2 maps
	  map_out = l_T2__map(datamtx,fit_vals,maptype,Dat);
	  map_out.Type = 'R2';
	  map_out.Linear = Dat.linear;
	  map_out.Map = 1./map_out.Map;
	  map_out.Map(isinf(map_out.Map))=0;
	  map_out.Map(isnan(map_out.Map))=0;
	  fext{1} = 'r2';
	  fext{2} = 's2';
	end
    
    
    %% Calculate perfusion-maps ----------------------------
  case {'perfusion','perf'}
    
    
    
     %% Calculate B1-maps ----------------------------
  case 'b1'
    maptype = 'B1';
	
    % Fit T1r maps
	if isstruct(DATA) && isfield(DATA,'KSPACE') && ...
		~isempty(DATA.KSPACE) && Dat.UseComplexData
	  datamtx = double(DATA.KSPACE);
	  datamtx=fftshift(fftshift(fft(fft(datamtx,[],1),[],2),1),2);
	  fprintf(1,'Calculating B1-maps using complex data...\n')
	  map_out = l_B1_map(datamtx,fit_vals,maptype,Dat);
	else
	  map_out = l_B1_map(datamtx,fit_vals,maptype,Dat);
	end
    map_out.Type = maptype;
    map_out.Linear = false;%Dat.linear;
	map_out.MaxIter = Dat.MaxIter;
    fext{1} = 'b1';
    fext{2} = 'mat';
    
  case 'mt'
    
    map_out=l_MTmap(datamtx,fit_vals,Dat);
    S0=[];
      
  otherwise
    error('Unknown map type "%s"!',maptype)
end

% Add name of the parent file to the structure
map_out.ParentFileName = ParentFileName;

%% Write maps to files
if ~isempty(Dat.filename)
  if nargout==0
	clear map_struct
  end
  
  [fp,fn,fe]=fileparts(Dat.filename);
  if isempty(fp)
	if ~isempty(ParentFileName)
	  if iscell(ParentFileName)
		[p,n,e]=fileparts(ParentFileName{1});
		fpath = [p,filesep];
	  else
		[p,n,e]=fileparts(ParentFileName);
		fpath = [p,filesep];
	  end
	else
	  fpath = [pwd,filesep];
	end
  else
    fpath = [fp,filesep];
  end
  if isempty(fn)
    fname = 'mapdata';
  else
    fname = fn;
  end
  for ii=1:size(map_out.Map,3)
    Data = map_out.Map(:,:,ii);
	Param = map_out;
	Param = rmfield(Param,'Map');
    
    % Save map
    filename = sprintf('%s%s_%03d.%s',fpath,fname,ii,fext{1});
    save(filename,'Data','Param','-mat')
    
    % Save "spin densities"
	if isfield(map_out,'S0') && Dat.SaveSpinDensities
	  Data = map_out.S0(:,:,ii);
	  filename = sprintf('%s%s_%03d.%s',fpath,fname,ii,fext{2});
	  save(filename,'Data','Param','-mat')
	end
  end
else
  map_struct = map_out;
end



%%%%%%%%%%%%%%%%%%%%%
% Calculate T1-map
%%%%%%%%%%%%%%%%%%%%%
function map_out=l_T1_map(data,TI,maptype,Dat)
% Calculate T1 relaxation times
%
% T1 functions:
% Inversion recovery  :   S=abs(S0*(1-2*exp(-TI/T1)))
% Saturation recovery :	S=S0*(1-1*exp(-TI/T1))
% 3 Parameter fit     :	S=S0*(1-A*exp(-TI/T1))
%
% where:
% S = measured signal, S0 = "spin density"
% TI = inversion time, T1 = T1 relaxation time

%% Data size
sz=size(data);

%% Fit values
TI = TI(:);

%% Data slice index matrix
IndMtx = reshape(1:sz(3),[length(TI) Dat.nMaps])';

%% Allocate space for parameters
map = zeros([sz(1) sz(2) Dat.nMaps]);
S0 = zeros([sz(1) sz(2) Dat.nMaps]);
A = zeros([sz(1) sz(2) Dat.nMaps]);

% Fit options for fminsearch
options = optimset('Display','off',...
  'MaxIter',Dat.MaxIter);

estimateInitVal = false;
if isempty(Dat.initVal)
  estimateInitVal = true;
end

% - Inversion recovery 
if strcmpi(maptype,'t1_ir')
  t1_map_type = 1;
elseif strcmpi(maptype,'t1_sr')
  % Saturation recovery
  t1_map_type = 2;
else
  % 3 Parameter fit
  t1_map_type = 3;  
end

% Variables for waitbar
nFits = (sz(1)*sz(2)*Dat.nMaps);
counter = 1;
meanTI = mean(TI);

% Calculate maps
for ii=1:Dat.nMaps
  if Dat.wbar && ii==1
    wbh=aedes_wbar(0,sprintf('Processing map %d/%d',ii,Dat.nMaps));
  elseif Dat.wbar
    aedes_wbar(counter/nFits,wbh,sprintf('Processing map %d/%d',ii,Dat.nMaps));
  end
  
  for kk=1:sz(1)
    for tt=1:sz(2)
      if ~isempty(Dat.Mask) && ~Dat.Mask(kk,tt,ii)
        if Dat.wbar
          aedes_wbar(counter/nFits,wbh);
        end
        counter = counter+1;
        continue
      end
      
      % Data values
      data_val = squeeze(data(kk,tt,IndMtx(ii,:)));
      data_val = data_val(:);
      
      % Initial values for the fit
      if estimateInitVal
        if t1_map_type == 3
          init_val = [max(data_val); meanTI; 2];
        else
          init_val = [max(data_val); meanTI];
        end
      else
        init_val = Dat.initVal;
      end
      
      % Nelder-Mead simplex iteration
      if t1_map_type == 1
        fhandle = @(x) norm(data_val - abs(x(1)*(1-2*exp(-TI./x(2)))));
      elseif t1_map_type == 2
        fhandle = @(x) norm(data_val - abs(x(1)*(1-1*exp(-TI./x(2)))));
      else
        fhandle = @(x) norm(data_val - abs(x(1)*(1-x(3)*exp(-TI./x(2)))));
      end
      th = fminsearch(fhandle,init_val,options);
      
      S0(kk,tt,ii) = th(1);
      map(kk,tt,ii) = th(2);
      if t1_map_type == 3
        A(kk,tt,ii) = th(3);
      end
      
      if Dat.wbar
        aedes_wbar(counter/nFits,wbh);
      end
      counter = counter+1;
    end
  end
end
if Dat.wbar
  close(wbh)
end

map(map<0)=0;
map(isinf(map))=0;
map(isnan(map))=0;

map_out.Map = map;
map_out.S0 = S0;
map_out.Angle = A;
map_out.FitValues = TI;



%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate B1-map
%%%%%%%%%%%%%%%%%%%%%%%%%%
function map_out = l_B1_map(data,fit_vals,maptype,Dat)

%% Data size
sz=size(data);

%% Fit values
fit_vals = fit_vals(:);

%% Data slice index matrix
IndMtx = reshape(1:sz(3),[length(fit_vals) Dat.nMaps])';

%% Allocate space for parameters
map = zeros([sz(1) sz(2) Dat.nMaps]);
A = zeros([sz(1) sz(2) Dat.nMaps]);
phase = zeros([sz(1) sz(2) Dat.nMaps]);

% Fit options for fminsearch
options = optimset('Display','off',...
  'MaxIter',Dat.MaxIter);

estimateInitVal = false;
if isempty(Dat.initVal)
  estimateInitVal = true;
end


% Variables for waitbar
nFits = (sz(1)*sz(2)*Dat.nMaps);
counter = 1;


% Loop over maps
for ii=1:Dat.nMaps
  if Dat.wbar && ii==1
	wbh=aedes_wbar(0,sprintf('Calculating B1-map %d/%d',ii,Dat.nMaps));
  elseif Dat.wbar
	aedes_wbar(counter/nFits,wbh,sprintf('Calculating map %d/%d',ii,Dat.nMaps));
  end
  
  for tt=1:sz(1)
	for kk=1:sz(2)
	  % Data values and initial values for the iteration
	  if Dat.UseComplexData
		data_val = squeeze(data(kk,tt,IndMtx(ii,:)));
		data_val = real(data_val);
		data_val = data_val(:);
		
		if estimateInitVal
		  val1=length(find(diff(data_val<0)<0));
		  val2=length(find(diff(data_val<0)>0));
		  tot_t = fit_vals(end)-fit_vals(1);
		  omega = (((val1+val2)/2)/tot_t)*2*pi;
		  init_val = [2*max(data_val) omega 0.1];
		  init_val=init_val(:);
		else
		  init_val = Dat.initVal;
		end
		
	  else
		data_val = squeeze(data(kk,tt,IndMtx(ii,:)));
		data_val = data_val(:);
		
		if estimateInitVal
		  val1=length(find(diff(data_val<(max(data_val)*0.5))<0));
		  val2=length(find(diff(data_val<(max(data_val)*0.5))>0));
		  tot_t = fit_vals(end)-fit_vals(1);
		  omega = (((val1+val2)/4)/tot_t)*2*pi*0.9;
		  init_val = [max(data_val)-min(data_val) omega 0.1];
		  init_val=init_val(:);
		else
		  init_val = Dat.initVal;
		end
	  end
	  
	  % Function handle for fminsearch
	  if Dat.UseComplexData
		fhandle= @(x) norm(data_val - x(1)*cos(x(2)*fit_vals+x(3)));
	  else
		fhandle= @(x) norm(data_val - abs(x(1)*cos(x(2)*fit_vals+x(3))));	  
	  end
	  th = fminsearch(fhandle,init_val,options);
	  
	  map(kk,tt,ii) = abs(th(2))/(2*pi);
	  A(kk,tt,ii) = th(1);
	  phase(kk,tt,ii) = th(3);
	  
	  if Dat.wbar
		aedes_wbar(counter/nFits,wbh);
	  end
	  counter = counter+1;
	end
  end
end
if Dat.wbar
  close(wbh)
end

map(map<0)=0;
map(isinf(map))=0;
map(isnan(map))=0;

map_out.Map = map;
map_out.A = A;
map_out.phase = phase;
map_out.FitValues = fit_vals;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate T2-, T1rho- and T2rho-maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function map_out=l_T2_map(data,fit_vals,maptype,Dat)
% Calculate T1rho relaxation times
%
% T1rho functions:
%
% Linear:    log(S)=log(S0)-SL/T1rho
% Nonlinear: S=S0*exp(-SL/T1rho)
%
% where S=measured signal
%       S0 = "spin density"
%       SL = spinlock length
%       T1rho = T1rho relaxation time
%
% T2-functions:
%
% Linear:    log(S)=log(S0)-te/T2
% Nonlinear: S=S0*exp(-te/T2)
%
% where S=measured signal
%       S0 = "spin density"
%       te = echo time
%       T2 = T2 relaxation time

done=false;
map=[];
S0=[];


%% Make sure that fit values are in column vector
fit_vals=fit_vals(:);

%% Data size
sz=size(data);

%% Data slice index matrix
IndMtx = reshape(1:sz(3),[length(fit_vals) Dat.nMaps])';

map = zeros([sz(1) sz(2) Dat.nMaps]);
S0 = zeros([sz(1) sz(2) Dat.nMaps]);
A = zeros([sz(1) sz(2) Dat.nMaps]);

% Intensity values should not be less than zero
data(data<0)=0;

%% Use linearized form (fast)
if Dat.linear
  
  nFits = sz(2)*Dat.nMaps;
  counter = 1;
  
  H = [ones(size(fit_vals)) -fit_vals];
  for ii=1:Dat.nMaps
	if Dat.wbar && ii==1
	  wbh=aedes_wbar(0,sprintf('Processing map %d/%d',ii,Dat.nMaps));
	elseif Dat.wbar
	  aedes_wbar(counter/nFits,wbh,sprintf('Processing map %d/%d',ii,Dat.nMaps));
	end
	
    for kk=1:sz(2)
      tmp=squeeze(data(:,kk,IndMtx(ii,:))).';
      z=log(tmp);
      th=H\z;
      S0(:,kk,ii)=exp(th(1,:)');
      
      % Try to avoid "divide by zero" warnings
      
      map_tmp = th(2,:)';
      map_tmp(map_tmp==0)=eps;
      
      map(:,kk,ii)=1./map_tmp;
	  counter = counter+1;
	  if Dat.wbar
		aedes_wbar(counter/nFits,wbh);
	  end
    end
  end
  map(map<0)=0;
  map(isinf(map))=0;
  map(isnan(map))=0;
else
  %% Fit nonlinearly using Nelder-Mead Simplex (slow, but more accurate)
  estimateInitVal = false;
  if isempty(Dat.initVal)
	estimateInitVal = true;
  end
 
  % Fit options
  options = optimset('Display','off',...
	'MaxIter',Dat.MaxIter);

  nFits = sz(1)*sz(2)*Dat.nMaps;
  counter = 1;
  meanFV = mean(fit_vals);
  
  for ii=1:Dat.nMaps
	if Dat.wbar && ii==1
	  wbh=aedes_wbar(0,sprintf('Processing map %d/%d',ii,Dat.nMaps));
	elseif Dat.wbar
	  aedes_wbar(counter/nFits,wbh,sprintf('Processing map %d/%d',ii,Dat.nMaps));
	end
	
	for tt=1:sz(1)
	  for kk=1:sz(2)
		
		% Data for the fit
		z = squeeze(data(tt,kk,IndMtx(ii,:)));
		z=z(:);
		
		% Initial values for the fit
		if estimateInitVal
		  %init_val = [mean(z); mean(fit_vals); 1];
		  init_val = [max(z); meanFV];
		else
		  init_val = Dat.initVal;
		end
		
		
		%fhandle = @(x) sum((z - x(1)*exp(-fit_vals./x(2))+x(3)).^2);
		fhandle = @(x) norm(z - x(1)*exp(-fit_vals./x(2)));
		th = fminsearch(fhandle,init_val,options);
		 
		S0(tt,kk,ii) = th(1);
		map(tt,kk,ii) = th(2);
		%A(tt,kk,ii) = th(3);
	  
		if Dat.wbar
		  aedes_wbar(counter/nFits,wbh);
		end
		counter = counter+1;
	  end
	end
  end
 
end
if Dat.wbar
  close(wbh)
end


map_out.Map = map;
map_out.S0 = S0;
map_out.FitValues = fit_vals;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Perfusion-map
%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_Perfusion_map()



%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Diffusion-map
%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_Diffusion_map()



%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate ADC-map
%%%%%%%%%%%%%%%%%%%%%%%%%%
function map_out=l_ADC_map(data,b,opt,Dat)
% Calculate Apparent Diffusion Coefficients
% The fitted functions:
%
% Linear:    log(S)=log(S0)-b*ADC
% Nonlinear: S=S0*exp(-b*ADC)
%
% where S=measured signal
%       S0 = "spin density"
%       b = diffusion weighting
%       ADC = apparent diffusion coefficient

if ~exist('opt','var')
  % Default to linearized fit
  opt = '';
end

b=b(:);

%% Data size
sz=size(data);

%% Data slice index matrix
IndMtx = reshape(1:sz(3),[length(b) Dat.nMaps])';

ADC = zeros([sz(1) sz(2) Dat.nMaps]);
S0 = zeros([sz(1) sz(2) Dat.nMaps]);

% Intensity values should not be less than zero
data(data<0)=0;

%% Use linearized form (fast)
if isempty(opt) || strcmpi(opt,'linear')
  H = [ones(size(b)) -b];
  for ii=1:Dat.nMaps
    for kk=1:sz(2)
      tmp=squeeze(data(:,kk,IndMtx(ii,:))).';
      z=log(tmp);
      th=H\z;
      S0(:,kk,ii)=exp(th(1,:)');
      ADC(:,kk,ii)=th(2,:)';
    end
  end
  ADC(ADC<0)=0;
  ADC(isinf(ADC))=0;
  ADC(isnan(ADC))=0;
else
  %% Fit using nonlinear LS (slow, but more accurate)
  
end

map_out.Map = ADC;
map_out.S0 = S0;
map_out.FitValues = b;





%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate MT-maps
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [map]=l_MTmap(data,fit_vals,nMaps)

sz=size(data);
lims = [-700 700];



[sorted_vals,ind]=sort(fit_vals);
[mx,mx_ind]=max(abs(sorted_vals));
sorted_vals(mx_ind)=[];
sorted_vals = sorted_vals(3:end-2);


%% Allocate space for map
map=zeros(sz(1),sz(2));

lim1 = find(sorted_vals==lims(1));
lim2 = find(sorted_vals==lims(2));


if false

for ii=1:sz(1)
  for kk=1:sz(2)
    vox_data = squeeze(data(ii,kk,:));
    vox_data=vox_data(ind);
    vox_data=vox_data./vox_data(mx_ind);
    vox_data(mx_ind)=[];
    vox_data = vox_data(3:end-2);
    
    df=diff([vox_data(lim1) vox_data(lim2)]);
    map(ii,kk)=df;
  end
end

else

zind=find(sorted_vals==0);
fit=sorted_vals(zind-3:zind+3);
fit=fit(:);
%fit(4)=[];
H = [fit.^2 fit ones(size(fit))];

for ii=1:sz(1)
  for kk=1:sz(2)
    vox_data = squeeze(data(ii,kk,:));
    vox_data=vox_data(ind);
    vox_data=vox_data./vox_data(mx_ind);
    vox_data(mx_ind)=[];
    vox_data = vox_data(3:end-2);
    
    %% Do peak picking by fitting a polynomial
    z=vox_data(zind-3:zind+3);
    %z(4)=[];
    th=H\z;
    zz=fit(1):fit(end);
    
    [mn,mn_ind]=min(polyval(th,zz));
    shift_val=zz(mn_ind);
    %map(ii,kk)=shift_val;
    %continue
    
    
    %% Do a spline interpolation to the data
    %interp_freq = sorted_vals(1):1:sorted_vals(end);
    
    interp_data = pchip(sorted_vals,vox_data,...
                        lims-shift_val);
    
    
    %% Shift frequency values
    %interp_freq = interp_freq-shift_val;
    
    % Get indices to limits
    %lim1 = find(interp_freq==lims(1));
    %lim2 = find(interp_freq==lims(2));
    
    try
      df=diff([interp_data]);
    catch
      df=1;
    end
    if abs(df)<0.5
      map(ii,kk)=df;
    end
    
    
% $$$     if ii==70 && kk==42
% $$$       plot(fit,z,'*-',zz,polyval(th,zz),'r')
% $$$       shift_val
% $$$       pause
% $$$     end
% $$$     
% $$$     if abs(interp_freq(mn_ind))<200
% $$$     
% $$$       % Shift interpolated data
% $$$       interp_freq = interp_freq-interp_freq(mn_ind);
% $$$       
% $$$       % Get indices to limits
% $$$       lim1 = find(interp_freq==lims(1));
% $$$       lim2 = find(interp_freq==lims(2));
% $$$       
% $$$       df=diff([interp_data(lim1) interp_data(lim2)]);
% $$$       if abs(df)<0.5
% $$$         map(ii,kk)=df;
% $$$       end
% $$$     end
  end
  %disp(num2str(ii))
end

end
