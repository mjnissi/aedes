function [maps_out,tmap_th,H] = aedes_fmri_analysis(DATA,varargin)
%
% This function does a quick fMRI analysis in SPM style. No motion
% correction or slice timing at this time...
%

% Defaults
TR = [];
contr = [1];
onset = [];
durat = [];
smooth_kernel = [];
d_cutoff = [];      % Don't temporally filter data by default
qFDR = 0.05;        % Default threshold for FDR correction
globalNorm = false; % Don't do global normalization by default
omitVols = [];      % Don't omit any volumes from the analysis by default
show_wbar = true;   % Show waitbar by default
mask = [];          % Don't mask by default
regr = [];          % Don't use additional regressors by default
hrf_type = 'rodent'; % Use human HRF by default
multicomp = 'FDR';  % Use false discovery rate for multiple comparison correction

if rem(length(varargin),2)~=0
  error('parameters/values missing!')
end

% Parse varargin
for ii=1:2:length(varargin)
  param = lower(varargin{ii});
  value = varargin{ii+1};
  switch param
    case 'tr'
      TR = value;
    case {'onset','onsets'}
      onset = value;
			if ~iscell(onset)
				onset = {onset(:)};
			else
				onset = cellfun(@(x) x(:),onset,'UniformOutput',false);
			end
    case {'durat','durations','duration'}
      durat = value;
			if ~iscell(durat)
				durat = {durat(:)};
			else
				durat = cellfun(@(x) x(:),durat,'UniformOutput',false);
			end
    case 'smooth'
      smooth_kernel = value;
    case 'contrast'
      contr = value;
		case 'multicomp'
			multicomp = value;
    case {'fdrth','p-value'}
      qFDR = value;
    case {'hipass','hicutoff'}
      d_cutoff = value;
    case 'omitvolumes'
      omitVols = value;
    case 'wbar'
      if value==1
        show_wbar = true;
      else
        show_wbar = false;
      end
    case 'mask'
      if isempty(value)
				mask = [];
			elseif isnumeric(value) || islogical(value)
        mask = value;
      else
        error('Invalid mask!')
      end
    case 'regressor'
      regr = value;
		case 'hrf'
			if strcmpi(value,'human') || strcmpi(value,'rodent')
				hrf_type = lower(value);
			else
				error('Unknown HRF type. Valid values for HRF are "human" or "rodent".')
			end
  end
end

% Check that we have something to fit...
if isempty(onset) && isempty(regr)
	warning('No onsets or regressors defined. Only mean will be fitted!')
end

% Check that TR has been given
if isempty(TR)
	error('Repetition time (TR) has not been defined!')
end

% Check onsets and durations
for ii=1:length(onset)
	if ~isempty(onset{ii})
		if isempty(durat{ii})
			error('Durations for onsets have not been defined!')
		end
		if length(durat{ii})==1
			durat{ii} = repmat(durat{ii},1,length(onset{ii}));
			durat{ii} = durat{ii}(:);
		elseif length(durat{ii})~=length(onset{ii})
			error('Mismatch in the number of elements in onset and duration!')
		end
	end
end

% Check data
if isstruct(DATA) && isfield(DATA,'FTDATA')
  data = DATA.FTDATA;
elseif isnumeric(DATA) && ndims(DATA)>3
  data=DATA;
else
  error('Input data has to be a 4D numeric matrix or Aedes data structure!')
end

% Initialize mask
if isempty(mask)
  mask = true(size(data,1),size(data,2),size(data,3));
end

% Check regressors
if ~isempty(regr)
  if size(regr,1)~=size(data,4)
    error('The lengths of the regressors do not match data length!')
	end
  regr = detrend(regr);
end

fprintf(1,'\n******************************************\n');
fprintf(1,'Starting fMRI analysis.\n');
if isstruct(DATA) && isfield(DATA,'HDR') && isfield(DATA.HDR,'fpath')
  filename = [DATA.HDR.fpath,DATA.HDR.fname];
  fprintf(1,'File: %s\n\n',filename)
else
  fprintf(1,'\n');
end

% Omit volumes from data and stimulus function if requested
if ~isempty(omitVols)
 
	for ii=1:length(onset)
		fprintf(1,'Skipping requested volumes...\n');
		if ~isempty(onset{ii})
			% Calculate new onsets and durations
			ton = onset{ii};
			tof = onset{ii}+durat{ii}+1;
			tmp=zeros(1,size(data,4));
			tmp(ton)=1;tmp(tof)=-1;
			sf=cumsum(tmp);
			sf(omitVols)=[];
			tmp=diff([0 sf]);
			new_onset = find(tmp==1);
			new_durat = find(tmp==-1);
			if isempty(new_durat)
				% Block stays up...
				new_durat=length(tmp)-new_onset(end)-1;
			elseif length(new_durat) < length(new_onset)
				new_durat(1:end-1) = new_durat(1:end-1)-new_onset(1:end-1)-1;
				new_durat(end) = length(tmp)-new_onset(end)-1;
			else
				new_durat = new_durat-new_onset-1;
			end
			new_onset=new_onset(:);
			new_durat=new_durat(:);
			
			fprintf(1,['New onsets: ',num2str(new_onset(:)'),'\n']);
			fprintf(1,['New duration: ',num2str(new_durat(:)'),'\n']);
			onset{ii} = new_onset;
			durat{ii} = new_durat;
		end
	end
		
  if ~isempty(regr)
    regr(omitVols,:)=[];
  end
  data(:,:,:,omitVols) = [];
    
end

% Create stimulus functions
CF=[];
k = size(data,4); % Number of scans
for ii=1:length(onset)
	ons = onset{ii};
	dur = durat{ii};
	RR = l_GetStimulusFunction(ons,dur,TR,k,hrf_type);
	CF = [CF RR];
end

% Smooth data if requested
if ~isempty(smooth_kernel) && ~ismember(smooth_kernel,[1 1 1],'rows')
	if show_wbar
		wbh = aedes_calc_wait('Smoothing data...');
		drawnow
	end
	if all(smooth_kernel)
		smooth_data = aedes_fmri_smooth(data,smooth_kernel);
	else
		fprintf(1,'fMRI analysis warning: Could not smooth data!\n');
		smooth_data = data;
	end
	if show_wbar
		close(wbh);
	end
else
	smooth_data = data;
end


if ~isempty(d_cutoff)
  % Create filtering matrix
  nCosine = ceil((2*k*TR)/(d_cutoff + 1));
  S = sqrt(2/k)*cos([1:nCosine]'*pi*([1:k]/k)).';
  KKT = eye(size(S,1))-2*S*S'+S*S'*S*S';
else
  KKT = eye(k);
  S = 0;
end
H = [RR regr ones(k,1)];
H = H-S*(S'*H); % Filter design matrix
nParam = size(H,2);
maps_out = struct('pmap',[],'tmap',[]);
maps_out.pmap = zeros(size(smooth_data,1),size(smooth_data,2),size(smooth_data,3),nParam);
maps_out.tmap = zeros(size(smooth_data,1),size(smooth_data,2),size(smooth_data,3));

% Calculate parametric map(s)
nPlanes = size(smooth_data,3);
nCols = size(smooth_data,2);
nRows = size(smooth_data,1);
if length(contr)<nParam
	contr(nParam)=0; % Pad with zeros
end
c = repmat(contr,nRows,1);
HTH = pinv(H'*H);
R = eye(size(H,1))-H*HTH*H';

if show_wbar
  wbh = aedes_wbar(0,sprintf('Estimating parameters. Processing plane 0/%d',nPlanes));
  drawnow
end

% Process data in columns
for ii=1:nPlanes
  
  %fprintf(1,'Processing plane %d/%d\n',ii,nPlanes);
  for kk=1:nCols
    if show_wbar
      aedes_wbar(ii/nPlanes,wbh,sprintf('Estimating parameters. Processing plane %d/%d, column %d/%d',ii,nPlanes,kk,nCols));
    end
    col_data = squeeze(smooth_data(:,kk,ii,:)).';
    col_data = col_data-S*(S'*col_data);
    th=H\col_data;
    r = col_data-H*th;
    rr=diag(r'*r);
    rr=rr(:);
    sig2 = rr./trace(R*KKT);
    sig2 = repmat(sig2,1,nParam);
    T = diag(c*th)./sqrt(c.*sig2*HTH*H'*KKT*H*HTH*contr');
    T(find(mask(:,kk,ii)==0))=0;
		maps_out.pmap(:,kk,ii,:) = th.';
		maps_out.tmap(:,kk,ii)=T.';
	end
end
if show_wbar
  close(wbh);
end

if show_wbar
  wbh = aedes_calc_wait('Calculating threshold...');
  drawnow
end

% Set NaN:s to zeros
maps_out.tmap(isnan(maps_out.tmap)) = 0;

% Calculate effective degrees of freedom
dof = (trace(R*KKT).^2)/trace(R*KKT*R*KKT);

% p-values
pval_map = 1-aedes_tdist(maps_out.tmap,dof);

if strcmpi(multicomp,'FDR')
	% Perform FDR (False Discovery Rate) correction
	cV = 1;
	
	pValues = pval_map(:);
	tValues = maps_out.tmap(:);

	[pValuesSorted,sortInd] = sort(pValues);
	tValuesSorted = tValues(sortInd);
	nP = length(pValues);
	
	pFDR = [1:nP]'/nP*qFDR/cV; % FDR-correction
	thresFDRind = find(pValuesSorted<=pFDR,1,'last');
	if ~isempty(thresFDRind)
		tmap_th = tValuesSorted(thresFDRind);
	else
		tmap_th = [];
	end
	if ~isempty(tmap_th)
		fprintf(1,['FDR threshold at p<',num2str(qFDR),': %.3f\n'],tmap_th)
	else
		fprintf(1,['No significant voxels survive FDR threshold at p<',num2str(qFDR),'!\n'])
	end
else
	b_inv = betaincinv(2*abs(qFDR-0.5),0.5,dof/2,'lower');
	tmap_th = sqrt(dof*(b_inv/(1-b_inv)));
	nVox = length(find(pval_map<=tmap_th));
	fprintf(1,['Uncorrected threshold at p<',num2str(qFDR),': %.3f\n'],tmap_th)
	if nVox==0
		fprintf(1,['No significant voxels survive uncorrected threshold at p<',num2str(qFDR),'!\n'])
	end
end

if show_wbar
	close(wbh);
end

fprintf(1,'******************************************\n');

% Local function for calculating the HRF
function hrf = l_GetHRF(p,TR,T)

dt  = TR/T;
u   = [0:(p(7)/dt)] - p(6)/dt;
hrf = l_GammaFunc(u,p(1)/p(3),dt/p(3)) - l_GammaFunc(u,p(2)/p(4),dt/p(4))/p(5);
%hrf = hrf([0:(p(7)/TR)]*T + 1);
hrf = hrf'/sum(hrf);

function f=l_GammaFunc(x,h,l)

%           l^h * x^(h-1) exp(-lx)
%    f(x) = ----------------------
%                   gamma(h)

% Calculate shaped gamma function
f = exp((h-1).*log(x) + h.*log(l) - l.*x-gammaln(h));


function RR = l_GetStimulusFunction(ons,dur,TR,k,hrf_type)
		
	% Create stimulus function
	T = 16;
	sf = zeros(k*T,1);
	t_ons = ons*T;
	t_offs = (ons+dur+1)*T;
	sf(t_ons) = 1;
	sf(t_offs) = -1;
	sf = cumsum(sf);
	
	% Get HRF
	if strcmpi(hrf_type,'rodent')
		% Parameters for rodent HRF
		p = [2.3 16 0.34 1 6 0 8];
	else
		% Parameters for human (SPM canonical) HRF
		p   = [6 16 1 1 6 0 32];
	end
	hrf =  l_GetHRF(p,TR,T);
	
	% Convolve with HRF
	sf_conv = conv(sf,hrf);
	
	% Resample to scan times
	RR = sf_conv(1:T:k*T);
	
