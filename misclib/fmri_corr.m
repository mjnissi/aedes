function corrmap = fmri_corr(data,seed)

corrmap = [];

if ischar(data)
  data=aedes_data_read(data);
  data = data.FTDATA;
elseif isstruct(data)
  data=data.FTDATA;
end

% Get mean time-series from seed ROIs
seed_struct = l_Roi2seed(data,seed);


% Calculate Pearson correlation coefficients -------------
data_sz = ones(1,4);
data_sz(1:ndims(data)) = size(data);

for kk = 1:length(seed_struct)
 
  % Allocate space for ccc
  ccc = zeros(data_sz(1:3));
  
  T = length(seed_struct(kk).data);
  std_seed = std(seed_struct(kk).data,1);
  norm_seed = seed_struct(kk).data-mean(seed_struct(kk).data);
  seed_plane = reshape(reshape(repmat(norm_seed,...
	[data_sz(1),data_sz(2)]),data_sz(1)*T,[]).',...
	data_sz(1),data_sz(2),[]);
  
  % Loop over in-plane slices
  for ii=1:data_sz(3)
	tmp_data = double(squeeze(data(:,:,ii,:)));
	
	std_data = std(tmp_data,1,3);
	
	ccc(:,:,ii) = ((1/T)*sum(seed_plane.*tmp_data,3))./(std_data.*std_seed);
  end
  
  % Make sure that there are no Infs or nans ...
  ccc(find(isnan(ccc))) = 0;
  ccc(find(isinf(ccc))) = 0;
  
  % Fisher's z-transform to make values normally distributed
  zf = ccc;
  zf = 0.5*log((1+ccc)./(1-ccc));
  
  % Calculate Z-scores
  zscore = zf./(1/sqrt(T-3));
  
  % Construct output structure
  corrmap(kk).ccc = ccc;
  corrmap(kk).zf = zf;
  corrmap(kk).zscore = zscore;
  corrmap(kk).label = seed_struct(kk).label;
end

% Subfunctions -----------------------------------
function seed_struct = l_Roi2seed(data,ROI)
% Thus subfunction extracts mean seed time-series from ROI voxels 
%

seed_struct = [];

if isstruct(ROI)
  for ii=1:length(ROI)
    ind = repmat(ROI(ii).voxels{1}(:,:,:,1),[1 1 1 size(data,4)]);
    seed_struct(ii).data = mean(reshape(double(data(ind)),[],size(data,4)));
    seed_struct(ii).label = ROI(ii).label;
  end
else
  tmp = ROI(:).';
  seed_struct(1).data = double(tmp);
  seed_struct(1).label = 'custom seed';
end



