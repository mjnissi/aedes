function smooth_data_out = aedes_fmri_smooth(data,fwhm_sz,voxsize,out_file)

if nargin==4
  writeSmoothedData = true;
else
  writeSmoothedData = false;
end

if nargin<3 || isempty(voxsize)
  fwhmInPixels = true;
else
  fwhmInPixels = false;
end

% Check if data is an Aedes structure
if isstruct(data)
  data = data.FTDATA;
elseif ischar(data)
  data = aedes_data_read(data);
  data = data.FTDATA;
end

UseImFilter = true;

smooth_data = zeros(size(data));

% Check data dimensions
if any(ndims(data)==[4,3]) && size(data,3)~=1
  use3Dkernel = true;
else
  use3Dkernel = false;
end

% Calculate standard deviations using FWHM
if fwhmInPixels
	stds = fwhm_sz/sqrt(8*log(2));
else
  stds = (fwhm_sz/sqrt(8*log(2)))./voxsize;
end

% Calculate kernel size using STDs
kernel_sz = round(6*stds);
if kernel_sz(3)>floor((size(data,3)/2))
  kernel_sz(3)=floor((size(data,3)-1)/2);
end

% Construct the smoothing kernel
if use3Dkernel
  [x,y,z] = meshgrid(-kernel_sz(2):kernel_sz(2),...
	-kernel_sz(1):kernel_sz(1),...
	-kernel_sz(3):kernel_sz(3));
  s_kernel = exp(-(x).^2/(2*(stds(1)).^2)...
	-(y).^2/(2*(stds(2)).^2)...
	-(z).^2/(2*(stds(3)).^2));
  s_kernel = s_kernel/sum(s_kernel(:));
else
  [x,y] = meshgrid(-kernel_sz(2):kernel_sz(2),...
	-kernel_sz(1):kernel_sz(1));
  s_kernel = exp(-(x).^2/(2*(stds(1)).^2)...
	-(y).^2/(2*(stds(2)).^2));
  s_kernel = s_kernel/sum(s_kernel(:));
end

% Smooth the image data -----------------------
nVols = size(data,4);
if ~UseImFilter
  tmp_sz = size(data(:,:,:,1));
  %tmp_sz(3)=tmp_sz(3)+fwhm_sz(3);
  %ind_hi=ceil(tmp_sz/2)+floor(size(s_kernel)/2);
  %ind_lo=ind_hi-(size(s_kernel)-1);
  %ind_lo(ind_lo<1)=1;
  %tmp_kernel = zeros(tmp_sz);
  %tmp_kernel((ind_lo(1):ind_hi(1))+1,(ind_lo(2):ind_hi(2))+1,...
  %  (ind_lo(3):ind_hi(3)))=s_kernel;
	tmp_kernel = s_kernel;
	tmp_kernel(tmp_sz(1),tmp_sz(2),tmp_sz(3))=0; % Pad with zeros
end
for ii=1:nVols
  if ii==1
    fprintf(1,'Smoothing volume %d/%d',ii,nVols);
    bsz = length(sprintf('%d/%d',ii,nVols));
  else
    fprintf(1,repmat('\b',1,bsz));
    bsz = length(sprintf('%d/%d',ii,nVols));
    fprintf(1,'%d/%d',ii,nVols);
  end
  if UseImFilter
    tmp_data = data(:,:,:,ii);
    smooth_data(:,:,:,ii) = imfilter(tmp_data,s_kernel);
  else
    tmp_data = double(data(:,:,:,ii));
    %pad_sz=size(tmp_data)+fwhm_sz*2;
    %tmp_data(pad_sz(1),pad_sz(2),pad_sz(3))=0;
    %tmp_data = cat(3,tmp_data(:,:,1),tmp_data,tmp_data(:,:,end));
    tmp_smooth = fftshift(ifftn(fftn(tmp_data).*fftn(tmp_kernel)),3);
    smooth_data(:,:,:,ii) = tmp_smooth;%tmp_smooth(:,:,1:size(data,3));
  end
end
fprintf(1,'\n');

if nargout>0
  smooth_data_out = smooth_data;
end

% Write output files
if writeSmoothedData
  aedes_write_nifti(smooth_data,out_file)
end
