function filtered_data_out = fmri_filter(data,TR,varargin)

% Defaults
hipass = [];
lowpass = [];
detrending = false;
smoothpr = false;
reserve_mean = true;
out_file = '';

% Parse varargin
for ii=1:2:length(varargin)
  param = varargin{ii};
  value = varargin{ii+1};
  switch lower(param)
	case 'hipass'
	  hipass = value;
	case 'lowpass'
	  lowpass = value;
	case 'detrending'
	  if strcmpi(value,'on')
		detrending = true;
    end
    case 'smoothpr'
      if strcmpi(value,'on')
        smoothpr = true;
      end
	case 'reservemean'
	  if strcmpi(value,'off')
		reserve_mean = false;
	  end
	case 'filename'
	  out_file = value;
	otherwise
	  error('Unknown parameter "%s"',param)
  end
end

if isstruct(data)
  data = data.FTDATA;
elseif ischar(data)
  data = aedes_data_read(data);
  data = data.FTDATA;
end
filtered_data = double(data);

if isempty(hipass) && isempty(lowpass) && ~detrending
  return
end

% Construct hi-pass filter
if ~isempty(hipass)
  [B_hipass,A_hipass]=butter(5,hipass/((1/TR)/2),'high');
  %[bz_h,bp_h,bk_h]=butter(15,hipass/((1/TR)/2),'high');
  %[sos_var_high,g_high] = zp2sos(bz_h, bp_h, bk_h);
  %Hd_high = dfilt.df2sos(sos_var_high, g_high);
end

% Construct low-pass filter
if ~isempty(lowpass)
  [B_lowpass,A_lowpass]=butter(10,lowpass/((1/TR)/2),'low');
  %[bz_l,bp_l,bk_l]=butter(15,hipass/((1/TR)/2),'low');
  %[sos_var_low,g_low] = zp2sos(bz_l, bp_l, bk_l);
  %Hd_low = dfilt.df2sos(sos_var_low, g_low);
end


% Filter image data -----------------------------
if ~isempty(lowpass) || ~isempty(hipass) || detrending
  
  do_lp = true;
  do_hp = true;
  if isempty(lowpass)
	do_lp = false;
  end
  
  if isempty(hipass)
	do_hp = false;
  end
  
  % Waitbar
  wbh=aedes_wbar(0,'Filtering image data...');
  
  counter=1;
  sz = ones(1,4);
  sz(1:length(size(data)))=size(data);
  nVox = prod(sz(1:3));
  for ii=1:size(data,1)
	for kk=1:size(data,2)
	  for zz=1:size(data,3)
		tmp = squeeze(double(data(ii,kk,zz,:)));
		mean_tmp = mean(tmp);
		if do_hp
		  tmp = filtfilt(B_hipass,A_hipass,tmp);
      %tmp = filter(Hd_high,tmp);
		  if reserve_mean
			tmp = tmp+mean_tmp;
		  end
		end
		if do_lp
		  tmp = filtfilt(B_lowpass,A_lowpass,tmp);
      %tmp = filter(Hd_low,tmp);
		end
		if detrending
		  tmp = detrend(tmp,'linear');
		  if reserve_mean
			tmp = tmp+mean_tmp;
		  end
		end
		filtered_data(ii,kk,zz,:) = tmp;
		aedes_wbar(counter/nVox,wbh)
		counter=counter+1;
	  end
	end
  end
  close(wbh)
end

if nargout>0
  filtered_data_out = filtered_data;
end

% Write output to file
if ~isempty(out_file)
  aedes_write_nifti(filtered_data,out_file);
end

