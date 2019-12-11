function filter_fmri_data()

if nargin==0
  % Prompt for the input file
  [fname,fpath,findex] = uigetfile({'*.nii','4D NIfTI Files (*.nii)';...
	'*.*','All Files (*.*)'},'Select 4D NIfTI file');
  if isequal(fname,0)
	% Canceled
	return
  end
  infile = fullfile(fpath,fname);
end

if nargin==0
  % Prompt for the output file
  [fn,fp] = uiputfile({'*.nii','NIfTI files (*.nii)';...
	'*.*','All Files (*.*)'},'Select output file',...
	'4D_nifti.nii');
  if isequal(fn,0)
	% Canceled
	return
  end
  outfile = fullfile(fp,fn);
end

% Read the data
data = aedes_read_nifti(infile);
if length(size(data.FTDATA))<4
  error('Input data has to be 4D data!')
end

% High-pass cutoff
hp_co = 0.009;
lp_co = 0.085;
TR = 1.7;

% High pass filter
[B_hi,A_hi]=butter(8,hp_co/((1/TR)/2),'high');

% Low-pass filter
[B_lo,A_lo]=butter(15,lp_co/((1/TR)/2),'low');

% Waitbar
wbh=aedes_wbar(0,'Filtering image data...');

counter=1;
sz=size(data.FTDATA);
nVox = prod(sz(1:3));
for ii=1:size(data.FTDATA,1)
  for kk=1:size(data.FTDATA,2)
	for zz=1:size(data.FTDATA,3)
	  tmp = squeeze(double(data.FTDATA(ii,kk,zz,:)));
	  data.FTDATA(ii,kk,zz,:) = filtfilt(B_lo,A_lo,filtfilt(B_hi,A_hi,tmp))+mean(tmp);
	  aedes_wbar(counter/nVox,wbh)
	  counter=counter+1;
	end
  end
end
close(wbh)


% Write new 4D NIfTI  image
cwh = aedes_calc_wait(sprintf('Writing 4D NIfTI file:\n%s',outfile));
aedes_write_nifti(data,outfile);
delete(cwh)


