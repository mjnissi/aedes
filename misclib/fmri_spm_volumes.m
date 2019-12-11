function fmri_spm_volumes(infile,outdir,varargin)
% This function reads VNMR EPI data and writes volumes in separate NIfTI
% files for SPM analysis
%


if nargin<2
  error('Too few input arguments')
end

% Defaults
if infile(end)==filesep
  infile = infile(1:end-1);
end
[fp,fn,fe] = fileparts(infile);
if strcmpi(fn,'fid')
  [fp,fn,fe] = fileparts(fp);
  filename = fn;
  indir = fp;
else
  filename = fn;
  indir = infile;
end

% If the outout directory does not contain a path
% try to create the output directory into the input
% directory...
if ~any(outdir==filesep)
  [tmp,msg,msg_id] = mkdir(indir,outdir);
  outdir = [indir,filesep,outdir];
end

if ~isdir(outdir)
  error('The output folder %s does not exist!',outdir)
end

VoxelSize = [1 1 1.5 2.039];
XYZUnits = 'mm';
TimeUnits = 'sec';
OmitVolumes = [];
StartVolume = [];
EndVolume = [];

% Parse varargin
if rem(length(varargin),2)~=0
  error('Invalid property/value pairs!')
end

for ii=1:2:length(varargin)
  prop = varargin{ii};
  value = varargin{ii+1};
  
  switch lower(prop)
    case 'filename'
      filename = value;
    case 'voxelsize'
      VoxelSize = value;
    case 'xyzunits'
      XYZUnits = value;
    case 'timeunits'
      TimeUnits = value;
    case 'omitvolumes'
      OmitVolumes = value;
    case 'startvolume'
      StartVolume = value;
    case 'endvolume'
      EndVolume = value;
    otherwise
      error('Unknown property %s',prop)
  end
end

% Check that outdir contains the final file separator
if ~strcmp(outdir(end),filesep)
  outdir = [outdir,filesep];
end


% Read the fMRI data
[data,msg] = aedes_readfid(infile,'precision','single');
if isempty(data)
  fprintf(1,'FMRI_SPM_VOLUMES: Error reading data from %s...\n',infile)
  return
end

% Check if is EPI or RASER data
isEPI = false;
isRASER = false;
if isfield(data,'PROCPAR')
  if isfield(data.PROCPAR,'phaseres') && ...
      isfield(data.PROCPAR,'readres')
    isEPI = true;
  elseif isfield(data.PROCPAR,'teType')
    isRASER = true;
  end
else
  fprintf(1,'FMRI_SPM_VOLUMES: No PROCPAR found (%s)\n',infile)
  return
end


if isempty(StartVolume)
  StartVolume = 1;
end
if isempty(EndVolume)
  EndVolume = size(data.FTDATA,4);
end

% Number of volumes to write
nVols = length(StartVolume:EndVolume);

% Initialize waitbar
wbh = aedes_wbar(0,['Writing volumes...']);

% Write output images
count=0;
for ii=StartVolume:EndVolume
  count=count+1;
  
  aedes_wbar(count/nVols,wbh,'Writing volumes...');
  
  % Check if volume should be omitted
  if ~isempty(OmitVolumes) && any(ii==OmitVolumes)
    continue
  end
  
  % Get volume
  vol = data.FTDATA(:,:,:,ii);
  
  % Rotate and flip volumes VNMR EPI and RASER to right orientations
  %if isEPI
  %  vol = flipdim(aedes_rot3d(vol,1,3),2);
  %elseif isRASER
  %  vol = flipdim(vol,2);
  %end
  
  % Construct file name
  fname = sprintf('%s_%04d.nii',filename,ii);
  
  % Write the NIfTI file
  tmp=aedes_write_nifti(vol,[outdir,fname],...
    'VoxelSize',VoxelSize,...
    'XYZUnits',XYZUnits,...
    'TimeUnits',TimeUnits);
  if ~tmp
    fprintf(1,'Could not write file %s\n',[outdir,fname])
    return
  end
end
close(wbh)







  