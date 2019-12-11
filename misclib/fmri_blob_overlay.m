function [im_rgb,im,cmap,minmaxT]=fmri_blob_overlay(ImBlob,ImAnat,thold,varargin)
% FMRI_BLOB_OVERLAY - Display fmri blobs over an anatomical image.
%
% Synopsis:
%       [ImRGB,ImIDX,CMap,minmaxT]=fmri_blob_overlay(ImBlob,ImAnat,THold,
%                                          param1,value1,param2,value2,...) 
% 
% Description:
%       The function takes the fMRI activation blob image IMBLOB and the
%       corresponding anatomical image IMANAT and the T-map threshold THOLD
%       and returns the thresholded activation blobs overlayed on the
%       anatomical image as an RGB image (ImRGB), indexed image (ImIDX)
%       together with the corresponding stacked colormap (CMap), and the
%       min and max values of the T-map clim. IMBLOB will be reshaped to
%       match the size of IMANAT if needed.
%
%       ImBlob and ImAnat can be paths to the corresponding data files
%       (NIfTI or VNMR), matrices or Aedes DATA structures. Additional
%       options are given as parameter/value pairs:
%
%       Parameter:              Description:
%       **********              ************
%
%       'Deactivation'         'on','off': Select if deactivation blobs are also
%                               displayed (default = 'off')
%
%       'BlobCMin'              Min color limit for scaling activation
%                               blobs. (default = min(ImBlob))
%
%       'BlobCMax'              Max color limit for scaling activation
%                               blobs. (default = max(ImBlob))
%
%       'AnatCLim'              A two element vector [CLOW CHIGH] for
%                               scaling anatomical image. (default =
%                               min/max)
%
%       'BlobCmap'              String of a valid colormap name
%                               ('jet','hot',etc.) to be used with 
%                               activation blobs. (default = 'hot')
%
%       'BlobNegCmap'           String of a valid colormap name
%                               ('cool','jet',etc.) to be used with 
%                               deactivation blobs. (default = 'cool')
%  
%       'AnatCmap'               String of a valid colormap name
%                               ('jet','hot',etc.) to be used with 
%                               the anatomical image. (default = 'gray')
%
%       'VNMRfix'               'on'/'off'. If set to 'on' the activation
%                               blobs will be rotate 90 degrees and flipped
%                               in left/right orientation. (default = 'off')
%
%
% Examples:
%       [imRGB,imIDX,cmap,minmaxT]=fmri_blob_overlay(...
%           '/path/to/blobs.nii','/path/to/anatomical.nii',3.1); 
%       imagesc(imRGB,minmaxT),axis image,colorbar,colormap('hot')       
%
% See also:
%       FMRI_CORR, FMRI_SMOOTH, FMRI_FILTER

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
blobcmin = [];
blobcmax = [];
blobnegcmin = [];
blobnegcmax = [];
anatcmin = [];
anatcmax = [];
deactivations = false;
blobcmap = 'hot';
blobnegcmap = 'cool';
anatcmap = 'gray';
vnmrfix = false;
voxelsizemin = 0;

% Check number of input arguments
if nargin<3
  error('Too few input arguments!')
end

% Check that the length of varargin is even
if rem(length(varargin),2)~=0
  error('Param/value pair number mismatch!');
end

% Parse varargin ------------------------------
for ii=1:2:length(varargin)-1
  param = lower(varargin{ii});
  value = varargin{ii+1};
  
  switch param
    case 'deactivation'
      if strcmpi(value,'on')
        deactivations = true;
      end
    case 'blobclim'
      blobcmin = value(1);
      blobcmax = value(2);
    case 'blobcmin'
      blobcmin = value;
    case 'blobnegcmin'
      blobnegcmin = abs(value);
    case 'blobnegcmax'
      blobnegcmax = abs(value);
    case 'blobcmax'
      blobcmax = value;
    case 'anatclim'
      anatcmin = value(1);
      anatcmax = value(2);
    case 'blobcmap'
      blobcmap = value;
    case 'blobnegcmap'
      blobnegcmap = value;
    case 'anatcmap'
      anatcmap = value;
    case 'voxelsizemin'
      voxelsizemin = value;
    case 'vnmrfix'
      if islogical(value)
        vnmrfix = value;
      else
        if strcmpi(value,'on')
          vnmrfix = true;
        end
      end
    otherwise
      error('Unknown parameter %s!',param)
  end
end

% Check ImBlob and ImAnat -------------------------
if isstruct(ImBlob) && isfield(ImBlob,'FTDATA')
  ImBlob = ImBlob.FTDATA;
elseif ischar(ImBlob)
  data = aedes_data_read(ImBlob);
  ImBlob = double(data.FTDATA);
  clear data;
elseif ~isnumeric(ImBlob)
  error('First input argument not valid!')
end
if isstruct(ImAnat) && isfield(ImAnat,'FTDATA')
  ImAnat = ImAnat.FTDATA;
elseif ischar(ImAnat)
  data = aedes_data_read(ImAnat);
  ImAnat = data.FTDATA;
  clear data;
elseif ~isnumeric(ImAnat)
  error('Second input argument not valid!')
end


% Crop the anatomical image because the FOV in anatomical and blob image is
% not the same in VNMR EPI-data...
if vnmrfix
  %ImAnat = ImAnat(65:65+127,65:65+127,:,:);
  ImBlob = fliplr(rot90(ImBlob));
end

if length(thold)==2
  if isnan(thold(2))
    thold = thold(1);
  else
    deactivations = true;
  end
end

if deactivations
  tholdneg = abs(thold(2));
  thold = thold(1);
end

% Check image sizes ------------------------------------
ImBlobSz = zeros(1,4);
ImAnatSz = zeros(1,4);
[ImBlobSz(1),ImBlobSz(2),ImBlobSz(3),ImBlobSz(4)] = size(ImBlob);
[ImAnatSz(1),ImAnatSz(2),ImAnatSz(3),ImAnatSz(4)] = size(ImAnat);
if ~isequal(ImBlobSz,ImAnatSz)
  
  % Reshape 2D blobs to anatomical size
  if isequal(ImBlobSz(3:4),[1 1]) && isequal(ImAnatSz(3:4),[1 1])
    ImBlob=imresize(ImBlob,ImAnatSz(1:2));
    warning('Resizing blob image...')
  else
    error('Blob and anatomical images have to be of same size or 2D!')
  end
end

ImBlobNeg = -ImBlob;

% Threshold blobs
ind_lo = find(ImBlob<=thold);
ind_hi = find(ImBlob>thold);
minT = min(ImBlob(:));
maxT = max(ImBlob(:));

if deactivations
  ind_neg_lo = find(ImBlobNeg<=tholdneg);
  ind_neg_hi = find(ImBlobNeg>tholdneg);
  minTneg = min(ImBlobNeg(:));
  maxTneg = max(ImBlobNeg(:));
end
%ImBlob(ImBlob<=thold) = 0;
%ind = find(ImBlob~=0);

% Voxel size threshold
if voxelsizemin>0
  ImBloIndHi = ImBlob>thold;
  ImBloIndHi = bwareaopen(ImBloIndHi,voxelsizemin,8);
  ind_hi = find(ImBloIndHi);
end
if deactivations
  ImBloNegIndHi = ImBlobNeg>tholdneg;
  ImBloNegIndHi = bwareaopen(ImBloNegIndHi,voxelsizemin,8);
  ind_neg_hi = find(ImBloNegIndHi);
end

% Handle clims for blobs and anatomical image -----------
if ~isempty(blobcmin)
  ImBlob(ImBlob<blobcmin) = blobcmin;
  minT = blobcmin;
end
if ~isempty(blobcmax)
  ImBlob(ImBlob>blobcmax) = blobcmax;
  maxT = blobcmax;
end
ImBlob=(ImBlob-min(ImBlob(:)));

if deactivations
  if ~isempty(blobnegcmin)
    ImBlobNeg(ImBlobNeg<blobnegcmin) = blobnegcmin;
    minTneg = blobnegcmin;
  end
  if ~isempty(blobnegcmax)
    ImBlobNeg(ImBlobNeg>blobnegcmax) = blobnegcmax;
    maxTneg = blobnegcmax;
  end
  ImBlobNeg=(ImBlobNeg-min(ImBlobNeg(:)));
end

if ~isempty(anatcmin)
  ImAnat(ImAnat<anatcmin) = anatcmin;
end
if ~isempty(anatcmax)
  ImAnat(ImAnat>anatcmax) = anatcmax;
end

% Construct indexed image -------------------------------
ImBlob = round(ImBlob./(maxT-minT)*255);
if deactivations
  ImBlobNeg = round(ImBlobNeg./(maxTneg-minTneg)*255);
  ImBlobNeg = abs(ImBlobNeg-256);
end
ImAnat = round(ImAnat./max(ImAnat(:))*255);
ImRes = ImAnat;
if deactivations
  ImRes(ind_neg_hi) = ImBlobNeg(ind_neg_hi)+256;
  ImRes(ind_hi) = ImBlob(ind_hi)+2*256;
else
  ImRes(ind_hi) = ImBlob(ind_hi)+256;
end

% Construct joint colormap
if deactivations
  tmp = eval([blobcmap,'(256)']);
  %tmp2 = jet(512);
  %tmp2=tmp2(1:256,:);
  tmp2 = eval([blobnegcmap,'(256)']);
  ImBlobColormap = [tmp2;tmp];
else
  ImBlobColormap=eval([blobcmap,'(256)']);
end
ImAnatColormap=eval([anatcmap,'(256)']);
StackedColormap = [ImAnatColormap;ImBlobColormap];

if nargout>0
  im = ImRes;
  cmap = StackedColormap;
  im_rgb = ind2rgb(im,cmap);
  if deactivations
     minmaxT = [minT maxT -maxTneg -minTneg];
  else
    minmaxT = [minT maxT];
  end
end


