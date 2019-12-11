function ROI = aedes_copy_roi(ROI,roi_label,dir_ind,slice_ind,vol_ind,data_ind,copytype)
% AEDES_COPY_ROI - Copy ROI in X, Y, or Z direction
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


% Check input arguments
if nargin<5
  error('Too few input arguments')
end

% Return immediately if ROI is an empty structure
% or if ROI is not structure type
if isempty(ROI) || ~isstruct(ROI)
  return
end

%% Check if ROI is for normal or mixed type data
if length(ROI(1).voxels)>1
  isMixedType = true;
else
  isMixedType = false;
end


% Check if the ROI is identied as an index
% or label
if ischar(roi_label)
  roi_ind = find(strcmp({ROI(:).label},roi_label));
  
  if isempty(roi_ind)
    error(['ROI labeled "' roi_label '" was not found!'])
    return
  end
  
elseif isnumeric(roi_label)
  roi_ind = roi_label;
  if roi_ind<1 || roi_ind>length(ROI)
    error('ROI index exceeds ROI dimensions.')
    return
  end
else
  error('Second input argument not valid!')
  return
end

if isMixedType
  dir_ind = 3;
  
  % Current ROI size
  currentRoiSz = size(ROI(roi_ind).voxels{data_ind});
  currentRoiSlice = ROI(roi_ind).voxels{data_ind};
  
  % Copy "current slice" of the selected ROI
  for ii=1:length(ROI(roi_ind).voxels)
    if any(ii==slice_ind)
      % Resize the ROI slice before adding
      sz = size(ROI(roi_ind).voxels{ii});
      newRoiSlice = imresize(currentRoiSlice,sz,'nearest');
      newRoiInd = find(newRoiSlice);
      if strcmpi(copytype,'append')
        ROI(roi_ind).voxels{ii}(newRoiInd)=true;
      else
        ROI(roi_ind).voxels{ii}(:,:)=false;
        ROI(roi_ind).voxels{ii}(newRoiInd)=true;
      end
    end
  end
  
else
%   if dir_ind==1
%     dir_ind=3;
%   elseif dir_ind==3
%     dir_ind=1;
%   end
  
  % Calculate the projection that is copied
  roi_sz=size(ROI(roi_ind).voxels{1}(:,:,:,vol_ind));
  %tmp=zeros(sz([1 2 3]~=dir_ind),'uint8');
  %keyboard
  if dir_ind == 4
	tmp=find(ROI(roi_ind).voxels{1}(:,:,:,vol_ind)>0);
	[r,c,x]=ind2sub(size(ROI(roi_ind).voxels{1}(:,:,:,vol_ind)),...
	  tmp);
  else
	[r,c]=find(squeeze(sum(ROI(roi_ind).voxels{1}(:,:,:,vol_ind),dir_ind))>0);
  end

  % Generate indices
  lr = length(r);
  r=repmat(r,length(slice_ind),1);
  c=repmat(c,length(slice_ind),1);
  if dir_ind==4
	x=repmat(x,length(slice_ind),1);
  end
  slice_ind = slice_ind(:).';
  slice_ind=repmat(slice_ind,lr,1);
  slice_ind=slice_ind(:);
  v=ones(size(slice_ind))*vol_ind;
  
  
  % Copy to slices
  if dir_ind==3
    % Z direction
    ind = sub2ind(size(ROI(roi_ind).voxels{1}),r,c,slice_ind,v);
  elseif dir_ind==2
    % Y direction
    ind = sub2ind(size(ROI(roi_ind).voxels{1}),r,slice_ind,c,v);
  elseif dir_ind==1
    % X direction
    ind = sub2ind(size(ROI(roi_ind).voxels{1}),slice_ind,r,c,v);
  elseif dir_ind==4
		% V direction
		ind = sub2ind(size(ROI(roi_ind).voxels{1}),r,c,x,slice_ind);
  end
  
	if strcmpi(copytype,'append')
    ROI(roi_ind).voxels{1}(ind) = true;
  else
    ROI(roi_ind).voxels{1}(:) = false;
    ROI(roi_ind).voxels{1}(ind) = true;
  end
end


