function brain_masking(DATA,ROI,AddInfo)
% BRAIN_MASKING - Try to mask the (Rat) brain
%   (this is an Aedes plugin)
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


% Prompt for intensity threshold and erode/dilate size
answer = inputdlg({'Intensity threshold (range: 0 - 1)',...
  'Erode/dilate radius (in pixels)'},...
  'Parameters for masking',1,...
  {'0.35','5'});
if isempty(answer)
  % Canceled
  return
end
intensity_th=str2num(answer{1});
if isempty(intensity_th)
  autoThreshold = true;
else
  autoThreshold = false;
end
ErodeDilateSize=str2num(answer{2});

% Mask for morphological operations
se=strel('disk',ErodeDilateSize);
se2 = strel('disk',ErodeDilateSize);

% Initialize ROI
ROI=[];


for kk=1:length(DATA)
  %L = zeros(size(DATA{kk}.FTDATA(:,:,:,AddInfo.CurrentVol)));
  ROI(1).voxels{kk}=false(size(DATA{kk}.FTDATA(:,:,:,AddInfo.CurrentVol)));
  ROI(1).fpath = {''};
  ROI(1).fname = {''};
  ROI(1).label = 'brain_mask';
  ROI(1).color = [255 0 0];
  
  % Calculate index vector for center of mass
  data=double(DATA{kk}.FTDATA(:,:,1,AddInfo.CurrentVol));
  n1=1:size(data,1);
  n2=1:size(data,2);
  N2 = repmat(n2,size(data,1),1);
  N2 = N2(:);
  N1 = n1.';
  N1 = repmat(N1,size(data,2),1);
  C = [N1 N2];
  
  % Number of slices
  nSlices = size(DATA{kk}.FTDATA,3);
  
  %Initiate waitbar
  wbh=aedes_wbar(1/nSlices,sprintf('Estimating mask. Processing slice 1/%d',nSlices));
  
  
  % Work through the 3rd direction
  for ii=1:nSlices
    
    % Update waitbar
    aedes_wbar(ii/nSlices,wbh,sprintf('Estimating mask. Processing slice %d/%d',ii,nSlices))
    
    data_raw=double(DATA{kk}.FTDATA(:,:,ii,AddInfo.CurrentVol));
    data = data_raw;
    
    % Scale data
    data=data./max(data(:));
    
    % Saturate 1% of high pixel values  
    clim=stretchlim(data,[0 0.99]);
    data(data>clim(2))=clim(2);
    data=data./clim(2);
   
    % Median filter data to remove random noise
    data = medfilt2(data,[3 3]);
    
    % Open image  
    data = imerode(data,se);
    data = imdilate(data,se);
    data_ed = data;
    
    % Erode image for obtaining better separation for different structures
    data = imerode(data,se);
      
    % Scale data again
    data=data./max(data(:));
    data_final_erode(:,:,ii) = data;
    
    % Do automatic threshold estimation if requested
    if autoThreshold
      [intensity_th,em] = graythresh(data);
    end
     
    % Threshold image and make binary
    data(data<intensity_th)=0;
    data(data>=intensity_th)=1;
    data=logical(data);
    
    % Fill holes and remove clusters with low number of pixels
    data=bwfill(data,'holes');
    data=bwareaopen(data,20,4);
    
    % Check that we still have something left...
    if ~any(data(:)>0)
      continue
    end
    
    % Label separate regions in the image 
    [L,n]=bwlabel(data,4);
    if n>1
      % Select the region whose center of mass that is closest to the 
      % image center of mass
      tmp=regionprops(L,'Centroid');
      COM_regions=round(reshape([tmp(:).Centroid],2,[]).');
      COM_data = round(data_ed(:).'*C/sum(data_ed(:),'double'));
      COM_data = repmat(COM_data,size(COM_regions,1),1);
      d=sqrt(sum((COM_regions-COM_data).^2,2));
      [mn,ind]=min(d);
      L(L~=ind)=0;
      data=L>0;
    end
    
    % Account for the last erosion
    data=imdilate(data,se2);
    
    % Assign mask to ROI
    ROI(1).voxels{kk}(:,:,ii) = data; 
  end
  
  % Close waitbar
  close(wbh);
end

%assignin('base','data_ed',data_final_erode)

% Open data and mask in Aedes
if length(DATA)>1
  aedes(DATA,ROI);
else
  if size(DATA{1}.FTDATA,4)==1
    aedes(DATA,ROI);
  else
    aedes(DATA{1}.FTDATA(:,:,:,AddInfo.CurrentVol),ROI);
  end
end

