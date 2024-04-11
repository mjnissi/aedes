function ROIStruct = nifti2roi(fname)
%Converts a nifti file to a .roi -file for aedes
%   Assumes areas of same value = roi


colorarray = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 0.5 0; 1 0 0.5; 1 0.5 0.5; 0 0.5 1; 0.5 0.5 1];
colorarray = round(255*[colorarray; colororder]);

fig = aedes_read_nifti(fname);
roi_vals = unique(fig.FTDATA(:));
roi_vals(roi_vals == 0) = [];
if isempty(roi_vals) || length(roi_vals) > 100
    error('Too few or too many unique values in the nifti file')
end
for ii = 1:length(roi_vals)
    ROI(ii).voxels{1} = fig.FTDATA == roi_vals(ii);
    ROI(ii).fpath = fig.HDR.fpath;
    ROI(ii).fname = fig.HDR.fname;
    ROI(ii).label = num2str(roi_vals(ii));
    if ii< size(colorarray,1)
        ROI(ii).color = colorarray(ii,:);
    else
        ROI(ii).color = round(255*rand(1, 3));
    end
end

DateTime = char(datetime);
FileInfo.DataFileName{1} = fig.HDR.fname;
FileInfo.DataPathName{1} = fig.HDR.fpath;
%RotateFlip.Rotate = 0;
%RotateFlip.Flip = 0;

if nargout == 0
    savename = [extractBefore(fname,'.nii') '.roi'];
    save(savename,"ROI","DateTime","FileInfo");
else
    ROIStruct.ROI = ROI;
    ROIStruct.DateTime = DateTime;
    ROIStruct.FileInfo = FileInfo;
end
end