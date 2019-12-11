function network_analysis(DATA,ROI,AddInfo)
%
% This Aedes plugin calculates and displays partial correlations from ROIs
%



if AddInfo.isDataMixed
  % Make data a 4D-matrix if it is loaded in aedes as a stack of 2D or 3D images
  DATA2{1}=DATA{1};
  DATA2{1}.FTDATA = zeros([size(DATA{1}.FTDATA,1),...
    size(DATA{1}.FTDATA,2),size(DATA{1}.FTDATA,3),length(DATA)],...
    class(DATA{1}.FTDATA));
  for ii=1:length(DATA)
    DATA2{1}.FTDATA(:,:,:,ii)=DATA{ii}.FTDATA;
  end
  DATA=DATA2;
  CurrentVol=AddInfo.CurrentSlice;
else
  CurrentVol=AddInfo.CurrentVol;
end

CurrentVol=AddInfo.CurrentVol;

% Check that there are ROIs defined at all
if isempty(ROI)
  % No ROIs defined, nothing to calculate
  errordlg('No ROIs defined. At least three ROIs must be defined for calculating partial correlations.',...
    'ROIs not defined.','modal');
  return
end

% Check that at least 3 ROIs have been defined
if length(ROI)<3
	errordlg('At least three ROIs must be defined for calculating partial correlations.',...
    'Too few ROIs','modal');
  return
end

% Check that there are ROIs defined in this volume
RoiInds = [];
for ii=1:length(ROI)
  if any(any(any(ROI(ii).voxels{1}(:,:,:,CurrentVol))))
    RoiInds(end+1)=ii;
  end
end

% Return if there were no ROIs in the current volume
if isempty(RoiInds) || length(RoiInds)<3
  errordlg('Too few ROIs defined in the current volume. At least three ROIs must be defined for calculating partial correlations.',...
    'ROIs not defined in current volume.','modal');
  return
end

X = zeros(size(DATA{1}.FTDATA,4),length(RoiInds));

% Get time series from ROIs
for kk=RoiInds
	ind = repmat(ROI(kk).voxels{1}(:,:,:,CurrentVol),[1 1 1 size(DATA{1}.FTDATA,4)]);
	
	% Mean EPI time series
	ts_data = reshape(double(DATA{1}.FTDATA(ind)),[],size(DATA{1}.FTDATA,4));
	ts_data = ts_data.';
	ts_data = mean(detrend(ts_data),2);%+repmat(mean(ts_data),size(ts_data,1),1),2);
	%data_trend=aedes_trendest(double(ts_data),10000);
	X(:,kk) = ts_data;
end

% Calculate partial correlations
[PC,P,CC] = pcorr(X);

% Linear indexes to unique correlations
ind = find(tril(PC,-1)~=0).';

% Display results
fprintf('***********************************\n');
fprintf('* Networks between ROIs\n')
fprintf('***********************************\n');
for ii=ind
	[I,J] = ind2sub(size(PC),ii);
	if P(ii)<=0.01
		fprintf(2,'%s <-> %s (PCC=%.4f, CC=%.4f, PCC/CC=%.04f, P=%.4f)**\n',...
			ROI(J).label,ROI(I).label,PC(ii),CC(ii),PC(ii)/CC(ii),P(ii));
	elseif P(ii) <= 0.05
		fprintf(2,'%s <-> %s (PCC=%.4f, CC=%.4f, PCC/CC=%.04f, P=%.4f)*\n',...
			ROI(J).label,ROI(I).label,PC(ii),CC(ii),PC(ii)/CC(ii),P(ii));
	else
		fprintf(1,'%s <-> %s (PCC=%.4f, CC=%.4f, PCC/CC=%.04f, P=%.4f)\n',...
			ROI(J).label,ROI(I).label,PC(ii),CC(ii),PC(ii)/CC(ii),P(ii));
	end
end
fprintf('***********************************\n\n');




