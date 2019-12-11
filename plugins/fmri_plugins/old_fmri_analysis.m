function old_fmri_analysis(DATA,ROI,AddInfo)
%
% This Aedes plugin does a quick basic fMRI analysis.
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
end

% Prompt for design and contrast
prompt={'Stimulus onsets (scans):','Stimulus durations (scans):',...
  'Contrast vector','Smooth FWHM','TR','High-pass filter cutoff (sec.)',...
  'P-value for FDR correction','Volumes to omit from analysis'};
name='Enter parameters for fMRI analysis';
numlines=1;
defaultanswer={'30 75 120','15 15 15','1 0',...
  '2 2 1','2.039','256','0.05',''};

resp=inputdlg(prompt,name,numlines,defaultanswer);
if isempty(resp)
  % Canceled
  return
end

% Call fmri_analysis in /misclib
onset = str2num(resp{1});
durat = str2num(resp{2});
contr = str2num(resp{3});
if length(contr)~=2
  errordlg('Invalid contrast vector!','Error','modal');
  return
end
smooth_sz = str2num(resp{4});
TR = str2num(resp{5});
hipass = str2num(resp{6});
qFDR = str2num(resp{7});
omitVols = str2num(strrep(resp{8},'end',num2str(size(DATA{1}.FTDATA,4))));

% if isfield(DATA{1},'PROCPAR') && isfield(DATA{1}.PROCPAR,'readres')
%   % Remove reference image from VNMR EPI data
%   DATA{1}.FTDATA = DATA{1}.FTDATA(:,:,:,2:end);
% end

if length(durat)==1
  durat = ones(1,length(onset))*durat;
end

if ~isempty(ROI) && any(strcmpi({ROI(:).label},'mask'))
  roi_ind = find(strcmpi({ROI(:).label},'mask'));
  [maps_out,fdr_th] = fmri_analysis(DATA{1},'TR',TR,...
    'onsets',onset,'durations',durat,...
    'contrast',contr,'smooth',smooth_sz,'FDRth',qFDR,...
    'hipass',hipass,'omitvolumes',omitVols,...
    'mask',ROI(roi_ind).voxels{1}(:,:,:,1));
  ROI(roi_ind)=[];
else
  [maps_out,fdr_th] = fmri_analysis(DATA{1},'TR',TR,...
    'onsets',onset,'durations',durat,...
    'contrast',contr,'smooth',smooth_sz,'FDRth',qFDR,...
    'hipass',hipass,'omitvolumes',omitVols);
end

% Resize tmap to EPI time series size
tmap = maps_out.tmap;
%tmap = repmat(tmap,[1,1,1,size(DATA{1}.FTDATA,4)]);

overlay.FTDATA = tmap;
overlay.ImageOverlayCmapStr = 'hot';
overlay.ImageOverlayTholdDirPos = 1;
overlay.ImageOverlayAlpha = 1;
overlay.fMRIonsets = onset;
overlay.fMRIdurats = durat;
if ~isempty(fdr_th)
  overlay.ImageOverlayThold = fdr_th;
  overlay.ImOverlayMin = 0;
else
  overlay.ImageOverlayAlpha = 0;
end

% If some volumes were omitted, calculate new onsets and durations
if ~isempty(omitVols)
  % create stick function
  ton = onset;
  tof = onset+durat+1;
  tmp=zeros(1,size(DATA{1}.FTDATA,4));
  tmp(ton)=1;tmp(tof)=-1;
  sf=cumsum(tmp);
  sf(omitVols)=[];
  tmp=diff([0 sf]);
  new_onset = find(tmp==1);
  new_durat = find(tmp==-1)-new_onset-1;
  overlay.fMRIonsets = new_onset;
  overlay.fMRIdurats = new_durat;
  DATA{1}.FTDATA(:,:,:,omitVols)=[];
  %overlay.FTDATA(:,:,:,omitVols)=[];
else
  new_onset = onset;
  new_durat = durat;
end


% Estimate BOLD strength for ROIs
if ~isempty(ROI)
  nScans = size(DATA{1}.FTDATA,4);
  fprintf(1,'---------------------------------------\n');
  fh=figure;
	nRows = length(ROI);
	nCols = 1;
	ax_w = 0.9;
	ax_l = (1-ax_w)/2;
	ax_gap = 0.01;
	ax_h = (1-(length(ROI)+5)*ax_gap)/length(ROI);
	bgax = axes('parent',fh,...
		'units','normal','position',[0 0 1 1],...
		'xlim',[0 1],'ylim',[0 1],'visible','off');
%   if length(ROI)<=3
%     nRows = length(ROI);
%     nCols=1;
%   elseif length(ROI)==4
%     nRows = 2;
%     nCols=2;
%   elseif length(ROI)<10
%     nRows = 3;
%     nCols = ceil(length(ROI)/nRows);
%   else
%     nRows = 4;
%     nCols = ceil(length(ROI)/nRows);
%   end
  for kk=1:length(ROI)
    
    % Get time-series indices from ROI
    ind = repmat(ROI(kk).voxels{1}(:,:,:,1),[1 1 1 size(DATA{1}.FTDATA,4)]);

    % Mean EPI time series
    ts_data = reshape(double(DATA{1}.FTDATA(ind)),[],size(DATA{1}.FTDATA,4));
    ts_data = ts_data.';
    ts_data = detrend(ts_data)+repmat(mean(ts_data),size(ts_data,1),1);
    
    % Normalize mean to bold-%
    mean_ts_data = mean(ts_data,2);
    
    % estimate BOLD-% for the ROI mean
    [bold_prc,z,th,z_hat]=estimate_bold(new_onset,new_durat,TR,nScans,mean_ts_data);
    
    z_norm = (z./th(2)-1)*100;
    %z_hat_norm = (z_hat./th(2)-1)*100;
    
    z_hat_norm=aedes_trendest(double(z_norm),10);
    
    % Normalize raw timeseries to bold-%
    %ts_data_bold = (ts_data./repmat(mean(ts_data(1:20,:)),size(ts_data,1),1)-1)*100;
    
    % Plot results
    %ax=subplot(nRows,nCols,kk,'align','parent',fh,'fontsize',8);
		ax=axes('parent',fh,'units','normal',...
			'position',[ax_l 1-kk*ax_h-kk*ax_gap ax_w ax_h],...
			'yaxislocation','right','layer','top','box','on');
    line(1:length(z),z_norm,'color','k',...
      'parent',ax);
    line(1:length(z),z_hat_norm,...
      'color','r','linewidth',2,'parent',ax);
		text(ax_l-0.01,1-kk*ax_h-kk*ax_gap+ax_h/2,...
			['ROI: ',ROI(kk).label],'parent',bgax,...
			'rotation',90,'fontsize',8,...
			'horizontalalign','center',...
			'verticalalign','bottom');
    %title(['Time series for ROI: ',ROI(kk).label],'fontsize',8)
    ylabel(ax,'BOLD-%','fontsize',8);
    set(ax,'xlim',[0 length(z)],...
      'ylim',[min(z_norm)-min(z_norm)*0.05 ...
      max(z_norm)+max(z_norm)*0.05],'fontsize',8);
    for ii=1:length(new_onset)
      xdata = [new_onset(ii) new_onset(ii)+new_durat(ii) ...
        new_onset(ii)+new_durat(ii) new_onset(ii)];
      tmp = get(ax,'ylim');
      ydata = [tmp(1) tmp(1) tmp(2) tmp(2)];
      patch('parent',ax,...
        'xdata',xdata,'ydata',ydata,'FaceColor','r',...
        'FaceAlpha',0.3,'LineStyle','none');
		end
		if kk~=length(ROI)
			set(ax,'xticklabel',[])
		end
		set(ax,'layer','top')
		if kk==length(ROI)
			xlabel('Time (scans)','fontsize',8)
		end
    
    fprintf(1,'BOLD-%% %s: %.3f\n',ROI(kk).label,bold_prc);
	end
	if length(ROI)>=3
		tmp_pos1 = get(0,'screensize');
		tmp_pos2 = get(fh,'position');
		set(fh,'position',[tmp_pos2(1) 100 tmp_pos2(3)*1.5 tmp_pos1(4)-100]);
	end
  fprintf(1,'---------------------------------------\n');
else
  % Estimate BOLD strength for the time-series corresponding to the 
  % maximum value in T-map.
  nScans = size(DATA{1}.FTDATA,4);
  [mx,I_mx]=max(maps_out.tmap(:));
  [II_mx,JJ_mx,KK_mx]=ind2sub(size(maps_out.tmap),I_mx);
  [mn,I_mn]=min(maps_out.tmap(:));
  [II_mn,JJ_mn,KK_mn]=ind2sub(size(maps_out.tmap),I_mn);
  
  ts_data_max = double(squeeze(DATA{1}.FTDATA(II_mx,JJ_mx,KK_mx,:)));
  ts_data_min = double(squeeze(DATA{1}.FTDATA(II_mn,JJ_mn,KK_mn,:)));
  
  [bold_prc_max,z_max,th_max,z_hat_max]=estimate_bold(new_onset,new_durat,TR,nScans,ts_data_max);
  [bold_prc_min,z_min,th_min,z_hat_min]=estimate_bold(new_onset,new_durat,TR,nScans,ts_data_min);
  
  z_norm_max = (z_max./th_max(2)-1)*100;
  z_hat_norm_max=aedes_trendest(double(z_norm_max),10);
  %z_hat_norm_max = (z_hat_max./th_max(2)-1)*100;
  
  z_norm_min = (z_min./th_min(2)-1)*100;
  z_hat_norm_min=aedes_trendest(double(z_norm_min),10);
  %z_hat_norm_min= (z_hat_min./th_min(2)-1)*100;
  
  % Plot results
  figure;
  ax1=subplot(2,1,1,'align');
  line(1:length(z_max),z_norm_max,'color','k',...
    'parent',ax1);
  line(1:length(z_max),z_hat_norm_max,...
    'color','r','linewidth',2,'parent',ax1);
  title('Time series for Max T-value');
  ylabel('BOLD-%');
  set(ax1,'xlim',[0 length(z_max)],...
    'ylim',[min(z_norm_max)-min(z_norm_max)*0.05 ...
    max(z_norm_max)+max(z_norm_max)*0.05]);

  ax2=subplot(2,1,2,'align');
  line(1:length(z_min),z_norm_min,'color','k',...
    'parent',ax2);
  line(1:length(z_min),z_hat_norm_min,...
    'color','r','linewidth',2,'parent',ax2);
  title('Time series for Min T-value');
  ylabel('BOLD-%');
  set(ax2,'xlim',[0 length(z_min)],...
    'ylim',[min(z_norm_min)-min(z_norm_min)*0.05 ...
    max(z_norm_min)+max(z_norm_min)*0.05]);
  
  for ii=1:length(new_onset)
    xdata = [new_onset(ii) new_onset(ii)+new_durat(ii) ...
      new_onset(ii)+new_durat(ii) new_onset(ii)];
    tmp = get(ax1,'ylim');
    ydata = [tmp(1) tmp(1) tmp(2) tmp(2)];
    patch('parent',ax1,...
      'xdata',xdata,'ydata',ydata,'FaceColor','r',...
      'FaceAlpha',0.3,'LineStyle','none');
    tmp = get(ax2,'ylim');
    ydata = [tmp(1) tmp(1) tmp(2) tmp(2)];
    patch('parent',ax2,...
      'xdata',xdata,'ydata',ydata,'FaceColor','r',...
      'FaceAlpha',0.3,'LineStyle','none');
  end
  
  % Print percents to command window
  fprintf(1,'---------------------------------------\n');
  fprintf(1,'Max T BOLD-%%: %.3f\n',bold_prc_max);
  fprintf(1,'Min T BOLD-%%: %.3f\n',bold_prc_min);
  fprintf(1,'---------------------------------------\n');
end

% Show in Aedes as overlay
aedes(DATA,[],overlay);

if isempty(fdr_th)
  warndlg(['No significant voxels over FDR threshold at p < ',...
    num2str(qFDR),'!'],'No significant voxels!','modal')
end


function [bold_prc,z,th,z_hat]=estimate_bold(onset,durat,TR,num_scans,ts_data)
%% -------------------------------------------------

% Load Rat HRF (8 seconds long) NOTE: DON'T USE FOR HUMAN DATA!!!
bf = [
  0
  0.000008876945146
  0.000387621003097
  0.003012441669394
  0.011548214895067
  0.030056522080168
  0.061233628208021
  0.105349995559045
  0.160158500633323
  0.221528149631163
  0.284405279593501
  0.343761753554314
  0.395323573307322
  0.436002661868442
  0.464045998242252
  0.478965423658889
  0.481327079129854
  0.472473786978796
  0.454237528221769
  0.428680153860941
  0.397883140401152
  0.363793656845973
  0.328124931280142
  0.292303457117353
  0.257453130446147
  0.224406054242129
  0.193730677047637
  0.165769521647024
  0.140680556701329
  0.118477987483032
  0.099069734765398
  0.082290068881472
  0.067926763712502
  0.055742762013807
  0.045492744488986
  0.036935220196777
  0.029840852217657
  0.023997740489123
  0.019214335904674
  0.015320580887648
  0.012167779438633
  0.009627606069476
  0.007590575511228
  0.005964217696074
  0.004671136982604
  0.003647081075751
  0.002839102788446
  0.002203865389816
  0.001706118264261
  0.001317352436400
  0.001014633776781
  0.000779604140824
  0.000597636251764
  0.000457125954290
  0.000348904854607
  0.000265756797153
  0.000202022712087
  0.000153279812313
  0.000116082720651
  0.000087755728843
  0.000066226941806
  0.000049896490408
  0.000037532277385
  ];


%TR=2.039;
%onset=[30 75 120]';
%durat=[15 15 15]';
onset = onset(:);
durat = durat(:);

%bf = spm_hrf(1/16);

% Create stimulus function (32 bin offset)
%k = 165; % Number of scans
k = num_scans;
T     = 16;
dt    = TR/T;
u     = ones(size(onset));
if ~any(durat)
  u  = u/dt;
end
ton = round(onset*TR/dt) + 32;			% onsets
tof = round(durat*TR/dt) + ton + 1;			% offset
sf = zeros((k*T + 128),size(u,2));

for j = 1:length(ton)
  if numel(sf)>ton(j),
    sf(ton(j),:) = sf(ton(j),:) + u(j,:);
  end;
  if numel(sf)>tof(j),
    sf(tof(j),:) = sf(tof(j),:) - u(j,:);
  end;
end
sf        = cumsum(sf);					% integrate

% Convolve stimulus with the HRF
conv_sf = conv(sf,bf);

% Resample the convolved stimulus
R = conv_sf([0:(k-1)]*16+1+32);
R=R(:);

% Remove linear trend but keep the mean
mean_ts_data = mean(ts_data);
ts_data = detrend(ts_data)+mean_ts_data;

% Estimate baseline from the paradigm
bl_ind = find(sf==0);
bl_ind = bl_ind(1:min(length(bl_ind),20));
if isempty(bl_ind)
  baseline = mean(ts_data);
else
  baseline = mean(ts_data(bl_ind));
end

mean_ts_data_bold = (ts_data./baseline-1)*100;


% Fit block model to the data
%z=mean_ts_data_bold;
z=ts_data;
H = [R ones(size(R))];
th = H\z;
z_hat = H*th;

% Calculate results
bold_prc = ((th(1)*max(R))/th(2))*100;

