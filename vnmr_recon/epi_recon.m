function [kspace,data,msg_out]=epi_recon(kspace,Dat,procpar)
% This is a custom VNMR k-space reconstruction code for EPI data used by 
% aedes_readvnmr.

% If called without input arguments, return the sequence names that
% this code reconstructs
if nargin==0
  kspace = {'epi_se_rapid_sp3','epi','epi_fMRI','epip','epip_fMRI','epiT1rho'};
  return
end

data=[];
msg_out = '';

% Number of receivers
if isfield(procpar,'rcvrs') && ~isempty(procpar.rcvrs)
	nRcvrs = length(procpar.rcvrs{1}=='y');
else
	nRcvrs = 1;
end



% =====================================================
% Reconstruct EPI sequences meaured with the old INOVA
% =====================================================
if isfield(procpar,'readres') && isfield(procpar,'phaseres')
  
  % Number of slices
  tmp_ns=length(procpar.pss);
  
  if isfield(procpar,'navecho') && strcmpi(procpar.navecho{1},'y')
    tmp_nv = procpar.nv-procpar.nseg;
  else
    tmp_nv = procpar.nv;
  end
  kspace = reshape(kspace,[size(kspace,1) ...
    size(kspace,2)/tmp_nv tmp_nv]);
  kspace = permute(kspace,[1 3 2]);
  
  % Reshape to 4D matrix
  kspace = reshape(kspace,[size(kspace,1) size(kspace,2) ...
    tmp_ns size(kspace,3)/tmp_ns]);
	
	
	data = abs(fftshift(fftshift(fft(fft(kspace,[],1),[],2),1),2));
	data = permute(data,[2 1 3 4]);
	data = flipdim(flipdim(data,1),2);
	data(:,:,:,1) = []; % Remove reference image
	
	% ===================================
	% Reconstruct VNMRj and EPIP sequences
	% ===================================
else
	
	if isfield(procpar,'nphase') && ...
			isfield(procpar,'nread')
		isEPIP = true;
	else
		isEPIP = false;
	end
	
	% Number of navigator echos
	if isfield(procpar,'navigator') && strcmpi(procpar.navigator,'y')
		nNav = length(procpar.nav_echo);
	else
		nNav = 0;
	end
	
	% EPI phase correction
	doPhaseCorrection = true;
	if strcmpi(Dat.EPI_PC,'auto')
		% Select phase correction using procpar
		if isfield(procpar,'epi_pc')
			switch lower(procpar.epi_pc{1})
				case 'triple_ref'
					nRef = 3;
				case 'scaled_triple'
					% Issue a warning because unscaled triple ref is used anyway...
					warning(['Using unscaled TRIPLE REFERENCE phase correction although "%s" ',...
						'is indicated in PROCPAR.EPI_PC.'],procpar.epi_pc{1});
					nRef = 3;
				case {'linear','quadratic',...
						'center_pair','pairwise','first_pair'}
					% Issue a warning because pointise correction is used anyway...
					warning(['Using POINTWISE phase correction although "%s" ',...
						'is indicated in PROCPAR.EPI_PC.'],procpar.epi_pc{1});
					nRef = 1;
				case 'pointwise'
					nRef = 1;
				otherwise
					warning('Unknown EPI_PC "%s". Skipping phase correction,',procpar.epi_pc{1})
			end
		else
			% Don't use phase correction
			warning('Cannot find EPI_PC parameter in PROCPAR. Skipping phase correction...');
			doPhaseCorrection = false;
		end
	elseif strcmpi(Dat.EPI_PC,'off')
		nRef = 0;
		doPhaseCorrection = false;
	elseif strcmpi(Dat.EPI_PC,'triple')
		nRef = 3;
	elseif strcmpi(Dat.EPI_PC,'pointwise')
		nRef = 1;
	end
	
% 	% Number of reference images
% 	if isfield(procpar,'epiref_type') && ~isempty(procpar.epiref_type) && ...
% 			strcmpi(procpar.epiref_type,'triple')
% 		nRef = 3;
% 	else
% 		nRef = 1;
% 	end
	
	% Number of segments
	if isfield(procpar,'nseg') && ~isempty(procpar.nseg)
		nSeg = procpar.nseg;
	else
		nSeg = 1;
	end
		
	
	%procpar.image = [1 1 0 -2 -1 1 1 1 1];
	% Get index to first reference image
	ref_start = find(procpar.image == 0);
	ref_start = ref_start(1);
	nDummys = 0;
	if isempty(ref_start) && doPhaseCorrection
		kspace = [];
		data = [];
		msg_out = 'Reference image(s) not found.';
		return
	elseif ref_start(1) > 1
		nDummys = ref_start-1;
		warning('Rejecting %d volumes before first reference image as dummy scans...',nDummys);
	end
	
	% Number of volumes
	nVols = length(procpar.image(nDummys+1:end))-nRef;
	
	% Number of slices
	ns = procpar.ns;
	
	if isEPIP
		nv = procpar.nphase; % Phase sampling
		np = procpar.nread/2; % Read sampling
	else
		nv = procpar.nv;
		np = procpar.np/2;
	end
	
	% Get orientation
	orient = procpar.orient{1};
	
	% Detect partial fourier EPI
	if isfield(procpar,'fract_ky') && procpar.fract_ky~=nv/2
		nv_full = nv;
		nv = nv/2+procpar.fract_ky;
	elseif isEPIP && isfield(procpar,'kzero') && procpar.kzero~=0
		nv_full = nv;
		nv = nv/2+procpar.kzero;
	else
		nv_full = nv;
	end
	
	if isEPIP && (nv_full==64 & np==64)
		% The information about navigator echoes in the procpar doesn't seem to
		% be consistent. This code tried to guess the number of navigators from
		% the k-space size
		nNav = size(kspace,1)/np-nv;
	end
	
	% Calculate indexes for navigator and data echoes
	if nNav > 0
		ind = reshape([1:(nv+nNav*nSeg)],[],nSeg);
		NavInd = ind(1:nNav,:);
		DataInd = ind(nNav+1:end,:);
		NavInd = NavInd(:).';
		DataInd = DataInd(:).';
		
		% Warn that navigator echoes are not used...
		warning('Navigator correction has not been implemented. All navigator echoes are ignored.');
	else
		DataInd = [1:nv];
		NavInd = [];
	end
	
	if Dat.EPIPhasedArrayData
		data = zeros(nv_full,np,ns,nVols,nRcvrs,'single');
	else
		data = zeros(nv_full,np,ns,nVols,'single');
	end
	
	% Look for phase table
	pe_table = [];
	if isfield(Dat,'phasetable') && ~isempty(Dat.phasetable)
		pe_table = Dat.phasetable;
	end
	
	% Reshape kspace and detach reference images
	kspace = reshape(kspace,np,nv+nNav*nSeg,ns,nRcvrs,nVols+nRef+nDummys);
	
	
	ref1_ind = find(procpar.image==0,1);
	ref2_ind = find(procpar.image==-2,1);
	im1_ind = find(procpar.image==-1,1);
	
	if ( isempty(ref2_ind) | isempty(im1_ind) ) && nRef == 3
		warning(['Triple reference phase correction was requested, but ',...
			' < 3 reference images were found. Skipping phase correction.'])
		doPhaseCorrection = false;
	end
		
	if doPhaseCorrection

		ref_data = kspace(:,DataInd,:,:,[ref1_ind ref2_ind im1_ind]);
		if ~isempty(pe_table)
			ref_data(:,pe_table,:,:,:) = ref_data;
		end
		% Flip even kspace lines
		ref_data(:,2:2:end,:,:,:) = flipdim(ref_data(:,2:2:end,:,:,:),1);
		%kspace(:,:,:,:,1:nRef)=[];
		
		% Calculate phase corrections for all receivers
		phase_e = ref_data(:,:,:,:,1);
		rev_phase = ref_data(:,:,:,:,1);
		for ii=1:nRcvrs
			if nRef==1
				%ref_im = ref_data(:,:,:,ii,ref1_ind);
				ref_im = ref_data(:,:,:,ii,1);
				
				phase1 = exp(-sqrt(-1)*angle(fft(ref_im,[],1)));
				rev_phase = 1;
				phase_e(:,:,:,ii,:) = phase1;
				
				%tmp = fft(ref_im,[],1);
				%tmp_max = repmat(max(tmp),[size(ref_im,1),1,1]);
				%phase1 = conj(tmp./tmp_max);
				
				phase_e(:,:,:,ii,:) = phase1;
			elseif nRef==3
				%ref1_ind = find(procpar.image==0,1);
				%ref1 = ref_data(:,:,:,ii,ref1_ind);
				ref1 = ref_data(:,:,:,ii,1);
				phase1 = exp(-sqrt(-1)*angle(fft(ref1,[],1)));
				
				%ref2_ind = find(procpar.image==-2,1);
				%ref2 = flipdim(ref_data(:,:,:,ii,ref2_ind),1);
				ref2 = flipdim(ref_data(:,:,:,ii,2),1);
				phase2 = exp(-sqrt(-1)*angle(fft(ref2,[],1)));
				
				%im1_ind = find(procpar.image==-1,1);
				%im1 = flipdim(ref_data(:,:,:,ii,im1_ind),1);
				im1 = flipdim(ref_data(:,:,:,ii,3),1);
				
				
				rev_phase(:,:,:,ii,:) = fft(im1,[],1).*phase2;
				phase_e(:,:,:,ii,:) = phase1;
			end
		end
	end
	
	% Reconstruct in blocks
	kssz=size(kspace);
	blksz = Dat.EPIBlockSize; % Process EPI data in 100 volume blocks (default)
	nBlocks = ceil(nVols/blksz);
	lnum = length(num2str(nBlocks));
	lnumstr = num2str(lnum);
	bsl = lnum*2+1;
	fprintf(1,'Processing data in blocks of %d volumes\n',blksz)
	fprintf(1,['Processing block...%0',lnumstr,'d/%0',lnumstr,'d'],1,nBlocks);
	
	for ii=1:nBlocks
		fprintf(1,repmat('\b',1,bsl));
		fprintf(1,['%0',lnumstr,'d/%0',lnumstr,'d'],ii,nBlocks);
		tmp_data = [];
		if ii==nBlocks
			vol_inds = (ii-1)*blksz+1+nRef+nDummys:nVols+nRef+nDummys;
		else
			vol_inds = (ii-1)*blksz+1+nRef+nDummys:ii*blksz+nRef+nDummys;
		end
		
		% At the moment navigator echoes are not used...
		%tmp_kspace = kspace(:,nNav+1:end,:,:,vol_inds);
		tmp_kspace = kspace(:,DataInd,:,:,vol_inds);
		if ~isempty(pe_table)
			tmp_kspace(:,pe_table,:,:,:) = tmp_kspace;
		end
		
		% Flip even lines
		tmp_kspace(:,2:2:end,:,:,:) = flipdim(tmp_kspace(:,2:2:end,:,:,:),1);
		
		if doPhaseCorrection
			for kk=1:nRcvrs
				if nRef==3
					for tt=1:length(vol_inds)
						tmp_kspace(:,:,:,kk,tt) = ifft(rev_phase(:,:,:,kk,:)+fft(tmp_kspace(:,:,:,kk,tt),[],1).*phase_e(:,:,:,kk,:));
					end
				else
					for tt=1:length(vol_inds)
						tmp_kspace(:,:,:,kk,tt) = ifft(fft(tmp_kspace(:,:,:,kk,tt),[],1).*phase_e(:,:,:,kk,:));
					end
				end
			end
		end
		
		if nv~=nv_full
			tmp_kspace(:,nv_full,:,:,:) = single(0); 
		end
		tmp_data = fftshift(fftshift(fft(fft(tmp_kspace,[],1),[],2),1),2);
		
		if Dat.EPIPhasedArrayData
			data_block = abs(tmp_data);
			% Permute to correct orientation
			if strcmpi(orient,'trans90')
				data_block = permute(data_block,[2 1 3 5 4]);
				data_block = flipdim(data_block,2);
			else
				data_block = permute(data_block,[1 2 3 5 4]);
			end
			data(:,:,:,vol_inds-nRef,:) = data_block;
		else
			if nRcvrs == 1
				data_block = abs(tmp_data);
			else
				data_block = sqrt(sum(tmp_data.*conj(tmp_data),4));
			end
			% Permute to correct orientation
			if strcmpi(orient,'trans90')
				data_block = permute(data_block,[2 1 3 5 4]);
				data_block = flipdim(data_block,2);
			end
			data(:,:,:,vol_inds-nRef-nDummys) = squeeze(data_block);
		end
		
	end
	fprintf(1,'\n')
end
	











