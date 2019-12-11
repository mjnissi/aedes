function [kspace,data,msg_out]=cine_recon(kspace,Dat,procpar)
% This is a custom VNMR k-space reconstruction code for CINE data used by 
% aedes_readvnmr.

% If called without input arguments, return the sequence names that
% this code reconstructs
if nargin==0
  kspace = {'cineSSFP'};
  return
end

data=[];
msg_out = '';

ne = procpar.ne;
nv = procpar.nv;
np = procpar.np;
pss = procpar.pss;

kspace=permute(reshape(kspace,np/2,1,1,ne,nv,length(pss)),[1 5 6 4 2 3]);

% Sort data using phasetable ------------------
if Dat.Sorting && ~isempty(Dat.phasetable)
  Dat.phasetable = Dat.phasetable.';
  kspace(:,Dat.phasetable(:),:,:,:)=kspace;
end

% FFT data
data = abs(fftshift(fftshift(ifft(ifft(kspace,[],1),[],2),1),2));
