function [kspace,data,msg_out]=mge3d_recon(kspace,Dat,procpar)
% This is a custom VNMR k-space reconstruction code for MGE3D data used by 
% aedes_readvnmr.

% If called without input arguments, return the list of sequence names that
% this code reconstructs
if nargin==0
  kspace = {'mge3d'};
  return
end

data=[];
msg_out = '';

% Reshape to 4D matrix
kspace = reshape(kspace,[procpar.np/2 procpar.ne procpar.nv ...
	procpar.nv2]);
kspace = permute(kspace,[3 1 4 2]);
data = abs(fftshift(fftshift(fftshift(fft(fft(fft(kspace,[],1),[],2),[],3),1),2),3));
