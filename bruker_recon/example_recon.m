function [kspace,data,msg_out]=example_recon(kspace,Dat,hdr)
% This is an example of a custom Bruker k-space reconstruction plugin.

% If called without input arguments, return the pluse program names that
% this code reconstructs (matched with hrd.acqp.PULPROG)
if nargin==0
  kspace = {'my_pulse_program1.ppg','my_pulse_program2.ppg'};
  return
end

data=[];
msg_out = '';

% Get parameters
NI = hdr.acqp.NI; % Number of objects
NSLICES = hdr.acqp.NSLICES; % Number of slices
NR = hdr.acqp.NR; % Number of repetitions
phase_factor = hdr.acqp.ACQ_phase_factor; % scans belonging to a single image
im_size = hdr.acqp.ACQ_size;im_size(1)=im_size(1)/2;
order = hdr.acqp.ACQ_obj_order;

% Sort k-space (example...)
kspace = reshape(kspace,im_size(1),phase_factor,NI,im_size(2)/phase_factor,NR);
kspace = permute(kspace,[1 2 4 3 5]);
kspace = reshape(kspace,im_size(1),im_size(2),NI,NR);

% Do fourier transform
data = abs(fftshift(fftshift(fftshift(fft(fft(fft(kspace,[],1),[],2),[],3),1),2),3));

