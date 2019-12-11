function spectrsum(DATA,fpath)
% SPECTRSUM - 
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


if nargin==0
  DATA = aedes_readfid('','return',3);
  if isempty(DATA)
    % Canceled
    return
  end
  DATA.KSPACE = DATA.KSPACE(4:end,:);
  DATA_out = DATA;
  fpath = DATA.HDR.fpath;
elseif nargin==1
  DATA.KSPACE = DATA.KSPACE(4:end,:);
  DATA_out = DATA;
  fpath = DATA.HDR.fpath;
else
  DATA.KSPACE = DATA.KSPACE(4:end,:);
  DATA_out = DATA;
end


%% Get data around the coline? peak
PeakLimits = [1035 1060];

data = abs(fftshift(fft(linebroad(DATA.KSPACE,0.6672,1)),1));
%data = abs(fftshift(fft(DATA.KSPACE),1));
data = data(PeakLimits(1):PeakLimits(2),:);

[mx,Ind] = max(data);

% Order using average index
TrueInd = Ind(1);
tmp=[];
cols = lines(16);
% Shift spectra
for ii=1:length(Ind)
  shift = -(Ind(ii)-TrueInd);
  if shift~=0
    freqdom = fftshift(fft(DATA.KSPACE(:,ii)));
    freqdom = freqdom(shiftind(length(freqdom),shift));
    DATA_out.KSPACE(:,ii) = ifft(ifftshift(freqdom));
  end
end
DATA_out.tmp = DATA_out.KSPACE;
DATA_out.KSPACE=sum(DATA_out.KSPACE,2);
DATA.tmp = DATA.KSPACE;
DATA.KSPACE = sum(DATA.KSPACE,2);


%% Plot shifted and unshifted sums
fh=figure('units','normal','position',[0.15 0.4 0.7 0.4]);
h=subplot(1,3,1);
plot(abs(fftshift(fft(DATA.KSPACE))),'b')
hold on
plot(abs(fftshift(fft(DATA_out.KSPACE))),'r')
hold off

h=subplot(1,3,2);
%tmp=abs(fftshift(fft(linebroad(DATA.KSPACE,0.6672,5))));
tmp=abs(fftshift(fft(DATA.tmp),1));
plot(tmp(PeakLimits(1):PeakLimits(2),:))
h=subplot(1,3,3);
%tmp=abs(fftshift(fft(linebroad(DATA_out.KSPACE,0.6672,5))));
tmp=abs(fftshift(fft(DATA_out.tmp),1));
plot(tmp(PeakLimits(1):PeakLimits(2),:))
pause

%% Accept/reject
resp = questdlg('Accept or reject?','Accept or reject?',...
                'Accept','Reject','Cancel','Accept');
if strcmpi(resp,'Accept')
  %makelcmraw(DATA_out,fpath)
elseif strcmpi(resp,'Reject')
  %makelcmraw(DATA,fpath)
end
if ishandle(fh)
  close(fh)
end
