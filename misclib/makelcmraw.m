function ok = makelcmraw(DATA,outdir,outname)
% MAKELCMRAW - Write data in LCModel RAW file format
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


%% Default RAW file name and path
FilePath = [pwd,filesep];
FileName = 'out.RAW';

%% Parse input arguments
if nargin==0
  %% Ask for a file
  [DATA,msg] = aedes_readfid('','return',3);
  if isempty(DATA)
    error(msg)
    return
  end
elseif nargin==2
  FilePath = outdir;
  if ~(FilePath(end)==filesep)
    FilePath = [FilePath,filesep];
  end
elseif nargin==3
  FilePath = outdir;
  if ~(FilePath(end)==filesep)
    FilePath = [FilePath,filesep];
  end
  FileName = outname;
  [fp,fn,fe]=fileparts(FileName);
  FileName = [fn,'.RAW'];
end

%% Get default RAW file name
if any(nargin==[0 1 2])
  try
    fpath=DATA.HDR.fpath;
    ind = find(fpath==filesep);
    if ind(end)==length(fpath)
      tmp=fpath(ind(end-1)+1:end-1);
    else
      tmp=fpath(ind(end)+1:end);
    end
    FileName = [tmp(1:end-4),'.RAW'];
  catch
    % Failed, issue a warning and use out.RAW instead
    warning('Could not determine filename. Using "out.RAW"')
    FilePath = [pwd,filesep];
    FileName = 'out.RAW';
  end
end


%% Open file for writing
fid = fopen([FilePath,FileName],'w');
if fid<0
  error('Could not open file "%s" for writing',FileName)
end

%% Write SEQPAR
fprintf(fid,' $SEQPAR\n');
fprintf(fid,' hzpppm=%3.3f\n',DATA.PROCPAR.sfrq);
fprintf(fid,' echot=%4.1f\n',DATA.PROCPAR.te*1000);
fprintf(fid,' seq=''%s''\n','STEAM');
fprintf(fid,' $END\n');

%% Write NMID
fprintf(fid,' $NMID\n');
fprintf(fid,' id=''%s''\n',FileName(1:min(length(FileName),20)));
fprintf(fid,' fmtdat=''(2e14.5)''\n');
fprintf(fid,' tramp=1.\n');
fprintf(fid,' volume=1.\n');
fprintf(fid,' $END\n');

%% Write data
real_data = real(DATA.KSPACE);
real_data = real_data(:);
imag_data = imag(DATA.KSPACE);
imag_data = imag_data(:);
tmp=strrep(strrep(strrep(sprintf('   %1.5e   %1.5e\n',[real_data ...
                    imag_data]'),'e+0','e+'),'e-0','e-'),' -','-');
%fprintf(fid,'  %1.5e  %1.5e\n',[real_data imag_data]');
fprintf(fid,'%s',tmp);


%% Close file
fclose(fid);
