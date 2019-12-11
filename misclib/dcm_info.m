function dcm_info(pathname,add_params)
% Print some information from DICOM files exported from PACS

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

% Default parameters
Param = {'SeriesDescription',...
  'AcqTime'};
Param_length = [25 8];
SortByAcqTime = true;

if nargin==0 || isempty(pathname)
  try
    tmp_dir=getpref('Aedes','GetDcmDir');
  catch
    tmp_dir = '';
  end
  pathname=uigetdir(tmp_dir);
  if isequal(pathname,0)
    % Canceled
    return
  end
  setpref('Aedes','GetDcmDir',pathname)
elseif nargin>2
  error('Too many input argumetns')
end

% Make sure that pathname is a cell array
if ischar(pathname)
  pathname = {pathname};
end

if nargin==2
  Param = {Param{:} add_params{:}};
  Param_length = [Param_length ones(1,length(add_params))*15];
end
total_width = 7+sum(Param_length)+2*length(Param_length)+6;

% DCM file to read (the first one)
fname = 'i0000_0000b.dcm';

for ii=1:length(pathname)
  % Get the im_(number) folders
  s=dir(pathname{ii});
  fpath = {s([s(:).isdir]).name};

  % Find the im_(number folders)
  tmp=regexp(fpath,'^im_\d{1,2}$');
  ind=~cellfun(@isempty,tmp);
  if any(ind)
    fp={fpath{find(ind)}};
    
    % Sort fp
    [tmp_fp,ind2]=sort(regexprep(fp,'(^.*_)(\d{1})$','$10$2'));
    fp=fp(ind2);
    
    %fprintf(1,'***************************************************\n')
    fprintf(1,'%s\n',repmat('*',1,total_width))
    fprintf(1,'%s\n\n',pathname{ii})
    fprintf(1,'%7s',' ')
    fprintf(1,'%s  ',l_trunkStr(Param{1},Param_length(1)))
    fprintf(1,'%s  ',l_trunkStr(Param{2},Param_length(2)))

    
    for tt=3:length(Param)
      fprintf(1,'%s  ',l_trunkStr(Param{tt},Param_length(tt)));
    end
    fprintf(1,'%s\n',l_trunkStr('nFiles',6))
    fprintf(1,'%s\n',repmat('-',1,total_width))

    dir_info={};
    AcqTime=[];
    for kk=1:length(fp)
      % Get the number of files in the current folder
      s=dir(fullfile(pathname{ii},fp{kk}));
      fn={s(~[s(:).isdir]).name};
      nFiles = length(find(~cellfun(@isempty,regexp(fn,'^i\d{4}_\d{4}b.dcm$'))));

      % Read the dicom headers
      try
        hdr=dicominfo(fullfile(pathname{ii},fp{kk},fname));
      catch
        dir_info{kk}=sprintf('Could not read: %s\n',fullfile(pathname{ii},fp{kk},fname));
        continue
      end
      dir_info{kk}='';
      dir_info{kk} = [dir_info{kk} sprintf('%5s: ',fp{kk})];
      dir_info{kk} = [dir_info{kk} sprintf('%s  ',l_trunkStr(hdr.SeriesDescription,25))];
      dir_info{kk} = [dir_info{kk} sprintf('%s:%s:%s  ',...
        hdr.AcquisitionTime(1:2),...
        hdr.AcquisitionTime(3:4),...
        hdr.AcquisitionTime(5:6))];
      AcqTime(kk)=str2num([hdr.AcquisitionTime(1:2),...
        hdr.AcquisitionTime(3:4),...
        hdr.AcquisitionTime(5:6)]);
      %fprintf(1,'%5s: ',fp{kk});
      %fprintf(1,'%s  ',l_trunkStr(hdr.SeriesDescription,25));
      %fprintf(1,'%s:%s:%s  ',hdr.AcquisitionTime(1:2),...
      %  hdr.AcquisitionTime(3:4),...
      %  hdr.AcquisitionTime(5:6));

      for tt=3:length(Param)
        if not(isfield(hdr,Param{tt}))
          str = 'N/A';
        else
          if isnumeric(hdr.(Param{tt}))
            str = num2str(hdr.(Param{tt})(:).');
          else
            str = num2str(hdr.(Param{tt}));
          end
        end
        dir_info{kk} = [dir_info{kk} ...
          sprintf('%s  ',l_trunkStr(str,Param_length(tt)))];
      end
       dir_info{kk} = [dir_info{kk} ...
         sprintf('(%d)\n',nFiles)];
    end
    if SortByAcqTime
      [tmp,sort_ind]=sort(AcqTime);
      dir_info=dir_info(sort_ind);
    end
    for jj=1:length(dir_info)
      fprintf(dir_info{jj});
    end
    fprintf('%s\n\n\n',repmat('*',1,total_width));
  else
    % Print warning and continue to next iteration
    fprintf(1,'\nWarning: Cannot find im_# folders from\n');
    fprintf(1,'%s\n',pathname{ii});
    fprintf(1,'Skipping folder...\n\n\n')
    continue
  end

end

function str_out=l_trunkStr(str,len)

str_out = fliplr(sprintf(['%',num2str(len),'s'],fliplr(str(1:min(end,len)))));

