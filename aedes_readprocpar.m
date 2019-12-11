function [procpar,msg] = aedes_readprocpar(filename)
% AEDES_READPROCPAR - Read VNMR (Varian) PROCPAR-file
%   
%
% Synopsis: 
%       [procpar,msg]=aedes_readprocpar(filename)
%
% Description:
%       The function reads the VNMR parameter file and returns a
%       procpar-structure with parameters as structure fields. The input
%       argument is a string containing either the full path to the
%       procpar file or the .fid-directory containing the procpar file.
%
%       If an error occurs, the first output argument will be empty
%       (i.e. procpar=[];) and the second output argument msg will
%       contain the error message. Msg is an empty string if no error has
%       occurred.
%
%       NOTE: AEDES_READPROCPAR adds an additional array field in the
%       procpar-structure. The parameter procpar.arrayjoint does NOT
%       exist in the original procpar-file. The arrayjoint field contains
%       only a flag that can be used to determine if multiparameter
%       arrayed acquisition is arrayed jointly (1=parameters arrayed
%       jointly, 0=parameters not arrayed jointly).
%
% Examples:
%       procpar=aedes_readprocpar('C:\MyData\MyAcquisition.fid\procpar')
%
%       or
%
%       procpar=aedes_readprocpar('C:\MyData\MyAcquisition.fid')
%
% See also:
%       AEDES_READFID, AEDES_DATA_READ, AEDES, AEDES_READTAB

% This function is a part of Aedes - A graphical tool for analyzing 
% medical images
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


% Default values for output arguments
procpar=[];
msg='';


%% Ask for a file if not given as an input argument
if nargin==0
  [fn,fp,fi]=uigetfile({'PROCPAR;procpar','Varian PROCPAR files';...
                   '*.*','All Files (*.*)'},...
                       'Select PROCPAR file');
  if all(fn==0) | all(fp==0)
    % Cancel pressed
    return
  end
  filename = [fp,fn];
end

% Parse filename
[fpath,fname,fext]=fileparts(filename);
if ~strcmpi(fname,'procpar')
  if isempty(fname)
    fpath = [fpath,filesep];
  else
    fpath = [fpath,filesep,fname,fext,filesep];
  end
  fname = 'procpar';
else
  fpath = [fpath,filesep];
end

%% Open file
fid=fopen([fpath,fname],'r','ieee-be');
if fid<0
  msg=['Could not open procpar file "',fpath,fname,'" for reading'];
  return
end


%% Read whole file into cell array
try
	try
		C = textscan(fid,'%s','delimiter','\n','BufSize',1024*1024);
	catch
		% Bufsize returns an error in R2014b->
		C = textscan(fid,'%s','delimiter','\n');
	end
  procpar_str=C{1};
catch
  msg={'Error while reading procpar file.','',...
       lasterr};
  fclose(fid);
  return
end

%% Close file
fclose(fid);

%% Parse procpar file lines
%try
  nonlabelchars='123456789';
  field_chars = '_1234567890abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ';
  nonfield_chars = char(setdiff(0:256,double(field_chars)));
  for ii=1:length(procpar_str)
    if any(procpar_str{ii}(1)==nonlabelchars)
      %% Read info
      str=deblank(procpar_str{ii});
      ind=find(str==' ');
      if isempty(ind)
        ind=0;
      end
      str=str(ind(1)+1:end);
      
      if str(1)=='"'
        ind2=find(str=='"');
        for kk=1:2:(length(ind2))
          procpar.(label){end+1} = strrep(str(ind2(kk):ind2(kk+1)),'"','');
        end
      else
        procpar.(label) = str2num(str);
      end
    elseif procpar_str{ii}(1)=='"' % Read string from line
      str=deblank(procpar_str{ii});
      procpar.(label){end+1} = strrep(str,'"','');
    elseif procpar_str{ii}(1)=='0' % Empty line, end of block
      continue
    else
      %% Read label
      ind=find(procpar_str{ii}==' ');
      label=procpar_str{ii}(1:ind-1);
      
      %% Make sure that the characters in the label are compatible with 
      %% Matlab structure fields. If not, replace them with underscore...
      ind2=ismember(label,nonfield_chars);
      if any(ind2)
        label(ind2)='_';
      end
      procpar.(label)={};
    end
  end
  
  %% Parse array parameter. If more than one parameters are arrayed, they
  %  are separated with commas. If parameters are jointly arrayed they
  %  are enclosed in brackets.
  if isfield(procpar,'array') && ~isempty(procpar.array{1})
    str = procpar.array{1};
    
	% Add Matlab cell characters
	str = strrep(strrep(str,'(','{'),')','}');
	str = ['{',str,'}'];
	
	% Add string characters around words
	str=regexprep(str,'(\w+)(,|\})','''$1''$2');
	
	% Evaluate to formulate s cell
	procpar.array = eval(str);
	
  end
%catch
%  msg=['Error while parsing procpar lines'];
%  procpar=[];
%  return
%end


