function [tablib,msg] = aedes_readphasetable(filename)
% AEDES_READPHASETABLE - Read VNMR table file
%   
%
% Synopsis: 
%       [tab,msg]=aedes_readphasetable(filename)
%
% Description:
%       The function reads and interprets a Varian phase table file
%       "filename" given as a string input argument. The phase order is
%       returned in the output argument "tab". The second output
%       argument contains the possible error message. If no error
%       occurred, the output argument "msg" is empty.
%
% Examples:
%       [tab,msg]=aedes_readphasetable('64alt2k');
%
% See also:
%       AEDES_READFID, AEDES_READPROCPAR, AEDES

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


msg='';
tablib=[];

%% Open table file
fid=fopen(filename,'r','ieee-be');
if fid<0
  tablib=[];
  msg=['Could not open file "' filename '" for reading.'];
  return
end

%% Scan file
C=textscan(fid,'%s','delimiter','\n');
if isempty(C)
  fclose(fid)
  tablib=[];
  msg=['Error while reading table file "' filename '".'];
  return
end

%% Close file
fclose(fid);
tab_str=C{1};

numchars='1234567890';
count_param=1;
count_val=1;
param={};
tablib={};
vals={};
tmp=[];
for ii=1:length(tab_str)
  str=deblank(tab_str{ii});
  
  % Skip empty lines
  if isempty(str) || all(str==32)
    continue;
  end
  
  % Check parameter name
  ind=find(str=='=');
  if isempty(ind)
    %tmp(count_val,:)=str2num(str);
    tmp=[tmp; l_ParseLine(str)];
    %count_val=count_val+1;
  else
    
    %% Write buffered values
    if ~isempty(param)
      vals{end+1}=tmp;
      count_val=1;
      tmp=[];
    end
    param{count_param}=str(1:ind-2);
    count_param=count_param+1;
    
    % Check if the line contains only the parameter
    if ~(str(end)=='=')
      %tmp(count_val,:)=str2num(str(ind+1:end));
      tmp=[tmp l_ParseLine(str(ind+1:end))];
      %count_val=count_val+1;
    end
  end
end

%% Write buffered values
if ~isempty(param)
  vals{end+1}=tmp;
end

%% Construct tablib cell array
tablib=cell(length(param),2);
tablib(:,1)=param;
tablib(:,2)=vals;



%%%%%%%%%%%%%%%%%%%%%
% Parse line
%%%%%%%%%%%%%%%%%%%%%
function out=l_ParseLine(str)


%% Check if possible shorthand indentifiers are found
repall=[];
vals = [];
if any(str=='{') || any(str=='(') || any(str=='[')
  
  ind=find(str=='}');
  if ~isempty(ind)
    repall=str2num(str(ind+1:end));
    str=str(1:ind);
    str=strrep(str,'}','');
    str=strrep(str,'{','');
  end
  
  if ~any(str=='(') && ~any(str=='[')
    vals = str2num(str);
  else
    
    ind_paren=find((str==')')+(str=='('));
    ind_paren=reshape(ind_paren,2,length(ind_paren)/2).';
    ind_brackets=find((str==']')+(str=='['));
    ind_brackets=reshape(ind_brackets,2,length(ind_brackets)/2).';
    for ii=1:length(str)
      if any(ii==ind_paren(:,1))
        ind=find(ii==ind_paren);
        vals = [vals l_ParseParen(str,ind_paren(ind,1),ind_paren(ind,2))];
      elseif any(ii==ind_brackets(:,1))
        ind=find(ii==ind_brackets);
        vals = [vals l_ParseBrackets(str,ind_brackets(ind,1),ind_brackets(ind,2))];
      end
    end
  end
else
  out = str2num(str);
  return
end

if ~isempty(repall)
  tmp=repmat(vals,repall,1);
  out = tmp(:).';
else
  out = vals;
end

function vals=l_ParseParen(str,ind1,ind2)

%% Parse shorthand parenthesis

% Get multiplier
tmp=find(str(ind2+1:end)==' ');
if isempty(tmp)
  mplier = str2num(str(ind2+1:end));
else
  mplier = str2num(str(ind2+1:ind2+tmp(1)));
end

% Get values
tmp=str2num(str(ind1+1:ind2-1));
vals=repmat(tmp,1,mplier);


function vals=l_ParseBrackets(str,ind1,ind2)


%% Parse shorthand brackets
% Get multiplier
tmp=find(str(ind2+1:end)==' ');
if isempty(tmp)
  mplier = str2num(str(ind2+1:end));
else
  mplier = str2num(str(ind2+1:ind2+tmp(1)));
end

% Get values
tmp=str2num(str(ind1+1:ind2-1));
tmp2=repmat(tmp,mplier,1);
vals=tmp2(:).';
