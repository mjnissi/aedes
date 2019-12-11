function jdx = aedes_readjcamp(filename)
% AEDES_READJCAMP - Read JCAMP DX format files (Bruker parameter files)
%   
%
% Synopsis: 
%       jdx=aedes_readjcamp(filename)
%
% Description:
%       The function reads the JCAMP DX files and returns a
%       structure with parameters as structure fields. The input
%       argument is a string containing the full path to the file.
%
% Examples:
%       jdx=aedes_readjcamp('C:\path\to\jcamp_dx_file')
%
% See also:
%       AEDES_READBRUKER, AEDES_DATA_READ, AEDES

jdx = [];

% Prompt for a file if not given as an input argument
if nargin == 0
	[fn,fp] = uigetfile({'*.*','All Files (*.*)'},'Open a JCAMP DX file');
	if isequal(fn,0)
		return
	end
	filename = [fp,fn];
elseif nargin > 1
	error('Too many input arguments.');
end

% Open the file for reading
fid = fopen(filename,'r');
if fid < 0
	error('Could not open file "%s" for reading.',filename);
end

% Check that the file is a JCAMP DX file
str = fread(fid,20,'char=>char');
if isempty(regexp(str.','^\s*##TITLE'))
	fclose(fid);
	error('File "%s" is not a valid JCAMP DX file.',filename)
end
fseek(fid,0,-1); % Rewind file

C = fread(fid,inf,'char');
fclose(fid);

% Remove carriage returns
C(C==13)=[];

% Convert to string
C = char(C.');

% Remove comment lines
C = regexprep(C,'\$\$([^\n]*)\n','');

% Remove unnecessary line breaks
f = @l_RemoveLineBreaks;
C=regexprep(C,'^(\s*[^#].*?)(?=\n\s*#)','${f($1)}','lineanchors');
C=regexprep(C,'(\([^\)]+?)\n(.*?\))','${f([$1,$2])}','lineanchors');
CC = regexp(C,'\s*##','split');
CC(1)=[];

% Parse the file line-by-line
for ii=1:length(CC)
	
	str = CC{ii};
	if strncmp(str,'END=',4)
		continue
	end
	
	% The commented regexp sometimes fails with long strings...
	%param = regexp(str,'^(.*)=','tokens','once');
	ind = find(str==61); % Find '=' chars...
	if isempty(ind)
		param='';
	else
		param=str(1:ind(1)-1);
	end
	%param = strrep(param{1},'$','');
	param = strrep(param,'$','');
	param = l_CheckParameter(param);
	
	if any(str==sprintf('\n'))
		% Get size
		sz = regexp(str,'=\s*\((.*?)\)\s*\n','tokens','once');
		sz = str2num(['[',sz{1},']']);
		
		% Parse value
		value = regexp(str,'\n(.*)$','tokens','once');
		value = value{1};
		value = l_CheckValue(value,sz);
	else
		value = regexp(str,'=\s*(.*)','tokens','once');
		value = value{1};
		value = l_CheckValue(value);
	end
	
	% Add to structure
	jdx.(param) = value;
	
end





% ==========================
% - Subfunctions -
% ==========================

% - Remove linebreaks
function out = l_RemoveLineBreaks(str)

out = strrep(str,sprintf('\n'),'');


% - Check parameter value --------------------------
function out = l_CheckValue(val,sz)

if nargin == 1
	sz = 0;
end

% Remove insignificant whitespace
val = strtrim(val);

if isempty(val)
	out = val;
	return
end

% Handle strings and string lists
if val(1) == '<' && val(end) == '>'
	val(val=='<')='''';
	val(val=='>')='''';
	out = eval(['{',val,'}']);
	if length(out) == 1
		out = out{1};
	end
	return
end

% Handle cell matrices
if val(1) == '(' && val(end) == ')'
	nRows = length(find(val==')'));
	
	% Nested tables are not supported. This is a workaround for nested tables
	% and everything is read in a single lined table...
	if nRows ~= sz && sz>0
		nRows=sz;
	end
	
	val(1) = '';
	val(end) = '';
	val(val=='(')='';
	val(val==')')=',';
	val(val=='<')='';
	val(val=='>')='';
	
	% Split using the commas
	val_split = regexp(val,',\s+','split');
	val_out = cell(size(val_split));
	
	% Try to convert to numbers
	for ii = 1:length(val_split)
		num = str2double(val_split{ii});
		if isnan(num)
			val_out{ii} = val_split{ii};
		else
			val_out{ii} = num;
		end
	end
	
	
	out = reshape(val_out,[],nRows).';
	return
end

% Check if the string contains only numbers before tryin to convert to a
% number. str2num uses eval command and if the string matches to a
% function name strange things can happen...
tmp2 = regexp(val,'[^\d\.\seE-+]');
if ~isempty(tmp2)
	out = val;
	return
end

% Convert value to numeric if possible
tmp = str2num(val);
if ~isempty(tmp) && isreal(tmp)
	if length(sz)>1
		tmp = reshape(tmp,sz(2),sz(1),[]);
		tmp = permute(tmp,[2 1 3]);
	end
	out = tmp;
	return
end

out = val;

% - Check parameter strings -------------------------
function out = l_CheckParameter(param)

alphabets = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ';
numbers = '1234567890';

% Remove insignificant whitespace
param = strtrim(param);

if isempty(param)
	out = 'EMPTY_PARAM';
	return
end

% Check parameter starts with a valid structure field character
if ~any(param(1)==alphabets)
	param = ['PAR_',param];
end

% Check that the parameter string does not contain any illegal characters
% (for Matlab structure fields)
ind = ~ismember(param,[alphabets,numbers,'_']);
if any(ind)
	param(ind) = '_';
end

out = param;
