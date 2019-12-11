function dataformat = aedes_getdataformat(filename)
% AEDES_GETDATAFORMAT - Try to determine the data format of a file
%
%
% Synopsis:
%       dataformat=aedes_getdataformat(filename);
%
% Description:
%       Returns an identifier string corresponding to the data format of
%       the file FILENAME. If the data format cannot be determined, an
%       empty string is returned.
%
%       The possible identifier strings are:
%
%       'vnmr'        <-> Varian FID-file
%       'bruker_raw'  <-> Bruker FID-file
%       'bruker_reco' <-> Bruker reconstructed 2dseq file
%       'nifti'       <-> NIfTI or Analyze 7.5 format file
%       'sur'         <-> S.M.I.S. SUR-File
%       'dcm'         <-> DICOM File
%       'spect/ct'    <-> Gamma Medica SPECT/CT File
%       'mat'         <-> Matlab MAT-File
%       'roi'         <-> Aedes ROI-File
%       'fdf'         <-> Varian FDF-File
%
% Examples:
%
% See also:
%       AEDES, AEDES_DATA_READ

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

% Check number of arguments
if nargin==0
	error('Too few input arguments')
elseif nargin>1
	error('Too many input arguments')
end

dataformat = '';

[f_path,f_name,f_ext] = fileparts(filename);
if isempty(f_path)
	f_path = [pwd,filesep];
else
	f_path=[f_path,filesep];
end

% Check if is gzipped NIfTI
if strcmpi(f_ext,'.gz') && length(f_name)>3 && ...
		strcmpi(f_name(end-3:end),'.nii')
	f_name = f_name(1:length(f_name)-4);
	f_ext = '.nii.gz';
end

if isempty(f_ext)
	if strcmpi(f_name,'fid')
		% Check if the file is a Bruker or Varian FID file
		if exist([f_path,'procpar'],'file') == 2
			dataformat = 'vnmr';
		else
			dataformat = 'bruker_raw';
		end
		return
	elseif strcmpi(f_name,'2dseq')
		dataformat = 'bruker_reco';
		return
	else
		% Check if the file is a DICOM file which can many times be without
		% file extension
		fid = fopen(filename,'r');
		if fid < 0
			return
		end
		
		% Seek over the possible DICOM preamble
		status = fseek(fid,128,-1);
		if status == -1
			% Unknown data format
			fclose(fid);
			return
		end
		
		% Try to read the 4 byte DICOM prefix
		[str,count] = fread(fid,4,'char');
		if count~=4
			fclose(fid);
			return
		end
		str = char(str.');
		if strcmp(str,'DICM')
			dataformat = 'dcm';
		end
		fclose(fid);
		return
	end
elseif strcmpi(f_ext,'.fid')
	dataformat = 'vnmr';
else
	if strcmpi(f_ext,'.xxm')
		dataformat='spect/ct';
	elseif any(strcmpi(f_ext,{'.nii','.nii.gz','.hdr','.img'}))
		dataformat = 'nifti';
		if strcmpi(f_ext,'.hdr')
			% Check if data is in Analyze/NIfTI or SPECT/CT format
			
			% Try to open the file for reading
			fid = fopen(filename,'r');
			if fid<0
				return
			end
			
			% Read 10 characters from the start
			[ident_str,count] = fread(fid,10,'char');
			if count~=10
				dataformat = '';
				fclose(fid);
				return
			end
			
			if strcmp(char(ident_str).','!INTERFILE')
				dataformat = 'spect/ct';
			else
				dataformat = 'nifti';
			end
		end
	elseif any(strcmpi(f_ext(2:end),...
			{'t1r','s1r','t2r','s2r','t1','t2','s1','s2','df','sf','r2','b1'}))
		% This file is probably Matlab MAT-file but could also be in the old
		% S.M.I.S. SUR-Format...
		fid = fopen(filename,'r');
		if fid < 0
			return
		end
		
		% The first 6 characters in MAT-File should be MATLAB
		[str,count]=fread(fid,6,'char');
		if count~=6
			fclose(fid);
			return
		end
		str = char(str).';
		if strcmpi(str,'MATLAB')
			dataformat = 'mat';
		else
			dataformat = 'sur';
		end
	elseif strcmpi(f_ext,'.sgl')
		dataformat = 'swift_sgl';
	else
		dataformat = lower(f_ext(2:end)); % Remove the dot
		
		% Check if file is a DICOM file
		fid = fopen(filename,'r');
		if fid < 0
			return
		end
		
		% Seek over the possible DICOM preamble
		status = fseek(fid,128,-1);
		if status == -1
			% Unknown data format
			fclose(fid);
			return
		end
		
		% Read the 4 byte DICOM prefix
		[str,count] = fread(fid,4,'char');
		fclose(fid);
		if count~=4
			return
		end
		str = char(str.');
		if strcmp(str,'DICM')
			dataformat = 'dcm';
		end
		
	end
	
	
end

