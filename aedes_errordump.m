function aedes_errordump(errorstruct)
% AEDES_ERRORDUMP - Write the lasterror structure into file
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


debug=false;

if debug
  rethrow(errorstruct)
else
  % Construct file name and path for the errordump file
  error_dir = [prefdir(1),filesep];
  error_fname =['aedes_errordump_',datestr(now,30),'.mat'];
  
  % Get revision information
  [rev,repo]=aedes_revision;
  
  % Save the errorstruct
  save([error_dir,error_fname],'errorstruct','rev','repo','-mat')
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Error dialog
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  hh=errordlg({['Oops... An unhandeled error has occurred in Aedes! ',...
	'This is probably a bug. ',...
	'An error dump was written into the following file:'],...
	['"',error_dir,error_fname,'"'],'','',...
	'The returned error message was:',...
	errorstruct.message,errorstruct.identifier,'',...
	['If you get this error dialog every time you try to close Aedes ' ...
	'main window, run the command "aedes_killfigs" in the Matlab ' ...
	'workspace. This will delete ALL currently open Matlab figures ' ...
	'with brute force!!!']},...
	'UNHANDELED ERROR!','modal');
end
