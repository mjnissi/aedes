function [pc,p,cc] = pcorr(X)
% PCORR - Calculate partial correlation coefficients 
%   
%
% Synopsis: 
%        
%
% Description:
%        
% Examples:
%        
%
% See also:
%        

% This function is a part of Aedes - A graphical tool for analyzing 
% medical images
%
% Copyright (C) 2011 Juha-Pekka Niskanen <Juha-Pekka.Niskanen@uef.fi>
% 
% Department of Applied Physics, Department of Neurobiology
% University of Eastern Finland, Kuopio, FINLAND
%
% This program may be used under the terms of the GNU General Public
% License version 2.0 as published by the Free Software Foundation
% and appearing in the file LICENSE.TXT included in the packaging of
% this program.
%
% This program is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
% WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

t=perms(1:size(X,2));
t=t(1:2:end,:);
T = size(t,1);

Z_ind = t(:,1:end-2);
X_ind = t(:,end-1);
Y_ind = t(:,end);

dz = size(Z_ind,2);
n = size(X,1);

cc = diag(ones(1,size(X,2)));
pc = diag(ones(1,size(X,2)));
p = zeros(size(cc));

for ii=1:T
	x=X(:,X_ind(ii));
	y=X(:,Y_ind(ii));
	z=[ones(size(y,1),1) X(:,Z_ind(ii,:))];
	xx = x-z*(z\x);
	yy = y-z*(z\y);
	C = cov(xx,yy);
	C=C./(std(xx)*std(yy));
	coef = C(2);
	
	% Correlation coefficients
	C2 = cov(x,y);
	C2=C2./(std(x)*std(y));
	cc(X_ind(ii),Y_ind(ii)) = C2(2);
	cc(Y_ind(ii),X_ind(ii)) = C2(2);
	
	% Partial correlation coefficients
	pc(X_ind(ii),Y_ind(ii)) = coef;
	pc(Y_ind(ii),X_ind(ii)) = cc(X_ind(ii),Y_ind(ii));
	
	% P-values
	df = max(n - dz - 2,0); % degrees of freedom
	t = sign(coef) .* Inf;
	k = (abs(coef) < 1);
	t(k) = coef(k) ./ sqrt(1-coef(k).^2);
	t = sqrt(df).*t;
	
	pval1 = 2*aedes_tdist(-abs(t),df); % Two-tailed
	%pval2 = tdist(-t,df);       % greater than...
	%pval3 = tdist(t,df);        % lower than...
	
	p(X_ind(ii),Y_ind(ii)) = pval1;
	p(Y_ind(ii),X_ind(ii)) = p(X_ind(ii),Y_ind(ii));
end

