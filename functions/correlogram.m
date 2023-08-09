function  corr_coeffs = correlogram( X, Y, order, is_unbiased )
% Copyright (C) 2005 Daniele de Rigo
% 
% 
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 2 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program; if not, write to the Free Software
%     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
% 
% 
% corr_coeffs = correlogram( X, Y, order, is_unbiased )
% 
% Given two signals  X  and  Y  (with the same length N),
% compute the series of the correlation coefficient
% between the two sub-sequences  Xi = X(1:(N-i))
% and  Yi = Y((1+i):N)  looping i from 0 to min(order,N-2)
% 
% the correlation coefficient is computed as:
% 
%        E[ (Xi - E[Xi]) (Yi - E[Yi]) ]
%  _____________________________________________
%      _______________________________________
%   \ | 
%    \|  E[(Xi - E[Xi])^2]   E[(Yi - E[Yi])^2]
% 
% If  is_unbiased  is false, then all the   E[.] 
% operators will be computed as   mean(.)  .
% Else the inner   E[.]   operators will always
% be computed as   mean(.)  , but the extern ones
% (i.e. those whose argument contains another E[.])
% will consider the freedom degrees reduction by
% computing itself as   sum(.)/(length(.)-1).
% 
% 
%  X            first  signal vector
%  Y            second signal vector (must have the same length
%               as  X  : if Y == X, then it will be computed
%               the autocorrelogram) [default: X]
%  order        max correlation coefficient order (must be
%               non negative integer) [default: length(X)-2] 
%  is_unbiased  boolean: if true compute an unbiased correlation
%               coefficient [default: 0]

if(nargin<1) 
  disp(  'too few arguments'  )
  disp(  'usage: correlogram( X, Y, order, is_unbiased )'  )
  return;
end;

if(nargin>4)
  disp(  'too many arguments'  )
  disp(  'usage: correlogram( X, Y, order, is_unbiased )'  )
  return;
end;

if(min(size(X))~=1)
  disp(  'error: the first argument must be a vector'  )
  return;
end;
N=length(X);

if(N<2)
  disp(  'error: the first argument must be a vector with at least 2 elements'  )
  return;
end;

if(nargin<2)
  Y=X;
end;

if( any(size(Y)~=size(X)) )
  disp([  'error: the second argument (size:['  num2str(size(Y,1))  ','  num2str(size(Y,2))  ']) must have the same size as the first one (size:['  num2str(size(X,1))  ','  num2str(size(X,2))  '])'  ])
  return;
end;

X=reshape(X,N,1); Y=reshape(Y,N,1);

if(nargin<3)
  order=N-2;
else
  %if 0 end; %-------------------------...etc...
  order=min(order,N-2);
end;

if(nargin<4)
  is_unbiased=0;
end;

corr_coeffs=zeros(1,order+1);
X2=X.^2; Y2=Y.^2;
if(is_unbiased) 
  for k=0:order
    iX=1:(N-k); iY=(1+k):N; Nk=N-k-1; N2k=(N-k)/Nk;
    mXY=X(iX)'*Y(iY)/Nk; 
    mX =mean( X(iX)); mY =mean( Y(iY)); 
    mX2=mean(X2(iX)); mY2=mean(Y2(iY)); 
    corr_coeffs(k+1)=(mXY-N2k*mX*mY)/sqrt((mX2-N2k*mX^2)*(mY2-N2k*mY^2)); 
  end;
else
  for k=0:order
    iX=1:(N-k); iY=(1+k):N; Nk=N-k;
    mXY=X(iX)'*Y(iY)/Nk;
    mX =mean( X(iX)); mY =mean( Y(iY)); 
    mX2=mean(X2(iX)); mY2=mean(Y2(iY)); 
    corr_coeffs(k+1)=(mXY-mX*mY)/sqrt((mX2-mX^2)*(mY2-mY^2)); 
  end;
end;


if(nargout==0)
  n=N:-1:(N-order);
  std_conf_95=n.^(-.5)*1.96;    % symmetrical confidence: 0.95 
  std_conf_995=n.^(-.5)*2.8070; % symmetrical confidence: 0.995
  lag=0:order;
  plot(lag,corr_coeffs,'b',lag,std_conf_95,'m',lag,std_conf_995,'r',lag,-std_conf_95,'m',lag,-std_conf_995,'r');
  title(   'Correlogram'  );
  xlabel(  'lag (delay steps)'  );
  ylabel(  'correlation'  );
  legend(  'correlation'  ,  '95 % confidence'  ,  '99.5% confidence'  );
  clear corr_coeffs;
end;


