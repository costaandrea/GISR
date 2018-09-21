function KZ=getKz(alpha,beta,Cd,U,Kint,z,N)
%function KZ=getKz(alpha,beta,Cd,U,Kint,z,c0,c1)

% KZ=Kint+(1/6).*beta.^(-1).*Cd.*U.^2.*(c0+c1.*z).^(-1/2).*log(1+exp(1).^(2.* ...
%   alpha)).^(-1).*(1+tanh(alpha+(-1).*beta.^(-1).*U.^(-1).*z.*(c0+c1.*z).^( ...
%   1/2))); 
%%%Substitute c0+c1.*z with N
KZ=Kint+(1/6).*beta.^(-1).*Cd.*U.^2.*(N.^2).^(-1/2).*log(1+exp(1).^(2.* ...
  alpha)).^(-1).*(1+tanh(alpha+(-1).*beta.^(-1).*U.^(-1).*z.*(N.^2).^(1/2)));

