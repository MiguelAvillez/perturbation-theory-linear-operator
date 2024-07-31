%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script uses the nomenclature, formulations and solutions from:
%   M. Avillez and D. Arnas, "Constructing Linear Operators Using Classical 
%   Perturbation Theory", TODO
% 
% Summary:
%   Computes the time derivative of the orbital elements subject to the J2
%   perturbation, with the differential equations listed in section IV.A.3
%   of the paper.
%
% Inputs:
%   state: [beta; x; y; p; raan; ctt; stt]
%       beta: sqrt(R/(sma * (1-ex^2-ey^2))), with sma the semi-major axis,
%           ex the x-eccentricity, and ey the y-eccentricity
%       x: ex/j2
%       y: ey/j2
%       p: cos(inc) / beta, with inc the inclination
%       raan: right ascension of the ascending node
%       ctt: cos(theta), with theta the argument of latitude
%       stt: sin(theta), with theta the argument of latitude
%   mu: gravitational parameter
%   j2: J2 coefficient of the gravity model
%   R: Radius of the central planet
%
% Outputs:
%   dstateDt: Time derivative of the state
%
%
% Authors: Miguel Avillez and David Arnas
% Modified: May 2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dstateDt = statePmJ2TimeDerivative(~, state, mu, j2, R)

beta = state(1);
x = state(2);
y = state(3);
p = state(4);
ctt = state(6);
stt = state(7);

dbetaDt = 3.*beta.^8.*ctt.*j2.*mu.^(1/2).*(1+(-1).*beta.^2.*p.^2).*R.^(-3/2).* ...
  stt.*(1+j2.*(ctt.*x+stt.*y)).^3;
dxDt = (3/2).*beta.^7.*mu.^(1/2).*R.^(-3/2).*stt.*(1+j2.*(ctt.*x+stt.*y)) ...
  .^3.*((-2).*beta.^2.*j2.*p.^2.*stt.*y+(-1).*ctt.*(1+(-1).*beta.^2.* ...
  p.^2).*(4.*ctt+3.*j2.*x+j2.*(ctt.^2+(-1).*stt.^2).*x+2.*ctt.*j2.* ...
  stt.*y)+((-1)+3.*(1+(-1).*beta.^2.*p.^2).*stt.^2).*(1+j2.*(ctt.*x+ ...
  stt.*y)));
dyDt = (-3/2).*beta.^7.*mu.^(1/2).*R.^(-3/2).*(1+j2.*(ctt.*x+stt.*y)).^3.*( ...
  (-2).*beta.^2.*j2.*p.^2.*stt.^2.*x+ctt.^2.*j2.*((-1)+5.*(1+(-1).* ...
  beta.^2.*p.^2).*stt.^2).*x+2.*ctt.^3.*j2.*(1+(-1).*beta.^2.*p.^2).* ...
  stt.*y+ctt.*((-1)+7.*(1+(-1).*beta.^2.*p.^2).*stt.^2).*(1+j2.*stt.* ...
  y));
dpDt = 0;
dooDt = (-3).*beta.^8.*j2.*mu.^(1/2).*p.*R.^(-3/2).*stt.^2.*(1+j2.*(ctt.*x+ ...
  stt.*y)).^3;
dcttDt = (-1).*beta.^3.*mu.^(1/2).*R.^(-3/2).*stt.*(1+j2.*(ctt.*x+stt.*y)) ...
  .^2.*(1+3.*beta.^6.*j2.*p.^2.*stt.^2.*(1+j2.*(ctt.*x+stt.*y)));
dsttDt = beta.^3.*ctt.*mu.^(1/2).*R.^(-3/2).*(1+j2.*(ctt.*x+stt.*y)).^2.*(1+ ...
  3.*beta.^6.*j2.*p.^2.*stt.^2.*(1+j2.*(ctt.*x+stt.*y)));

dstateDt = [dbetaDt; dxDt; dyDt; dpDt; dooDt; dcttDt; dsttDt];

end