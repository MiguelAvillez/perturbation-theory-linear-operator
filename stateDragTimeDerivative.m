%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script uses the nomenclature, formulations and solutions from:
%   M. Avillez and D. Arnas, "Constructing Linear Operators Using Classical 
%   Perturbation Theory", TODO
% 
% Summary:
%   Computes the time derivative of the orbital elements due to the drag 
%   perturbation (does not include the point-mass gravity term), using 
%   the differential equations listed in section V. of the paper.
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
%   cd: Drag coefficient of the spacecraft
%   S: cross-sectional area of the spacecraft
%   m: mass of the spacecraft
%   rho: atmospheric density
%
% Outputs:
%   dstateDt: Time derivative of the state
%
%
% Authors: Miguel Avillez and David Arnas
% Modified: May 2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dstateDt = stateDragTimeDerivative(~, state, mu, j2, R, cd, S, m, rho)

beta = state(1);
x = state(2);
y = state(3);
p = state(4);
ctt = state(6);
stt = state(7);

dbetaDt = (1/2).*beta.^2.*cd.*m.^(-1).*mu.^(1/2).*R.^(-1/2).*rho.*S.*(1+2.* ...
  ctt.*j2.*x+j2.^2.*x.^2+2.*j2.*stt.*y+j2.^2.*y.^2).^(1/2);
dxDt = (-1).*beta.*cd.*j2.^(-1).*m.^(-1).*mu.^(1/2).*R.^(-1/2).*rho.*S.*( ...
  ctt+j2.*x).*(1+2.*ctt.*j2.*x+j2.^2.*x.^2+2.*j2.*stt.*y+j2.^2.* ...
  y.^2).^(1/2);
dyDt = (-1).*beta.*cd.*j2.^(-1).*m.^(-1).*mu.^(1/2).*R.^(-1/2).*rho.*S.*( ...
  stt+j2.*y).*(1+2.*ctt.*j2.*x+j2.^2.*x.^2+2.*j2.*stt.*y+j2.^2.* ...
  y.^2).^(1/2);
dpDt = (-1/2).*beta.*cd.*m.^(-1).*mu.^(1/2).*p.*R.^(-1/2).*rho.*S.*(1+2.* ...
  ctt.*j2.*x+j2.^2.*x.^2+2.*j2.*stt.*y+j2.^2.*y.^2).^(1/2);
dooDt = 0;
dcttDt = 0;
dsttDt = 0;

dstateDt = [dbetaDt; dxDt; dyDt; dpDt; dooDt; dcttDt; dsttDt];

end