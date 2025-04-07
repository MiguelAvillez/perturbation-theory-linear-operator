%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script uses the nomenclature, formulations and solutions from:
%   M. Avillez and D. Arnas, "Constructing Linear Operators Using Classical 
%   Perturbation Theory", Journal of Guidance, Control, and Dynamics, 2025. 
%   https://doi.org/10.2514/1.G008683
% 
% Summary:
%   Computes the J2-perturbed frequency of the solution, using the
%   equations listed in section IV.C. of the paper
%
% Usage:
%   w = computeFrequency(state, mu, R, j2, expansionOrder)
%   [w, wArray] = computeFrequency(state, mu, R, j2, expansionOrder)
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
%   R: Radius of the central planet
%   j2: J2 coefficient of the gravity model
%   expansionOrder: order of the power expansion. Allowed values: 1, 2
%
% Outputs:
%   w: perturbed frequency
%   wArray: array with the values of the frequency for each order, e.g.
%   [w0, w1, w2].
%
%
% Authors: Miguel Avillez and David Arnas
% Modified: May 2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = computeFrequency(state, mu, R, j2, expansionOrder)

% Check order of the expansion
if expansionOrder ~= 2
    error("Invalid expansion order (%d).\n", expansionOrder);
end

beta = state(1);
p = state(4);
ctt = state(6);
stt = state(7);

% Compute frequency w
w0 = beta.^3.*mu.^(1/2).*R.^(-3/2);
w1 = (3/4).*beta.^7.*mu.^(1/2).*R.^(-3/2).*((-2)+ctt.^2.*(3+(-3).* ...
  beta.^2.*p.^2)+(-3).*stt.^2+beta.^2.*p.^2.*(8+3.* ...
  stt.^2));
wArray = [w0, w1];
if expansionOrder >= 2
    w2 = (3/32).*beta.^11.*mu.^(1/2).*R.^(-3/2).*(51+39.*ctt.^4.*((-1)+ ...
      beta.^2.*p.^2).^2+102.*stt.^2+39.*stt.^4+(-2).* ...
      beta.^2.*p.^2.*(98+276.*stt.^2+39.*stt.^4)+beta.^4.* ...
      p.^4.*(325+450.*stt.^2+39.*stt.^4)+(-6).*ctt.^2.*((-1) ...
      +beta.^2.*p.^2).*((-17)+(-39).*stt.^2+3.*beta.^2.* ...
      p.^2.*(25+13.*stt.^2)));
    wArray = [wArray, w2];
else
    w2 = 0;
end

w = w0 + j2 * w1 + j2^2 * w2;

% Create cell array for the appropriate number of outputs
nOutputs = nargout;
varargout = cell(1, nargout);
% Select outputs: either [w] or [w, wArray]
varargout{1} = w;
if nOutputs == 2
    varargout{2} = wArray;
elseif nOutputs > 2
    error("Invalid number of outputs.\n")
end

end