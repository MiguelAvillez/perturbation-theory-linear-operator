%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script uses the nomenclature, formulations and solutions from:
%   M. Avillez and D. Arnas, "Constructing Linear Operators Using Classical 
%   Perturbation Theory", TODO
% 
% Summary:
%   Propagates the approximate linear system generated based on the
%   Lindstedt-Poincaré approximation, subject to the J2 and drag (order of
%   magnitude of J2^2, constant density) perturbations, as described in 
%   section V.A. of the paper.
%
%   1) The initial Keplerian elements defining the orbit to be
%       approximated are specified. 
%   2) The frequency of the solution and initial extended state are computed. 
%   3) The matrix defining the approximated system and the array specifying the 
%       monomials associated with the matrix are loaded from an external file.
%   4) The approximate linear system and the full nonlinear system are propagated, 
%       and the error of the former with respect to the latter is computed
%       and plotted.
%
% Dependencies:
%   basisFunctions2extendedState.m
%   computeFrequency.m
%   extendedState2basisFunctions.m
%   statePmJ2TimeDerivative.m
%   stateDragTimeDerivative.m
%
% Inputs:
%   Constants:
%       R: Radius of the central planet
%       mu: gravitational parameter
%       j2: J2 coefficient of the gravity model
%   Spacecraft and atmosphere constants:
%       cd: drag coefficient
%       S: cross-sectional area
%       m: mass
%       rho: atmospheric density
%   Initial Keplerian orbital elements:
%       sma: semi-major axis
%       ex: x-eccentricity
%       ey: y-eccentricity
%       inc: inclination
%       raan: right ascension of ascending node
%       tt: argument of latitude
%   Script variables:
%       numRevolutions: number of revolutions to compute
%       numPointsPerRevolution: number of points per computed revolution
%
% Outputs:
%   Plot of orbital elements error
%   M: matrix describing the system linearly
%   mons: array representing the monomials associated with M. Each row
%       represents one monomial, with each coefficient being the exponent
%       of the associated element of the extended state. E.g.: If the
%       extended state is [x,y,z], a row [1,0,2] represents x*z^2.
%   w: perturbed frequency associated with the Lindstedt-Poincaré method
%
%
% Authors: Miguel Avillez and David Arnas
% Modified: May 2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

initialTime = 0;

% Define constants
R = 6378.15e3; % central planet radius, m
mu = 3.986e14; % gravitational parameter, m^3/s^2
j2 = 1.08262668e-3; % J2 coefficient of the gravity model, -
deg2rad = pi/180;
rad2deg = 1/deg2rad;

% Sun-synchronous frozen orbit
sma = 7077.722e3; % semi-major axis [m]
ex = 4.5742e-4; % x-eccentricity [-]
ey = 0; % y-eccentricity [-]
inc = 98.186 * deg2rad; % inclination [rad]
raan = 42 * deg2rad; % right ascension of ascending node [rad]
tt = 0 * deg2rad; % argument of latitude [rad]

% Satellite
cd = 2.2; % drag coefficient, -
S = 8.5; % cross-sectional area, m^2
m = 1285; % mass, kg

% Constant atmospheric density
rho = 2e-14; % kg/m^3

% Number of revolutions
numRevolutions = 20;
% Number of points per revolution used for plotting
numPointsPerRevolution = 100;

% Order of the expansion. Allowed values: 2
expansionOrder = 2;

%% Compute initial state with proposed orbital elements

beta = sqrt(R/(sma * (1-ex^2-ey^2)));
x = ex/j2;
y = ey/j2;
p = cos(inc) / beta;
ctt = cos(tt);
stt = sin(tt);

% Initial state with orbital elements
stateInitial = [beta; x; y; p; raan; ctt; stt];

%% Compute auxiliary quantities

% Check order of the expansion
if expansionOrder ~= 2
    error("Invalid expansion order (%d).\n", expansionOrder);
end

% Compute frequency w
[w, wArray] = computeFrequency(stateInitial, mu, R, j2, expansionOrder);
% Extract value of the frequency at each order of J2
w0 = wArray(1);
w1 = wArray(2);
w2 = wArray(3);

% Perturbed orbital period (to define propagation time)
orbitalPeriod = 2*pi/w; 

%  Construct extended initial state
extendedStateInitial = [
    % Initial state of order 0
    beta; p; raan; ctt; stt; 
    % Initial state of order 1
    0; x; y; 0; 0; 0; 
    % Initial state of order 2
    zeros(7,1); 
    % Define constant auxiliary variables k_beta, k_{\omega_1},
    % k_{\omega_2}, and k_\rho
    1/beta; w1/w0; w2/w0; cd./j2.^2.*m.^(-1).*R.*rho.*S]; 

%% Load operator matrix "M", and matrix with monomials "mons"
M = load(sprintf("Data/j2DragLP_m_order%d.txt", expansionOrder));
mons = load(sprintf("Data/j2DragLP_mons_order%d.txt", expansionOrder));

%% Propagate full dynamics and linear system

basisFunctionsInitial = extendedState2basisFunctions(mons, extendedStateInitial);

% Propagation settings
numPoints = numPointsPerRevolution * numRevolutions;
timeHistory = linspace(initialTime, numRevolutions * orbitalPeriod, numPoints);

tolerance = 1e-13;
odeOptions = odeset('RelTol', tolerance,'AbsTol', tolerance);

% Propagate full dynamics 
derivativeFunction = @(t,x) statePmJ2TimeDerivative(t, x, mu, j2, R) + stateDragTimeDerivative(t, x, mu, j2, R, cd, S, m, rho);
[~, propagatedStateHistory] = ode89(derivativeFunction, timeHistory, stateInitial, odeOptions);

% Propagate approximate linear system 
[~, basisFunctionsHistory] = ode89(@(t,x) M * x * w, timeHistory, basisFunctionsInitial, odeOptions);
% Compute evolution of the state
matrixStateHistory = zeros(length(timeHistory), length(stateInitial));
for i = 1:length(timeHistory)
    extendedState = basisFunctions2extendedState(mons, basisFunctionsHistory(i,:));

    matrixStateHistory(i,1) = extendedState(1) + j2 * extendedState(6); % beta
    matrixStateHistory(i,2) = extendedState(7); % x
    matrixStateHistory(i,3) = extendedState(8); % y
    matrixStateHistory(i,4) = extendedState(2); % p
    matrixStateHistory(i,5) = extendedState(3) + j2 * extendedState(9); % raan
    matrixStateHistory(i,6) = extendedState(4) + j2 * extendedState(10); % ctt
    matrixStateHistory(i,7) = extendedState(5) + j2 * extendedState(11); % stt

    if expansionOrder == 2
        matrixStateHistory(i,1) = matrixStateHistory(i,1) + j2^2 * extendedState(12); % beta
        matrixStateHistory(i,2) = matrixStateHistory(i,2) + j2 * extendedState(13); % x
        matrixStateHistory(i,3) = matrixStateHistory(i,3) + j2 * extendedState(14); % y
        matrixStateHistory(i,4) = matrixStateHistory(i,4) + j2^2 * extendedState(15); % p
        matrixStateHistory(i,5) = matrixStateHistory(i,5) + j2^2 * extendedState(16); % raan
        matrixStateHistory(i,6) = matrixStateHistory(i,6) + j2^2 * extendedState(17); % ctt
        matrixStateHistory(i,7) = matrixStateHistory(i,7) + j2^2 * extendedState(18); % stt
    end
end

%% Plot orbital elements error

yLabels = ["\beta", "X", "Y", "p", "\Omega", "\theta"];
unitLabels = ["[-]", "[-]", "[-]", "[-]", "[rad]", "[rad]"];
fontSize = 14;

figure(1)
set(gcf,'Position', [100 100 800 425])
for var = 1:length(yLabels)

    % Compute error in first five orbital elements (beta, x, y, p, raan)
    if var < 6
        error = matrixStateHistory(:,var) - propagatedStateHistory(:,var);
    % Compute error in argument of latitude (instead of its sine and cosine)
    else
        matrixTt = atan2(matrixStateHistory(:,7), matrixStateHistory(:,6));
        propagatedTt = atan2(propagatedStateHistory(:,7), propagatedStateHistory(:,6));
        error = matrixTt - propagatedTt;
    end

    subplot(2,3,var)
    plot(timeHistory/orbitalPeriod, error, "LineWidth", 2)

    grid;
    xlabel('Time [ T ]')
    ylabel("Error " + yLabels(var) + " " + unitLabels(var))
    set(gca, 'FontSize', fontSize);

end