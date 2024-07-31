%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script uses the nomenclature, formulations and solutions from:
%   M. Avillez and D. Arnas, "Constructing Linear Operators Using Classical 
%   Perturbation Theory", TODO
% 
% Summary:
%   Given some basis functions, computes the value of the extended state,
%   according to the definition of the basis function specified through
%   mons. 
%
% Inputs:
%   mons: array representing a set of monomials. Each row
%       represents one monomial, with each coefficient being the exponent
%       of the associated element of the extended state. E.g.: If the
%       extended state is [x,y,z], a row [1,0,2] represents x*z^2.
%   extendedState: Value of the extended state
%
% Outputs:
%   basisFunctions: Value of the basis functions
%
%
% Authors: Miguel Avillez and David Arnas
% Modified: May 2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function extendedState = basisFunctions2extendedState(mons, basisFunctions)

% Number of dimensions (size of extended state)
nd = size(mons, 2);
% Initialize extended state 
extendedState = zeros(1,nd);

for i = 1:size(mons,1)
    nonzero_ids = find(mons(i,:));
    % Check if monomial is an extended state element or combination
    % If a combination, then continue
    if length(nonzero_ids) ~= 1 || mons(i,nonzero_ids(1)) ~= 1
        continue
    % If it is a state element save it to matrix
    else
        nonzero_id = nonzero_ids(1);
        extendedState(nonzero_id) = basisFunctions(i);
    end
end

end