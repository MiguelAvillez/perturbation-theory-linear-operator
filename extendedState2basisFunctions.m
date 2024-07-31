%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script uses the nomenclature, formulations and solutions from:
%   M. Avillez and D. Arnas, "Constructing Linear Operators Using Classical 
%   Perturbation Theory", TODO
% 
% Summary:
%   Given an extended state, computes the value of the basis functions,
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
function basisFunctions = extendedState2basisFunctions(mons, extendedState)

if length(extendedState) ~= size(mons,2)
    error("Size of initial state (%d) and number of dimensions (%d) is inconsistent.\n", ...
        length(extendedState), size(mons,2));
end

basisFunctions = zeros(size(mons,1),1);

% Loop over basis functions
for i = 1:size(mons,1)
    if nnz(mons(i,:)) == 0
        ic = 0;
    else
        ic = 1;
        % Loop over extended state elements
        for j = 1:size(mons,2)
            if mons(i,j) ~= 0
                ic = ic * extendedState(j)^mons(i,j);
            end
        end
    end
    basisFunctions(i) = ic;
end

end