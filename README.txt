This repository contains the source-code describing the formulations and solutions in
	M. Avillez and D. Arnas, "Constructing Linear Operators Using Classical 
    Perturbation Theory", Journal of Guidance, Control, and Dynamics, 2025. 
	https://doi.org/10.2514/1.G008683

Each of the main files (listed below) imports the linear operator "M" and the array "mons"
(describing the basis functions associated with "M") from data files. These are then used
to propagate a frozen sun-synchronous orbit. Finally, the error of the approximate linear
solution is computed and plotted.

-----------------------------------------------------------------------------------
Main files:
- main_j2_simplePowerExpansion.m
	- Reproduces the solution from section IV.B. of the paper
	- Perturbations: J2
	- Solution obtained using a simple power expansion (without frequency control)
	- Solution provided for order 1 and 2. Associated data files:
		- j2SPE_m_order1.txt
		- j2SPE_m_order2.txt
		- j2SPE_mons_order1.txt
		- j2SPE_mons_order2.txt

- main_j2_lindstedtPoincare.m
	- Reproduces the solution from section IV.C. of the paper
	- Perturbations: J2
	- Solution obtained using the Lindstedt-Poincare method
	- Solution provided for order 1 and 2. Associated data files:
		- j2LP_m_order1.txt
		- j2LP_m_order2.txt
		- j2LP_mons_order1.txt
		- j2LP_mons_order2.txt

- main_j2Drag_lindstedtPoincare.m
	- Reproduces the solution from section V.A. of the paper
	- Perturbations: J2, drag (constant density, magnitude on the order of J2^2)
	- Solution obtained using the Lindstedt-Poincare method (i.e. with constant frequency)
	- Solution provided for order 2. Associated data files:
		- j2DragLP_m_order2.txt
		- j2DragLP_mons_order2.txt

- main_j2Drag_modifiedLindstedtPoincare.m
	- Reproduces the solution from section V.B. of the paper
	- Perturbations: J2, drag (constant density, magnitude on the order of J2^2)
	- Solution obtained using the modified Lindstedt-Poincare method (i.e. with varying frequency)
	- Solution provided for order 2. Associated data files:
		- j2DragMLP_m_order2.txt
		- j2DragMLP_mons_order2.txt

-----------------------------------------------------------------------------------
Full list of files:
- basisFunctions2extendedState.m
- computeFrequency.m
- extendedState2basisFunctions.m
- main_j2_lindstedtPoincare.m
- main_j2_simplePowerExpansion.m
- main_j2Drag_lindstedtPoincare.m
- main_j2Drag_modifiedLindstedtPoincare.m
- stateDragTimeDerivative.m
- statePmJ2TimeDerivative.m
- Data
	- j2DragLP_m_order2.txt
	- j2DragLP_mons_order2.txt
	- j2DragMLP_m_order2.txt
	- j2DragMLP_mons_order2.txt
	- j2LP_m_order1.txt
	- j2LP_m_order2.txt
	- j2LP_mons_order1.txt
	- j2LP_mons_order2.txt
	- j2SPE_m_order1.txt
	- j2SPE_m_order2.txt
	- j2SPE_mons_order1.txt
	- j2SPE_mons_order2.txt