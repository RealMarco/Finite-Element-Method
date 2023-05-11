% Date & Time: 2020/06/26 12:05
% Project: The Constant strain triangular Element
% Step 6:  Calculating the element stress vectors & Calculating the element principal stressess
% Instruction£º This function calculates the element stresses using the Elasticity Matrix D, the Geometric Matrix B(:,:,i) 
%               and the element displacement vector u. It returns the 3x1 element stress vector in the form [sigma_x sigma_y tao_xy]T
                            
function sigma=CSTElementStresses(D,B,u)
sigma= D * B * u;
