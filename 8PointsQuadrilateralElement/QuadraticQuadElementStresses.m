% Date & Time: 2020/07/06 18:52
% Project: The Eight Points Quadrilateral Element
% Step 5: Calculating the element stress vectors & Calculating the element principal stresses
% Instruction£ºThis function calculates the element stresses using the Elasticity Matrix D, the Geometric Matrix B(:,:,i) 
%              and the element displacement vector u. It returns two 3*1 results (in the form of [sigma_x sigma_y tao_xy]T)¨C the general quadratic stress functions in ¦Î and ¦Ç, 
%              and the numerical values of the stresses at the centroid of the element.

function [sigma_sym, centroid_sigma] =QuadraticQuadElementStresses(D,B,u)
syms ksi eta
% symbolic stress
sigma_sym = D*B*u;
% numerical stress at the centroid of the element
centroid_sigma = subs(sigma_sym, {ksi, eta},{0,0});
centroid_sigma= double(centroid_sigma );
