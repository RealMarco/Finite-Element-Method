% Date & Time: 2020/07/06 19:18
% Project: The Eight Points Quadrilateral Element
% Step 5: Calculating the element stress vectors & Calculating the element principal stresses
% Instruction£ºThis function calculates the element principal stresses using the element stress vector sigma. 
%              It returns a 3¡Á1 vector in the form [sigma1 sigma2 theta]T where sigma1 and sigma2 are the principal stresses for the element 
%              and theta is the principal angle.

function p=QuadraticQuadElementPStresses(sigma)
R = (sigma(1) + sigma(2))/2; 
Q = ((sigma(1) - sigma(2))/2)^2 + sigma(3)*sigma(3); 
M = 2*sigma(3)/(sigma(1) - sigma(2)); 
s1 = R +sqrt(Q); 
s2 = R -sqrt(Q); 
theta = (atan(M)/2)*180/pi; 
p = [s1 ; s2 ; theta];
