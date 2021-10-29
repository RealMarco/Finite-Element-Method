% Date & Time: 2020/06/26 09:04
% Project: The Constant strain triangular Element
% Step 3: Calculating the element Stiffness Matrix
% Instruction£ºThis function calculates the element stiffness matrix for each linear triangle with the Elasticity Matrix D, thickness t, 
%               and coordinates (xi, yi) for the first node, (xj, yj) for the second node, and (xm, ym) for the third node. 
%              It returns the 6¡Á6 element stiffness matrix k and the 3x6 Geometric Matrix B.

function [k,B]=CSTElementStiffness(D, t, xi, yi, xj, yj, xm, ym)
Area=0.5*(xi*(yj-ym)+xj*(ym-yi)+xm*(yi-yj));  % Calculating the trianguler element area
% Calculating the Element Geometric Matrix (aka Element Strain Matrix) B
bi = yj-ym;
bj = ym-yi; 
bm = yi-yj; 
ci = xm-xj;
cj = xi-xm; 
cm = xj-xi;
B=[bi 0  bj 0  bm 0;
   0  ci 0  cj 0  cm;
   ci bi cj bj cm bm]/(2*Area) ; 

% Calculating the element Stiffness Matrix k
k= t*Area* transpose(B) * D * B;

