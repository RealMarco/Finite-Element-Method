% Date & Time: 2020/07/06 10:43
% Project: The Eight Points Quadrilateral Element
% Step 2: Calculating the element Stiffness Matrix
% Instruction£ºThis function calculates the element stiffness matrix for each linear triangle with the Elasticity Matrix D, thickness t, 
%              and coordinates (x1, y1),(x2, y2),(x3, y3),(x4,y4) for the 1st, 2nd, 3rd, 4th element corner node respectively.  
%              It returns the 16¡Á16 element stiffness matrix k and the 3x16 Geometric Matrix B.
function [k,B]= QuadraticQuadElementStiffness(D, thickness, x1, y1, x2, y2, x3, y3, x4, y4)

syms ksi eta; 
x5 = (x1 + x2)/2;  y5 = (y1 + y2)/2; 
x6 = (x2 + x3)/2;  y6 = (y2 + y3)/2;
x7 = (x3 + x4)/2;  y7 = (y3 + y4)/2; 
x8 = (x4 + x1)/2;  y8 = (y4 + y1)/2; 

N1 = (1-ksi)*(1-eta)*(-ksi-eta-1)/4; 
N2 = (1+ksi)*(1-eta)*(ksi-eta-1)/4; 
N3 = (1+ksi)*(1+eta)*(ksi+eta-1)/4; 
N4 = (1-ksi)*(1+eta)*(-ksi+eta-1)/4; 
N5 = (1-eta)*(1+ksi)*(1-ksi)/2; 
N6 = (1+ksi)*(1+eta)*(1-eta)/2; 
N7 = (1+eta)*(1+ksi)*(1-ksi)/2; 
N8 = (1-ksi)*(1+eta)*(1-eta)/2;

x = N1*x1 + N2*x2 + N3*x3 + N4*x4 + N5*x5 + N6*x6 + N7*x7 + N8*x8; 
y = N1*y1 + N2*y2 + N3*y3 + N4*y4 + N5*y5 + N6*y6 + N7*y7 + N8*y8; 

xk = diff(x,ksi); xe = diff(x,eta); 
yk = diff(y,ksi); ye = diff(y,eta); 
J_old = xk*ye - yk*xe; 
J = simplify(J_old); 

N1k = diff(N1,ksi); N1e = diff(N1,eta); 
N2k = diff(N2,ksi); N2e = diff(N2,eta); 
N3k = diff(N3,ksi); N3e = diff(N3,eta); 
N4k = diff(N4,ksi); N4e = diff(N4,eta); 
N5k = diff(N5,ksi); N5e = diff(N5,eta);
N6k = diff(N6,ksi); N6e = diff(N6,eta); 
N7k = diff(N7,ksi); N7e = diff(N7,eta);
N8k = diff(N8,ksi); N8e = diff(N8,eta); 

B11 = ye*N1k - yk*N1e; 
B12 = 0; 
B13 = ye*N2k - yk*N2e; 
B14 = 0; 
B15 = ye*N3k - yk*N3e; 
B16 = 0; 
B17 = ye*N4k - yk*N4e; 
B18 = 0; 
B19 = ye*N5k - yk*N5e; 
B110 = 0; 
B111 = ye*N6k - yk*N6e; 
B112 = 0; 
B113 = ye*N7k - yk*N7e;
B114 = 0; 
B115 = ye*N8k - yk*N8e; 
B116 = 0; 
B21 = 0; 
B22 = xk*N1e - xe*N1k; 
B23 = 0; 
B24 = xk*N2e - xe*N2k; 
B25 = 0; 
B26 = xk*N3e - xe*N3k; 
B27 = 0; 
B28 = xk*N4e - xe*N4k; 
B29 = 0; 
B210 = xk*N5e - xe*N5k; 
B211 = 0; 
B212 = xk*N6e - xe*N6k; 
B213 = 0; 
B214 = xk*N7e - xe*N7k; 
B215 = 0; 
B216 = xk*N8e - xe*N8k; 
B31 = xk*N1e - xe*N1k; 
B32 = ye*N1k - yk*N1e; 
B33 = xk*N2e - xe*N2k; 
B34 = ye*N2k - yk*N2e; 
B35 = xk*N3e - xe*N3k; 
B36 = ye*N3k - yk*N3e; 
B37 = xk*N4e - xe*N4k; 
B38 = ye*N4k - yk*N4e; 
B39 = xk*N5e - xe*N5k; 
B310 = ye*N5k - yk*N5e; 
B311 = xk*N6e - xe*N6k; 
B312 = ye*N6k - yk*N6e; 
B313 = xk*N7e - xe*N7k; 
B314 = ye*N7k - yk*N7e; 
B315 = xk*N8e - xe*N8k; 
B316 = ye*N8k - yk*N8e; 

B_old= [B11 B12 B13 B14 B15 B16 B17 B18 B19 B110 B111 B112 B113 B114 B115 B116; 
     B21 B22 B23 B24 B25 B26 B27 B28 B29 B210 B211 B212 B213 B214 B215 B216; 
     B31 B32 B33 B34 B35 B36 B37 B38 B39 B310 B311 B312 B313 B314 B315 B316]/J;

B = simplify(B_old); 

BD = transpose(B)*D*B*J; 
z = thickness* int(int(BD, eta, -1, 1), ksi, -1, 1); 
k = double(z);






