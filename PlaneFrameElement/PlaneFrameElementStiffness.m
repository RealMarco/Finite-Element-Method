% Date & Time: 2020/06/20 20:09
% Project: The Plane Frame 
% Step2 Calculate the element Stiffness Matrix
% Instruction£ºThis function calculates the ELEMENT stiffness matrix for each plane frame element with 
%              modulus of elasticity E, cross sectional area A, moment of inertia J, length L and coordinates x1,y1,x2,y2. 
%              It returns the 6¡Á 6 element stiffness matrix in GLOBAL coordinate system k, 
%              coordinate Transformation Matrix lamda, element stiffness matrix in LOCAL coordinate system k_


function [ k,lamda,k_ ] = PlaneFrameElementStiffness(E,A,J,L,x1,y1,x2,y2)
x=x2-x1 ; y=y2-y1;  
C = (1*x + 0*y)/sqrt(x^2+y^2); % C=cos(theta)
S = (1*y - 0*x)/sqrt(x^2+y^2); % S=sin(theta)
lamda = [C   S   0   0   0   0;
        -S   C   0   0   0   0;
         0   0   1   0   0   0;
         0   0   0   C   S   0;
         0   0   0  -S   C   0;
         0   0   0   0   0   1];   % lamda is the Coordinate Transformation Matrix

c0 = E*A/L; 
c1 = E*J/L;
c2 = E*J/(L*L);
c3 = E*J/(L*L*L);
k_ = [c0        0       0       -c0     0       0;
       0        12*c3   6*c2     0    -12*c3    6*c2;
       0        6*c2    4*c1     0    -6*c2     2*c1;
     -c0        0       0        c0     0       0;
       0       -12*c3  -6*c2     0     12*c3   -6*c2;
       0        6*c2    2*c1     0    -6*c2     4*c1]; % Calculate stiffness matrix of each element in LOCAL coordinate system.

lamda_T =transpose(lamda);
k=lamda_T * k_ * lamda; % Calculate stiffness matrix of each element in GLOBAL coordinate system














