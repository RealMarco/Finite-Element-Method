 % Date & Time: 2020/06/20 11:00
% Project: The Plane Frame 
% Instruction: This is the MAIN fuction of the Plane Frame Element project.

% Functions used in the Plane Frame project (PFE - Plane Frame Element, aka, Beam Element with axial forces): 
% BeamLength(x1,y1,x2,y2)每This function returns the element length given the coordinates of the ?rst node(x1,y1) and the coordinates of the second node (x2, y2).
% PlaneFrameElementStiffness(E,A,J,L,x1,y1,x2,y2)每This function calculates the element stiffness matrix for each plane frame element with modulus of elasticity E, cross sectional area A, moment of inertia J, length L and coordinates x1,y1,x2,y2. 
%                                                  It returns the 6℅ 6 element stiffness matrix in GLOBAL coordinate system k, coordinate Transformation Matrix lamda,element stiffness matrix in LOCAL coordinate system k_
% PlaneFrameAssemble(k,nodes_num,beam) 每 This function assembles every element stiffness matrix k(:,:,i) into the global stiffness matrix K. It returns the 3n℅3n (n=number of nodes) global stiffness matrix K
% PlaneFrameElementNodalInternalForces(k_ ,lamda,U,beam)每This function calculates the nodal internal force vector of each PFE using the element stiffness Matrices in LOCAL coordinate system k_, 
%                                                        coordinate Transformation Matrices lamda,the displacement vector U, and the beam vector.It returns the 6℅number_of_beams matrix, and f(:,i) stands for the nodal internal force vector of beam i.   
% PlaneFrameElementInternalForceDiagram(beam,node,lamda,F)每This function plots the axial force diagram for the element 
%                                                           with object node, beam, coordinate Transformation Matrices lamda,and internal force vectors F .
clear;


%% ------Step1: Discretize the Frame manually & Input the details of each beam element-----
% Choose the inputting mode and Input data
disp('Hello, it is the Plane Frame Element (aka Beam Element with axial forces) project. Units: cm, N');
disp('------Step1: Discretize the Frame manually & Input the details of each beam element-----');
input_mode= input('Please input 0 to input data by the keyboard,\nor input 1 to input data and show the results of Example 1 (from the page 79 of the textbook Variational principle and Finite Element Method) automatically,\nor input 2 to input data and show the results of Example 2 (from the page 143 of Reference) automatically: \n');

jun=true;
while jun
    if input_mode == 0
        %% Input data by the keyboard.
        fprintf('\nInput data by the keyboard.\nPlease note that You should discretize the structure, number the elements and number the nodes manually.\n');
  
%         U=0.0011*ones(nodes_num*3,1);   % '0.0011' is the judgment criteria of unknown displacement
%         P=1.111*ones(nodes_num*3 ,1);   % '1.111' is the judgment criteria of unknown external forces
%         
%         fprintf('\nInput the coordinates, displacements and external forces of the node %i.\n',i)  % "%i" means output a integer

        % Input the coordinates of the nodes
        nodes= input('\nInput the Nodal Coordinate Matrix which is in the size of n℅2, while n is the total nodal number of the structure. \nEach row vector is in the form [x coordinate, y coordinate] :\n ');
        % Input the connectivity list of beams  
        beam= input('\nInput the b℅2 Connectivity List of the Beam Elements, while b is the total number of beams. \nEach row vector is in the form [serial number of the first node , serial number of the second node] :\n');
        
        % Input the displacements and external forces of the nodes
        % U is the nodal displacement vector
        % P is the nonal external force vector
        U= input('\nInput the Nodal Displacement Vector which is in the size of 3n℅1, while n is the total nodal number of the structure. \nThe horizontal displacement u(cm), vertical displacement v(cm) and rotation Phi(rad) at node i must be queried by U(3*i-2), U(3*i-1), U(3*i) respectively. \nUse 0.0011 to represent the unknown displacement:\n ');
        P= input('\nInput the Nodal Force Vector which is in the size of 3n℅1, while n is the total nodal number of the structure. \nThe horizontal force Px(N), vertical force Py(N) and moment M(N﹞cm) at node i must be queried by P(3*i-2), P(3*i-1), P(3*i) respectively. \nUse 1.111 to represent the unknown force:\n ');

        nodes_num= size(nodes,1);
        beams_num= size(beam,1);
        
        % Input the elastic constants
        E = input('\nInput Modulus of elasticity E (N/cm^2): ');
        A = input('Input Cross sectional area A (cm^2): ');
        J = input('Input Moment of inertia J (cm^4): ');
        disp('The length L(cm) of each beam may be diverse and will be calculated from the x,y coordinates automatically.');
    
        jun=false;
    elseif input_mode == 1
        %% Example 1 - from the page 79 of the textbook Variational principle and Finite Element Method
        nodes_num=4;  
        beams_num = 3;
        nodes = [0 0;
                0 100;
                100 100;
                100 0];
        beam = [1 2; 2 3; 3 4];

        E = 2e7; A=10; J=25;

        U=[0;0;0; 0.0011;0.0011;0.0011;0.0011;0.0011;0.0011; 0;0;0 ];
        P=[1.111;1.111;1.111; 10000;0;0;0;0;0; 1.111;1.111;1.111 ];
    
        jun= false;
    elseif input_mode == 2
        %% Example 2 - from the page 143 of Reference
        nodes_num= 4;
        beams_num=3;
        nodes=[0 0;
              0 300;
            400 300;
            400 0];
        beam=[1 2;2 3;3 4];
        
        E=2.1e7; A=200; J=5000;
        
        U=[0;0;0; 0.0011;0.0011;0.0011;0.0011;0.0011;0.0011; 0;0;0 ];
        P=[1.111;1.111;1.111; -20000;0;0;0;0;1200000; 1.111;1.111;1.111 ];
        
        jun=false;
    else
        disp('Input error, only 0,1 or 2 is viable. Please input again.');
        input_mode= input('Please input 0 to input data by the keyboard,\nor input 1 to input data and show the results of Example 1 (from the page 79 of the textbook Variational principle and Finite Element Method) automatically,\nor input 2 to input data and show the results of Example 2 (from the page 143 of Reference) automatically: \n');
    end
end


%% ------Step2: Calculate the element Stiffness Matrix------
L=zeros(beams_num,1);       % beam length vector
k=zeros(6,6,beams_num);     % element stiffness Matrices in GlOBAL coordinate system.
lamda=zeros(6,6,beams_num); % Coordinate Transformation Matrices
k_=zeros(6,6,beams_num);    % element stiffness Matrices in LOCAL coordinate system.
for i = 1:beams_num
    x1= nodes(beam(i,1),1);
    y1= nodes(beam(i,1),2);  
    x2= nodes(beam(i,2),1);
    y2= nodes(beam(i,2),2);
    L(i)= BeamLength(x1,y1,x2,y2);  % Calculate the beam length 
    
    [k(:,:,i), lamda(:,:,i),k_(:,:,i)] =PlaneFrameElementStiffness(E,A,J,L(i),x1,y1,x2,y2);
end

%% ------Step3: Assemble element stiffness Matrices into the Frame(global) Stiffness Matrix ------
K = PlaneFrameAssemble(k,nodes_num,beam); 

%% ------Step4: Constraint handling (Import  boundary conditions) and Solve the system of linear equations------ 

% Calculate unknown Nodal Displacements with Known external forces 
tag=find(P~=1.111);   %(p.s.'1.111' is the judgment criteria) find the indices of known external forces      
K_partition= K(tag,tag);  %partition the frame stiffness Matrix K.

P_known=P(tag);   % find the known external forces
U_unknown=K_partition\P_known;  %the backslash operator ※\§ is used for Gaussian elimination
% U_unknown=inv(K_partition)*P_known
U(tag)=U_unknown;      % back to the U vector

% Calculate unknown reactions with Nodal Displacements 
tag2=find(P==1.111);  %find the indices of unknown external forces
K_partition2= K(tag2,:);  %
P_unknown= K_partition2 * U;  % Calculate unknown reactions
P(tag2)=P_unknown;   % back to the P vector 

%% ------Step5: Calculate nodal internal forces of each beam------
F=PlaneFrameElementNodalInternalForces(k_ ,lamda,U,beam);  % F is a 6* number of beams matrix, and F(:,i) stands for the 
                                                           % nodal internal force vector of beam i.                                                           

%% ------Step6: Plot the internal force diagrams------
% PlaneFrameElementInternalForceDiagram(beam,node,lamda,F) % version 1
PlaneFrameElementInternalForceDiagram2(beam,nodes,lamda,F)  % version 2






