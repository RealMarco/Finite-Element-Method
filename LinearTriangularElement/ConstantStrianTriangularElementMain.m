% Date & Time: 2020/06/20 11:00
% Project: The Constant strain triangular Element
% Instruction: This is the MAIN fuction of the Constant strain triangular Element project.

% Functions used in the Constant Strian Triangular Element (Abbreviated as CST Element, aka the Linear Triangular Element) project:
% CSTElementArea(xi,yi,xj,yj,xm,ym)每This function returns the element area given the coordinates of the ?rst node(xi,yi),
%                                    the coordinates of the second node (xj, yj), and the coordinates of the third node (xm, ym).
% CSTElementStiffness(D, t, xi, yi, xj, yj, xm, ym) 每 This function calculates the element stiffness matrix for each linear triangle with the Elasticity Matrix D, 
%                                    thickness t, and coordinates (xi, yi) for the first node, (xj, yj) for the second node, and (xm, ym) for the third node. 
%                                    It returns the 6℅6 element stiffness matrix k and the 6x6 Geometric Matrix B.
% CSTElementStiffnessMatrixAssemble(k,elements_num, elements)每This function assembles the element stiffness matrices k into the global stiffness matrix K. 
%                                                              It returns the 2n℅2n (n represents the number of nodes) global stiffness matrix K.
% CSTElementStresses(D,B,u) 每 This function calculates the element stresses using the Elasticity Matrix D, the Geometric Matrix B(:,:,i) 
%                              and the element displacement vector u. It returns the 3x1 element stress vector in the form [sigma_x sigma_y tao_xy]T.
                            
% CSTElementPrincipalStresses(sigma)每This function calculates the element principal stresses using the element stress vector sigma. 
%                                     It returns a 3℅1 vector in the form [sigma1 sigma2 theta]T where sigma1 and sigma2 are the principal stresses for the element 
%                                     and theta is the principal angle.

clear;
tic

%% ------Step 1: Discretizing the structure (Input the details of each element)------
% Choose the inputting mode and Input data
disp('Hello, it is The Constant Strian Triangular Element ( aka Linear Triangular Element ) project. Units: m, kN');
disp('------Step1: Discretize the Frame manually & Input the details of each element-----');
input_mode= input('Please input 0 to input data by the keyboard,\nor input 1 to input data and show the results of Example 1 (from the page 109 of the textbook Variational Principle and Finite Element Method) automatically,\nor input 2 to input data and show the results of Example 2 (from the page 223 of Reference) automatically: \n');

jun=true;
while jun
    if input_mode == 0
        %% Input data by the keyboard.
        fprintf('\nInput data by the keyboard.\nPlease note that You should discretize the structure, number the elements and number the nodes manually.\n');
 
        % Input the coordinates of the nodes
        nodes= input('\nInput the Nodal Coordinate Matrix which is in the size of n℅2, while n is the total nodal number of the structure. \nEach row vector is in the form [x coordinate, y coordinate] :\n ');
        % Input the connectivity list of elements  
        elements= input('\nInput the b℅3 Connectivity List of the Beam Elements, while b is the total number of elements. \nEach row vector is in the form [serial number of node i, serial number of node j, serial number of node m] :\n');
        
        % Input the displacements and external forces of the nodes
        % U is the nodal displacement vector
        % P is the nonal external force vector
        U= input('\nInput the Nodal Displacement Vector which is in the size of 2n℅1, while n is the total nodal number of the structure. \nThe horizontal displacement u(m) and vertical displacement v(m) at node i must be queried by U(2*i-1), U(2*i) respectively. \nUse 0.0011 to represent the unknown displacement:\n ');
        P= input('\nInput the Nodal Force Vector which is in the size of 2n℅1, while n is the total nodal number of the structure. \nThe horizontal force Px(kN) and vertical force Py(kN) at node i must be queried by P(2*i-1), P(2*i) respectively. \nUse 1.111 to represent the unknown force:\n ');

        nodes_num= size(nodes,1);
        elements_num= size(elements,1);
        
        % Input the elastic constants
        E = input('\nInput Modulus of elasticity E (kN/m^2): ');
        mu = input('Input Poisson＊s ratio米: ');
        t = input('Input Thickness t (m): ');
        disp('The Area A(m^2) of each element will be calculated from the x,y coordinates automatically.');
        
        %% Choosing the Elasticity Matrix D according to the Question Type (plane stress or plane strain) 
        ques_type = input('\nQuestion Type - Please input 1 (stands for the plane stress question) or 0 (stands for the plane strain question): ');
        D=zeros(3,3);   %initialize
        ju=true;
        while ju
            if ques_type == 1   % the plane stress case
                D= (E/(1-mu*mu)) * [1 mu 0;
                                   mu 1  0;
                                    0 0  (1-mu)/2];
                ju=false;

            elseif ques_type == 0   %the plane strain case
                D= (E/((1+mu)*(1-2*mu))) *[1-mu mu   0;
                                           mu   1-mu 0;
                                           0    0    (1-2*mu)/2];
                ju=false;
            else
                ques_type= input('The inputted Question Type is wrong.\nPlease input 1 (stands for the plane stress case) or 0 (stands for the plane strain case):');
            end
        end

        jun=false;

    elseif input_mode == 1
        %% Example 1 - from the page 109 of the textbook Variational principle and Finite Element Method 
        span=18;  height=3; % meters 
        t=0.4;  % meters 
        E= 2e6;  % kN/m^2
        mu=0.167;  % 1
        
        % input nodes data
        nodes_num= 7*13;  
        nodes= -1*ones(nodes_num, 2);
        X_new= transpose(0:0.75:9);
        Y_new= transpose(-1.5:0.5:1.5);
        for i=1:length(Y_new)    %length(Y)=7
            for j=1:length(X_new)  %length(X)=13
                nodes((j-1)*length(Y_new)+i,1)= X_new(j);
                nodes((j-1)*length(Y_new)+i,2)= Y_new(i);
            end
        end
        
        %input element details
        indices=1:83;         
        index= find(mod(indices,7)~=0);
        indices= indices(index);
        elements_num = 144;
        j=0;
        elements= zeros(elements_num,3);
        for i=indices
            j=j+1;
            elements(j,:)=[i i+7 i+8];
            j=j+1;
            elements(j,:)=[i i+8 i+1];
        end 

        q= 10; % kN/m^2
        P=zeros(nodes_num*2,1);
        P(2*88)=1.111;
        P(2*(1:7) - 1) = 1.111;
        P([2*7 2*91])=  -q*t*span/48; 
        P(2*(14:7:84))= -q*t*span*2/48;

        U=0.0011*ones(nodes_num*2,1);
        U(2*88)=0;
        U(2*(1:7) - 1)=0;
        
        ques_type = 1;
        D= (E/(1-mu*mu)) * [1 mu 0;
                           mu 1  0;
                            0 0  (1-mu)/2];
        
        jun=false;
        
    elseif input_mode == 2
        %% Example 2 - from page 223 of Reference
        E=2.1e8; mu= 0.3; t=0.025; 
        nodes_num= 4;
        nodes=[0 0;
             0.5 0;
             0.5 0.25;
             0   0.25];
        elements_num= 2;
        elements = [1 3 4;1 2 3];
        
        U=[0;0;0.0011;0.0011;0.0011;0.0011;0;0];
        P=[1.111;1.111;9.375;0;9.375;0;1.111;1.111];
        
        ques_type = 1;
        D= (E/(1-mu*mu)) * [1 mu 0;
                           mu 1  0;
                            0 0  (1-mu)/2];
        
        jun=false;
        
    else
        disp('Input error, only 0,1 or 2 is viable. Please input again.');
        input_mode= input('Please input 0 to input data by the keyboard,\nor input 1 to input data and show the results of Example 1 (from the page 109 of the textbook Variational principle and Finite Element Method) automatically,\nor input 2 to input data and show the results of Example 2 (from the page 223 of Reference) automatically: \n');
    end
end

    

        
        

%% ------ Step 2: Calculating the trianguler element area ------
n1=elements(1,1);  % node 1
n2=elements(1,2);
n3=elements(1,3);
Area=  CSTElementArea(nodes(n1,1),nodes(n1,2),nodes(n2,1),nodes(n2,2),nodes(n3,1),nodes(n3,2));

%% ------ Step 3: Calculating the element Stiffness Matrix ------       
k= zeros(6,6, elements_num);  % the element Stiffness Matrices
B= zeros(3,6, elements_num);  %  the Element Geometric Matrices (aka Element Strain Matrices) 
for i=1:elements_num
    n1=elements(i,1);  % node 1
    n2=elements(i,2);
    n3=elements(i,3);
    [k(:,:,i),B(:,:,i)]= CSTElementStiffness(D, t, nodes(n1,1),nodes(n1,2),nodes(n2,1),nodes(n2,2),nodes(n3,1),nodes(n3,2));
end 
 
%% ------ Step 4: Assembling element stiffness matrices ------
K= CSTElementStiffnessMatrixAssemble(k,nodes_num, elements);

%% ------ Step 5: Applying the boundary conditions and calculate the unknown displacements and reactions ------
% u= zeros(6,elements_num);

% Calculating unknown Nodal Displacements with Known nodal external forces
tag=find(P~=1.111);   %(p.s.'1.111' is the judgment criteria) find the indices of known external forces      
K_partition= K(tag,tag);  %partition the frame stiffness Matrix K
P_known=P(tag);   % find the known external forces
U_unknown=K_partition\P_known;  %the backslash operator ※\§ is used for Gaussian elimination
% U_unknown=inv(K_partition)*P_known
U(tag)=U_unknown;      % back to the U vector

% Calculate unknown reactions with Nodal Displacements 
tag2=find(P==1.111);  %find the indices of unknown external forces
K_partition2= K(tag2,:);  %
P_unknown= K_partition2 * U;  % Calculate unknown reactions
P(tag2)=P_unknown;   % back to the P vector 

%% ------ Step 6: Calculating the element stress vectors & Calculating the element principal stressess ------
sigma = zeros(3,elements_num);
principal_sigma = zeros(3,elements_num);
for i=1:elements_num
    n1=elements(i,1);  % node 1
    n2=elements(i,2);
    n3=elements(i,3);
    u=U([2*n1-1 2*n1 2*n2-1 2*n2 2*n3-1 2*n3]);   % u is the element displacement vector in size of 6x1 
    sigma(:,i)= CSTElementStresses(D,B(:,:,i),u);    % sigma(:,i) is the 3x1 stress vector for element i
    principal_sigma(:,i) = CSTElementPrincipalStresses(sigma(:,i)); % principal_sigma(:,i) is the 3x1 principal stress vector for element i
end

%% ------ Step 7: Plot the diagram of structure, displacements, and stresses ------
centroid_x= zeros(1,elements_num);
centroid_y= zeros(1,elements_num);

% figure(1);axis equal;set(figure(1),'Name','Displacement Diagram (magnification factor = 100)');

if input_mode == 1 
    m_factor= 100;          % magnification factor of displacement
elseif input_mode == 2
    m_factor = 10000;
else            % input_mode ==0 
    m_factor = input('\nInput the magnification factor of displacement which will be applied to plot the displacement diagram.\n p.s. The default value is 100, just input 100 if you do not know what to input: ');   
end

hold on; axis equal;
title('Displacement Diagram (magnification factor = '+string(m_factor )+')');
for c= 1:elements_num
    i=elements(c,1);    % node 1
    j=elements(c,2);    % node 2
    m=elements(c,3);    % node 3
    
    centroid_x(c)= mean(nodes([i j m],1));
    centroid_y(c)= mean(nodes([i j m],2));
    
    displacements_x= U(2*[i j m]-1);
    displacements_y= U(2*[i j m]);
         
    displacements  =  m_factor*[displacements_x displacements_y];  
    original_loc = nodes([i j m],:);   % original_loc == original location of element nodes
    CSTElementDSDiagram(original_loc ,displacements);
    
    if input_mode ==1
        % Draw a mirror diagram of the Y-axis
        mirror_loc= original_loc.*[-1 1];
        mirror_displacements = displacements.*[-1 1];
        CSTElementDSDiagram(mirror_loc ,mirror_displacements);
    end
     
end

% shading interp;
% colorbar;
toc

%% Get the result and Plot the Stress Field of Example 1 
if input_mode==1
    % Plot the Stress Field 
    SIGMAX_NODES = CSTElementStressField(centroid_x,centroid_y,sigma,nodes);

    % output the 考x of line x=0
    centroid_sigmax=sigma(1,12:-2:2);
    centroid_y_value=(1.5-0.5/3):-0.5:(-1.5+1/3);
    p=polyfit(centroid_y_value,centroid_sigmax,1);   % polyfit
    syms x
    result = p(1)*x+ p(2);
    sigma0=double(subs(result,x, 1.5:-0.5:-1.5));
    disp('The 考x at line x=0: ');
    disp(sigma0);
end


