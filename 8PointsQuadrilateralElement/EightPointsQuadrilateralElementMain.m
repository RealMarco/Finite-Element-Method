 % Date & Time: 2020/07/05 16:55
% Project: The Isoparametric Eight Points Quadrilateral Element
% Instruction: This is the MAIN fuction of the 8 Points Quadrilateral Element (aka The Quadratic Quadrilateral Element) project.

% Functions used in the 8 Nodes Quadrilateral Element (aka The Quadratic Quadrilateral Element) project:
% QuadraticQuadElementArea(x1, y1, x2, y2, x3, y3, x4, y4) 每 This function returns the element area given the coordinates of the ?rst node(x1,y1),the coordinates of the second node (x2, y2), the coordinates of the third node (x3, y3), and the coordinates of the fourth node (x4, y4).
% QuadraticQuadElementStiffness(D, t, x1, y1, x2, y2, x3, y3, x4, y4) 每 This function calculates the element stiffness matrix for each linear triangle with the Elasticity Matrix D, thickness t, 
%                                                                        and coordinates (x1, y1),(x2, y2),(x3, y3),(x4,y4) for the 1st, 2nd, 3rd, 4th element corner node respectively.  
%                                                                        It returns the 16℅16 element stiffness matrix k and the 3x16 Geometric Matrix B.
% QuadraticQuadAssemble(k,nodes_num, elements) 每This function assembles the element stiffness matrices k into the global stiffness matrix K. 
%                                                It returns the 2n℅2n (n represents the number of nodes) global stiffness matrix K.
% QuadraticQuadElementStresses(D,B,u) 每  This function calculates the element stresses using the Elasticity Matrix D, the Geometric Matrix B(:,:,i) 
%                                         and the element displacement vector u. It returns two 3*1 results (in the form of [sigma_x sigma_y tao_xy]T) 每 the general quadratic stress functions in 缶 and 灰, 
%                                         and the numerical values of the stresses at the centroid of the element.
% QuadraticQuadElementPStresses(sigma)每 This function calculates the element principal stresses using the element stress vector sigma. 
%                                        It returns a 3℅1 vector in the form [sigma1 sigma2 theta]T where sigma1 and sigma2 are the principal stresses for the element 
%                                        and theta is the principal angle.

clear;
tic
%% ------ Step 1: Discretizing the Structure (Input the details of each element) ------
disp('Hello, it is The Isoparametric Eight Points Quadrilateral Element (aka The Isoparametric Quadratic Quadrilateral Element) project. Units: m, kN');
disp('------Step1: Discretize the Frame manually & Input the details of each element-----');
input_mode= input('Please input 0 to input data by the keyboard,\nor input 1 to input data and show the results of Example 1 (from the page 109 of the textbook Variational Principle and Finite Element Method) automatically,\nor input 2 to input data and show the results of Example 2 (from the page 266 of Reference) automatically: \n');

jun=true;
while jun
    if input_mode == 0
        %% Input data by the keyboard.
        fprintf('\nInput data by the keyboard.\nPlease note that You should discretize the structure, number the elements and number the nodes manually.\n');
 
        % Input the coordinates of the nodes
        nodes= input('\nInput the Nodal Coordinate Matrix which is in the size of n℅2, while n is the total nodal number of the structure. \nEach row vector is in the form [x coordinate, y coordinate] :\n ');
        % Input the connectivity list of elements  
        elements= input('\nInput the b℅8 Connectivity List of the Elements, while b is the total number of elements. \nEach row vector is in the form serial number of [node1, node2, node3, node4, node5, node6, node7, node8] :\n');
        
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
        jun=true;
        while jun
            if ques_type == 1   % the plane stress case
                D= (E/(1-mu*mu)) * [1 mu 0;
                                   mu 1  0;
                                    0 0  (1-mu)/2];
                jun=false;

            elseif ques_type == 0   %the plane strain case
                D= (E/((1+mu)*(1-2*mu))) *[1-mu mu   0;
                                           mu   1-mu 0;
                                           0    0    (1-2*mu)/2];
                jun=false;
            else
                ques_type= input('The inputted Question Type is wrong.\nPlease input 1 (stands for the plane stress case) or 0 (stands for the plane strain case):');
            end
        end

        jun=false;

        
    elseif input_mode == 1
        %% Example 1 - from the page 109 of the textbook Variational principle and Finite Element Method 
        span=18;  height=3; 
        t=0.4;  % meters 
        E= 2e6;  % kN/m^2
        mu=0.167;  % 1

        nodes_num= 73;  % input nodes data
        nodes= -1*ones(nodes_num, 2);
        X= transpose(0:1.5:9);
        Y= transpose(-1.5:0.5:1.5);
        for i=1:length(Y)    %length(Y)=7
            for j=1:length(X)  %length(X)=7
                nodes((j-1)*11+i,1)= X(j);
                nodes((j-1)*11+i,2)= Y(i);
            end
        end

        X2= transpose(0.75:1.5:8.25);
        Y2= transpose(-1.5:1:1.5);
        for i=1:length(Y2)
            for j= 1:length(X2)
                nodes((j-1)*11+7+i,1)=X2(j);
                nodes((j-1)*11+7+i,2)=Y2(i);
            end
        end

        elements_num= 18;
        elements= zeros(elements_num,8);
        j=0;
        for i= 1:11:56
            j=j+1;
            elements(j,:)= [i   i+11 i+13 i+2, i+7 i+12 i+8  i+1];
            j=j+1;
            elements(j,:)= [i+2 i+13 i+15 i+4, i+8 i+14 i+9  i+3];
            j=j+1;
            elements(j,:)= [i+4 i+15 i+17 i+6, i+9 i+16 i+10 i+5];
        end

        % boundary conditions and external forces
        q= 10; % kN/m^2
        P=zeros(nodes_num*2,1);
        P(2*(1:7)-1)=1.111;
        P(2*70) = 1.111;
        P([2*7 2*73])         =   -q*t*span/72; 
        P(2*(11:11:66))= -q*t*span*4/72;
        P(2*(18:11:62))=    -q*t*span*2/72;

        U=0.0011*ones(nodes_num*2,1);
        U(2*(1:7)-1)=0;
        U(2*70)=0;
        
        ques_type = 1;
        D= (E/(1-mu*mu)) * [1 mu 0;
                           mu 1  0;
                            0 0  (1-mu)/2];
        jun=false;
    elseif input_mode == 2
        %% Example 2 - from page 266 of Reference
        E=2.1e8; mu= 0.3; t=0.025; 
        nodes_num= 8;
        nodes=[0 0;
             0.5 0;
             0.5 0.25;
             0   0.25;
            0.25 0;
            0.5  0.125;
            0.25 0.25;
            0    0.125];
        elements_num= 1;
        elements = [1 2 3 4 5 6 7 8];
        
        U=[0;0; 0.0011;0.0011 ;0.0011;0.0011 ;0;0 ;0.0011;0.0011; 0.0011;0.0011; 0.0011;0.0011; 0; 0];
        P=[1.111;1.111;3.125;0;3.125;0;1.111;1.111; 0;0; 12.5;0; 0;0; 1.111;1.111];
        
        ques_type = 1;
        D= (E/(1-mu*mu)) * [1 mu 0;
                           mu 1  0;
                            0 0  (1-mu)/2];
        jun=false;
    else
        disp('Input error, only 0,1 or 2 is viable. Please input again.');
        input_mode= input('Please input 0 to input data by the keyboard,\nor input 1 to input data and show the results of Example 1 (from the page 109 of the textbook Variational principle and Finite Element Method) automatically,\nor input 2 to input data and show the results of Example 2 (from the page 266 of Reference) automatically: \n');
    end
    
end


%% ------ Step 2: Calculating the element Stiffness Matrix ------
k= zeros(16,16, elements_num);  % the element Stiffness Matrices
B= sym(zeros(3,16, elements_num));  %  the Element Geometric Matrices (aka Element Strain Matrices) 
syms ksi eta

for i=1:elements_num
    n1=elements(i,1);  % node 1
    n2=elements(i,2);
    n3=elements(i,3);
    n4=elements(i,4);
    [k(:,:,i),B(:,:,i)]= QuadraticQuadElementStiffness(D, t, nodes(n1,1),nodes(n1,2),nodes(n2,1),nodes(n2,2),nodes(n3,1),nodes(n3,2),nodes(n4,1),nodes(n4,2));
end 
 
%% ------ Step 3: Assembling element stiffness matrices ------
K= QuadraticQuadAssemble(k,nodes_num, elements);

%% ------ Step 4: Applying the boundary conditions and calculate the unknown displacements and reactions ------
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


%% ------ Step 5: Calculating the element stress vectors & Calculating the element principal stresses------
sigma = sym(zeros(3,elements_num));   
centroid_sigma= zeros(3,elements_num);
principal_sigma = zeros(3,elements_num);
for c=1:elements_num
    i=elements(c,1);    % node 1
    j=elements(c,2);    % node 2
    m=elements(c,3);    % node 3
    p=elements(c,4);    % node 4
    q=elements(c,5);    % node 5
    r=elements(c,6);    % node 6
    s=elements(c,7);    % node 7
    t=elements(c,8);    % node 8
    
    indices= [2*i-1 2*i, 2*j-1 2*j, 2*m-1 2*m, 2*p-1 2*p, 2*q-1 2*q, 2*r-1 2*r, 2*s-1 2*s, 2*t-1 2*t];
    u=U(indices);   % u is the element displacement vector in size of 16x1 
%     sigma_sym(:,c)=D*B(:,:,c)*u;
%     centroid_sigma(:,c) = subs(sigma_sym(:,c), {ksi, eta},{0,0});
%     centroid_sigma(:,c)= double(centroid_sigma );
    
    [sigma(:,c), centroid_sigma(:,c)]= QuadraticQuadElementStresses(D,B(:,:,c),u);    % 
    principal_sigma(:,c) = QuadraticQuadElementPStresses(centroid_sigma(:,c)); % principal_sigma(:,i) is the 3x1 principal stress vector for element i
end

%% ------ Step 6: Plot the diagram of structure, displacements, and stresses ------
if input_mode == 1 
    m_factor= 100;          % magnification factor of displacement
elseif input_mode == 2
    m_factor = 10000;
else            % input_mode ==0 
    m_factor = input('\nInput the magnification factor of displacement which will be applied to plot the displacement diagram.\n p.s. The default value is 100, just input 100 if you do not know what to input: ');   
end

hold on; axis equal;
title('Displacement Diagram (magnification factor = '+string(m_factor )+')& Stress 考x Heat Map');

for c= 1:elements_num
    i=elements(c,1);    % node 1
    j=elements(c,2);    % node 2
    m=elements(c,3);    % node 3
    p=elements(c,4);    % node 4
    q=elements(c,5);    % node 5
    r=elements(c,6);    % node 6
    s=elements(c,7);    % node 7
    t=elements(c,8);    % node 8
    
    displacements_x= U(2*[i q j r m s p t]-1);
    displacements_y= U(2*[i q j r m s p t]);
    displacements  =  m_factor*[displacements_x displacements_y];  
    original_loc = nodes([i q j r m s p t],:);   % original_loc == original location of element nodes
    sigma_e=zeros(8,1);         
    sigma_e(1)= double(subs(sigma(1,c), {ksi, eta},{-1,-1}));
    sigma_e(2)= double(subs(sigma(1,c), {ksi, eta},{0,-1}));
    sigma_e(3)= double(subs(sigma(1,c), {ksi, eta},{1,-1}));
    sigma_e(4)= double(subs(sigma(1,c), {ksi, eta},{1,0}));
    sigma_e(5)= double(subs(sigma(1,c), {ksi, eta},{1,1}));
    sigma_e(6)= double(subs(sigma(1,c), {ksi, eta},{0,1}));
    sigma_e(7)= double(subs(sigma(1,c), {ksi, eta},{-1,1}));
    sigma_e(8)= double(subs(sigma(1,c), {ksi, eta},{-1,0}));
    
    QuadraticQuadElementDSDiagram(original_loc ,displacements,sigma_e);
    
    if input_mode == 1
        % Draw a mirror diagram of the Y-axis
        mirror_loc= original_loc.*[-1 1];
        mirror_displacements = displacements.*[-1 1];
        QuadraticQuadElementDSDiagram(mirror_loc ,mirror_displacements,sigma_e);
    end
    
%     QuadraticQuadElementStressHeatmap();
end

shading interp;
colorbar;

toc

%% Get the the 考x of line x=0 in Example 1
if input_mode == 1
    sigma0=zeros(3,7);
    sigma0(:,1)= double(subs(sigma(:,3), {ksi, eta},{-1,1}));
    sigma0(:,2)= double(subs(sigma(:,3), {ksi, eta},{-1,0}));
    % sigma0(:,3)=  (double(subs(sigma_sym(:,3), {ksi, eta},{-1,-1}))+ double(subs(sigma_sym(:,2), {ksi, eta},{-1,1})))/2;
    sigma0(:,3)= double(subs(sigma(:,2), {ksi, eta},{-1,1}));
    sigma0(:,4)= double(subs(sigma(:,2), {ksi, eta},{-1,0}));
    sigma0(:,5)= double(subs(sigma(:,2), {ksi, eta},{-1,-1}));
    % sigma0(:,5)= (double(subs(sigma_sym(:,1), {ksi, eta},{-1,1}))+ double(subs(sigma_sym(:,2), {ksi, eta},{-1,-1})))/2;
    sigma0(:,6)= double(subs(sigma(:,1), {ksi, eta},{-1,0}));
    sigma0(:,7)= double(subs(sigma(:,1), {ksi, eta},{-1,-1}));
    disp('The stresses at line x=0 in the form of transpose([考x, 考y, 而xy]): /n');
    disp(sigma0);
end

