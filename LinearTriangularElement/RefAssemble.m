% Date & Time: 2020/07/06 17:24
% Project: The Eight Points Quadrilateral Element
% Step 3: Assembling element stiffness matrices
% Instruction£ºThis function assembles the element stiffness matrices k into the global stiffness matrix K. 
%              It returns the 2n¡Á2n (n represents the number of nodes) global stiffness matrix K.

function K= QuadraticQuadAssemble(k,n, elements)
K=zeros(2*n,2*n);
k_expanded=zeros(2*n,2*n,size(elements,1));  %Add zeros to expand the order of the matrix from 16*16 to 2n*2n
for c =1:size(elements,1)       
    i=elements(c,1);
    j=elements(c,2);
    m=elements(c,3);
    p=elements(c,4);
    q=elements(c,5);
    r=elements(c,6);
    s=elements(c,7);
    t=elements(c,8);
    
    % assign the k_expanded(:,:,c)
    indices= [2*i-1 2*i, 2*j-1 2*j, 2*m-1 2*m, 2*p-1 2*p, 2*q-1 2*q, 2*r-1 2*r, 2*s-1 2*s, 2*t-1 2*t];
    k_expanded(indices, indices, c)= k(:,:,c);
    
    K=K+k_expanded(:,:,c); % assemble element stiffness matrix
end