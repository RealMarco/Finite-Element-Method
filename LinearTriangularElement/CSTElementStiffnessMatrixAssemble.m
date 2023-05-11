% Date & Time: 2020/06/26 09:52
% Project: The Constant strain triangular Element
% Step 4: Assembling element stiffness matrices
% Instruction£ºThis function assembles the element stiffness matrices k into the global stiffness matrix K. 
%              It returns the 2n¡Á2n (n represents the number of nodes) global stiffness matrix K.
function K=CSTElementStiffnessMatrixAssemble(k,n,elements)
K=zeros(2*n,2*n);
k_expanded=zeros(2*n,2*n,size(elements,1));  %Add zeros to expand the order of the matrix from 6*6 to 2n*2n
for c =1:size(elements,1)       
    i=elements(c,1);
    j=elements(c,2);
    m=elements(c,3);
    
    % assign the k_expanded(:,:,m)
    indices= [2*i-1 2*i, 2*j-1 2*j, 2*m-1 2*m];
    k_expanded(indices, indices, c)= k(:,:,c);
    
    K=K+k_expanded(:,:,c); % assemble element stiffness matrix
end
