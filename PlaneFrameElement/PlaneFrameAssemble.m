% Date & Time: 2020/06/20 21:11
% Project: The Plane Frame 
% Step3: Assemble element stiffness matrixes and get the Frame Stiffness Matrix
% Instruction£ºThis function assembles every element stiffness matrix k(:,:,m) into the global stiffness matrix K. 
%              It returns the 3n¡Á3n (n=number of nodes) global stiffness matrix K

function K=PlaneFrameAssemble(k,n,beam)
K=zeros(3*n,3*n);
k_expanded=zeros(3*n,3*n,size(beam,1));  %Add zeros to expand the order of the matrix from 6*6 to 3n*3n
for m =1:size(beam,1)       
    i=beam(m,1);
    j=beam(m,2);
    
    % assign the k_expanded(:,:,m)
    indices= [3*i-2 3*i-1 3*i 3*j-2 3*j-1 3*j];
    k_expanded(indices, indices, m)= k(:,:,m);

    K=K+k_expanded(:,:,m); % assemble element stiffness matrix
end




