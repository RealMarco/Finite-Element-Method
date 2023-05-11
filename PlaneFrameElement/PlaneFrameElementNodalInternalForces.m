% Date & Time: 2020/06/21 16:39
% Project: The Plane Frame 
% Step5: Calculate nodal internal forces of each beam
% Instruction£ºThis function calculates the nodal internal force vector of each PFE using the 
%              element stiffness matrixes in LOCAL coordinate system k_, coordinate Transformation Matrixes lamda,
%              the displacement vector U, and the beam vector. It returns the 6¡Ánumber_of_beams matrix, 
%              and f(:,i) stands for the nodal internal force vector of beam i. 
function f=PlaneFrameElementNodalInternalForces(k_ ,lamda,U,beam)
f=zeros(6,size(beam,1));
for i = 1:size(beam,1)
    node1 =beam(i,1);
    node2 =beam(i,2);
    u= U([3*node1-2 3*node1-1 3*node1 3*node2-2 3*node2-1 3*node2]);
    f(:,i)=k_(:,:,i) * lamda(:,:,i) * u;  % nodal internal force vector of beam i
end
