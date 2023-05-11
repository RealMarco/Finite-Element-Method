% Date & Time: 2020/06/20 16:00
% Project: The Plane Frame 
% Step2 Calculate the element Stiffness Matrix
% Instruction£ºThis function returns the element length given the coordinates of the ?rst node(x1,y1) 
%              and the coordinates of the second node (x2, y2).
function L = BeamLength(x1,y1,x2,y2)
L = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));