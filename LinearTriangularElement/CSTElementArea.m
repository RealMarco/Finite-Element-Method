% Date & Time: 2020/06/25 18:06
% Project: The Constant strain triangular Element
% Step 2: Calculating the trianguler element area
% Instruction£ºThis function returns the element area given the coordinates of the ?rst node(xi,yi),
%                                    the coordinates of the second node (xj, yj), and the coordinates of the 3rd node (xm, ym).

function Area=CSTElementArea(xi,yi,xj,yj,xm,ym)
Area=0.5*(xi*(yj-ym)+xj*(ym-yi)+xm*(yi-yj));
