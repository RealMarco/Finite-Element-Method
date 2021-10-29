% Date & Time: 2020/07/07 21:03
% Project: The Constant strain triangular Element
% Instruction£ºPlot the Stress Field
function SIGMAX_NODES=CSTElementStressField(centroid_x,centroid_y,sigma,nodes)
figure(2)
element_num= size(centroid_x,2);
X_ori = reshape(centroid_x(2:2:element_num),6,12);   % meshgrid manually
Y_ori = reshape(centroid_y(2:2:element_num),6,12);   % meshgrid manually
element_sigmax= reshape(sigma(1,2:2:element_num),6,12);
nodes_x=nodes(:,1);
nodes_y=nodes(:,2);
X_new=reshape(nodes_x,7,13);
Y_new=reshape(nodes_y,7,13);
SIGMAX_NODES= interp2(X_ori,Y_ori,element_sigmax,X_new,Y_new,'makima');   % 2D interpolation
surf(X_new,Y_new,SIGMAX_NODES); 
hold on
surf(-X_new,Y_new,SIGMAX_NODES); 
shading interp
xlabel('Length') ;
ylabel('Height') ;
zlabel('Stress ¦Òx');
title('Stress Field');

% X_ori,Y_ori,element_sigmax,X_new,Y_new
