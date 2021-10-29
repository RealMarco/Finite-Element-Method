% Date & Time: 2020/07/07 08:01
% Project: The Eight Points Quadrilateral Element
% Step 6: Plot the diagram of structure, displacements, and stresses
% Instruction£ºPlot the displacement diagram of element and stress Heat Map
function  QuadraticQuadElementDSDiagram(original_loc, u,sigma_e)  % original_loc == original location of element nodes, u == element nodal displacements
new_loc = original_loc + u;  % new_loc == new location of element nodes
hold on 

fill(new_loc(:,1),new_loc(:,2),sigma_e,'FaceColor','interp');      % Plot stress heat map

for c=1:size(original_loc,1)-1
    % plot the line and point except the last ones 
    plot([original_loc(c,1) original_loc(c+1,1)],[original_loc(c,2) original_loc(c+1,2)],'k--','linewidth',0.1); 
    plot([new_loc(c,1) new_loc(c+1,1)],[new_loc(c,2) new_loc(c+1,2)],'r-');  
    plot(new_loc(c,1),new_loc(c,2),'ro','linewidth',0.1); 
end
% plot the last line and point
plot([original_loc(end,1) original_loc(1,1)],[original_loc(end,2) original_loc(1,2)],'k--','linewidth',0.1); 
plot([new_loc(end,1) new_loc(1,1)],[new_loc(end,2) new_loc(1,2)],'r-'); 
plot(new_loc(end,1),new_loc(end,2),'ro','linewidth',0.1); 
  
