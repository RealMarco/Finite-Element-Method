
figure(1);axis equal;set(figure(1),'Name','Displacement Diagram (magnification factor = 100)');
figure(2);axis equal;set(figure(2),'Name','Stress Heat Map');
figure(1); hold on 
for c= 1:10
    TestPlot(c,c);
end 
plot(11,11,'bo');
figure(2);hold on
plot(12,12,'bo');
plot([0 8],[1 7],'k-');
