% Date & Time: 2020/06/21 19:38
% Project: The Plane Frame 
% Step6: Plot the internal force diagrams
% Instruction£ºThis function plots the axial force diagram for the element with
%              object node, beam, coordinate Transformation Matrixes lamda,and internal force vectors F .
function y=PlaneFrameElementInternalForceDiagram(beam,node,lamda,F)
if size(beam,1)>0
    y=1;
else
    y=0;
end

figure
tab1 = uitab('Title','Axial Force Diagram (N)');
ax1 = axes(tab1);
axis off;
tab2 = uitab('Title','Shear Force Diagram (N)');
ax2 = axes(tab2);
axis off;
tab3 = uitab('Title','Bending moment Diagram (N¡¤cm)');
ax3 = axes(tab3);
axis off;

for i=1:size(beam,1)
    
    x1= node(beam(i,1),1);
    y1= node(beam(i,1),2);
    x2= node(beam(i,2),1);
    y2= node(beam(i,2),2);
    
    % control the order of magnitude in order to get a better diagram 
    f=F(:,i);
    % f=[-4288 5008 286328 4288 -5008 214428];
    t=log10(abs(f));
    scale_f=floor(t);
    f=f ./ (2*10.^(scale_f-1));
    
    % Plot the axial force diagram
    ua_=[0; -f(1); 0; 0 ; f(4);0];
    T=lamda(:,:,i);
    ua=T\ua_;
    
    Px1a=x1+ua(1);
    Py1a=y1+ua(2);
    Px2a=x2+ua(4);
    Py2a=y2+ua(5);
    
    hold on;
    plot(ax1,[x1 x2],[y1 y2],'b','linewidth',3);
   
    plot(ax1,[Px1a Px2a],[Py1a Py2a],'k','linewidth',1);
    plot(ax1,[x1 Px1a],[y1 Py1a],'k','linewidth',0.5);
    plot(ax1,[x2 Px2a],[y2 Py2a],'k','linewidth',0.5);
    text(ax1,(Px1a+Px2a)/2,(Py1a+Py2a)/2, num2str(F(4,i)) );
    
    %Plot the shear force diagram
    us_=[0; f(2); 0; 0 ; -f(5);0];
    T=lamda(:,:,i);
    us=T\us_;
    
    Px1s=x1+us(1);
    Py1s=y1+us(2);
    Px2s=x2+us(4);
    Py2s=y2+us(5);
    
    hold on;
    plot(ax2,[x1 x2],[y1 y2],'b','linewidth',3);
    plot(ax2,[Px1s Px2s],[Py1s Py2s],'k','linewidth',1);
    plot(ax2,[x1 Px1s],[y1 Py1s],'k','linewidth',0.5);
    plot(ax2,[x2 Px2s],[y2 Py2s],'k','linewidth',0.5);
    text(ax2,(Px1s+Px2s)/2,(Py1s+Py2s)/2, num2str(F(2,i)) );
    
    %Plot the bending moment diagram
    ub_=[0; -f(3); 0; 0 ; f(6);0];
    T=lamda(:,:,i);
    ub=T\ub_;

    Px1b=x1+ub(1);
    Py1b=y1+ub(2);
    Px2b=x2+ub(4);
    Py2b=y2+ub(5);

    hold on;
    plot(ax3,[x1 x2],[y1 y2],'b','linewidth',3);
    plot(ax3,[Px1b Px2b],[Py1b Py2b],'k','linewidth',1);
    plot(ax3,[x1 Px1b],[y1 Py1b],'k','linewidth',0.5);
    plot(ax3,[x2 Px2b],[y2 Py2b],'k','linewidth',0.5);
    text(ax3,Px1b,Py1b, num2str(-F(3,i)) );
    text(ax3,Px2b,Py2b, num2str(F(6,i)) );
    
end