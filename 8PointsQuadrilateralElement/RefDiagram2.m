% Date & Time: 2020/06/21 22:44
% Project: The Plane Frame 
% Step6: Plot the internal force diagrams
% Instruction£ºThis function plots the axial force diagram for the element with
%              object node, beam, coordinate Transformation Matrixes lamda,and internal force vectors F.
%Improvement:  Using different figures respectively instead of tabs in a figure.
function y=PlaneFrameElementInternalForceDiagram2(beam,node,lamda,F)
if size(beam,1)>0
    y=1;
else
    y=0;
end

figure(1);axis equal;axis off;  set(figure(1),'Name','Axial Force Diagram (N)');
figure(2);axis equal;axis off;  set(figure(2),'Name','Shear Force Diagram (N)');
figure(3);axis equal;axis off;  set(figure(3),'Name','Bending moment Diagram (N¡¤cm)');

% calculate the minimum order of magnitude of the 3 types of internal forces
t=log10(abs(F));
scale_f=floor(t); % order of magnitude
c1= min([scale_f(1,:),scale_f(4,:)]); % minimum order of magnitude of the axial forces
c2= min([scale_f(2,:),scale_f(5,:)]); % minimum order of magnitude of the shear forces
c3= min([scale_f(3,:),scale_f(6,:)]); % minimum order of magnitude of the bending moments
scale_min = [c1;c2;c3;c1;c2;c3];

for i=1:size(beam,1)
    
    x1= node(beam(i,1),1);
    y1= node(beam(i,1),2);
    x2= node(beam(i,2),1);
    y2= node(beam(i,2),2);
    
    % control the order of magnitude in order to get better diagrams 
    f=F(:,i);
    % F(:,1)=[-4288 ;5008 ;286328 ;4288; -5008 ;214428] in example 1
    f=f ./ (2*10.^(scale_min-1));
    
    % Plot the axial force diagram
    figure(1);
    ua_=[0; -f(1); 0; 0 ; f(4);0];
    T=lamda(:,:,i);
    ua=T\ua_;
    
    Px1a=x1+ua(1);
    Py1a=y1+ua(2);
    Px2a=x2+ua(4);
    Py2a=y2+ua(5);
    
    hold on;
    plot([x1 x2],[y1 y2],'b','linewidth',3);
    plot([Px1a Px2a],[Py1a Py2a],'k','linewidth',1);
    plot([x1 Px1a],[y1 Py1a],'k','linewidth',0.5);
    plot([x2 Px2a],[y2 Py2a],'k','linewidth',0.5);
    text((Px1a+Px2a)/2,(Py1a+Py2a)/2, num2str(F(4,i)),'VerticalAlignment','top','HorizontalAlignment','right' );
    
    %Plot the shear force diagram
    figure(2);
    us_=[0; f(2); 0; 0 ; -f(5);0];
    T=lamda(:,:,i);
    us=T\us_;
    
    Px1s=x1+us(1);
    Py1s=y1+us(2);
    Px2s=x2+us(4);
    Py2s=y2+us(5);
    
    hold on;
    plot([x1 x2],[y1 y2],'b','linewidth',3);
    plot([Px1s Px2s],[Py1s Py2s],'k','linewidth',1);
    plot([x1 Px1s],[y1 Py1s],'k','linewidth',0.5);
    plot([x2 Px2s],[y2 Py2s],'k','linewidth',0.5);
    text((Px1s+Px2s)/2,(Py1s+Py2s)/2, num2str(F(2,i)),'VerticalAlignment','top','HorizontalAlignment','right' );
    
    %Plot the bending moment diagram
    figure(3);
    ub_=[0; -f(3); 0; 0 ; f(6);0];
    T=lamda(:,:,i);
    ub=T\ub_;

    Px1b=x1+ub(1);
    Py1b=y1+ub(2);
    Px2b=x2+ub(4);
    Py2b=y2+ub(5);

    hold on;
    plot([x1 x2],[y1 y2],'b','linewidth',3);
    plot([Px1b Px2b],[Py1b Py2b],'k','linewidth',1);
    plot([x1 Px1b],[y1 Py1b],'k','linewidth',0.5);
    plot([x2 Px2b],[y2 Py2b],'k','linewidth',0.5);
    text(Px1b,Py1b, num2str(-F(3,i)),'VerticalAlignment','top','HorizontalAlignment','right');
    text(Px2b,Py2b, num2str(F(6,i)) ,'VerticalAlignment','top','HorizontalAlignment','right' );
    
end