% =======================================================================
% Post Plots results
% =======================================================================

nx=Nx; ny=Ny; h=(X2-X1)/Nx;
uu(1:nx+1,1:ny+1)=0.5*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
vv(1:nx+1,1:ny+1)=0.5*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));

% xx(1:nx+1,1:ny+1)=0.5*(xu(1:nx+1,2:ny+2)+xu(1:nx+1,1:ny+1));
% yy(1:nx+1,1:ny+1)=0.5*(yv(2:nx+2,1:ny+1)+yv(1:nx+1,1:ny+1));

w(1:nx+1,1:ny+1)=(u(1:nx+1,2:ny+2)-u(1:nx+1,1:ny+1)-...
    v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1))/(2*h);
figure; hold off,quiver(xi,yj,flipud(rot90(uu)),flipud(rot90(vv)),'r');hold on;
hold on; axis equal; contour(xP,yP,ro); grid on;
% contour(flipud(rot90(xx)),flipud(rot90(yy)),flipud(rot90(w)),100),axis equal
figure; 
plot(xFront(1:NFront+1),yFront(1:NFront+1),'k','linewidth',1);
axis([0.1,0.9,0.8,1.4]); axis equal; grid on
figure; plot(Nit,VelBubble(:,2)','k','linewidth',1);
figure; plot(Nit,XcBubble(:,2)','k','linewidth',1);
figure; plot(Nit,CirBubble','k','linewidth',1);

% =======================================================================

