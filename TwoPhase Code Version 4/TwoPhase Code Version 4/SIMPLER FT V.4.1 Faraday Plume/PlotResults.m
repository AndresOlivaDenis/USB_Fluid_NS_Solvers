% =======================================================================
% Plot of Results
% =======================================================================

if (Plot1==1) % =========================================================
    nx=Nx; ny=Ny; h=(X2-X1)/Nx;
    uu(1:nx+1,1:ny+1)=0.5*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
    vv(1:nx+1,1:ny+1)=0.5*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
    xx(1:nx+1,1:ny+1)=0.5*(xu(1:nx+1,2:ny+2)+xu(1:nx+1,1:ny+1));
    yy(1:nx+1,1:ny+1)=0.5*(yv(2:nx+2,1:ny+1)+yv(1:nx+1,1:ny+1));
    w(1:nx+1,1:ny+1)=(u(1:nx+1,2:ny+2)-u(1:nx+1,1:ny+1)-...
        v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1))/(2*h);
    hold off;  subplot(4,4,[1 2 5 6 9 10 13 14]);% subplot(3,3,[1 2 4 5 7 8]);
    quiver(xi,yj,flipud(rot90(uu)),flipud(rot90(vv)),'r');
    axis equal; hold on;  contour(xP,yP,ro);
    plot(xFront(1:NFront+1),yFront(1:NFront+1),'k','linewidth',2.5);pause(0.01)
    % contour(flipud(rot90(xx)),flipud(rot90(yy)),flipud(rot90(w)),100);
    hold off;
    % --------------------------------------------------------------------
    ErrorNit(it)=Error; niteraNit(it)=nitera; subplot(4,4,[3 4]);
    plot(Nit(1:it),ErrorNit(1:it),'ok',Nit(1:it),ones(it,1)*Tol,'--r',...
        'linewidth',2.5);
    xlabel('# Iteracion'); ylabel('Error'); subplot(4,4,[7 8]);
    plot(Nit(1:it),niteraNit(1:it),'ok',...
        Nit(1:it),ones(it,1)*nmax,'--r','linewidth',2.5);
    xlabel('# Iteracion'); ylabel('#SubInt');
    %         pause(0.01);
    % --------------------------------------------------------------------
    ErrorNitPE(it)=ErrorPE; niteraNitPE(it)=niteraPE;
    ErrorNitPC(it)=ErrorPC; niteraNitPC(it)=niteraPC; subplot(4,4,[11 12]);
    plot(Nit(1:it),ErrorNitPE(1:it),'og',Nit(1:it),ErrorNitPC(1:it),'ob',...
        Nit(1:it),ones(it,1)*TolPC,'--r',Nit(1:it),ones(it,1)*TolPE,'--r',...
        'linewidth',2.5);
    xlabel('# Iteracion'); ylabel('Error'); subplot(4,4,[15 16]);
    plot(Nit(1:it),niteraNitPE(1:it),'og',Nit(1:it),niteraNitPC(1:it),'ob',...
        Nit(1:it),ones(it,1)*nmaxP,'--r','linewidth',2.5);
    xlabel('# Iteracion'); ylabel('#SubInt');
    pause(0.005);
    % --------------------------------------------------------------------
elseif (Plot1==2) % ======================================================
    nx=Nx; ny=Ny; h=(X2-X1)/Nx;
    uu(1:nx+1,1:ny+1)=0.5*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
    vv(1:nx+1,1:ny+1)=0.5*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
    xx(1:nx+1,1:ny+1)=0.5*(xu(1:nx+1,2:ny+2)+xu(1:nx+1,1:ny+1));
    yy(1:nx+1,1:ny+1)=0.5*(yv(2:nx+2,1:ny+1)+yv(1:nx+1,1:ny+1));
    w(1:nx+1,1:ny+1)=(u(1:nx+1,2:ny+2)-u(1:nx+1,1:ny+1)-...
        v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1))/(2*h);
    hold off;
    quiver(xi,yj,flipud(rot90(uu)),flipud(rot90(vv)),'r');
    axis equal; hold on; % contour(xP,yP,ro,5);
    plot(xFront(1:NFront+1),yFront(1:NFront+1),'k','linewidth',1.5);pause(0.01)
    
    pause(0.005);
elseif (Plot1==3) % =====================================================
    nx=Nx; ny=Ny; h=(X2-X1)/Nx;
    uu(1:nx+1,1:ny+1)=0.5*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
    vv(1:nx+1,1:ny+1)=0.5*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
    xx(1:nx+1,1:ny+1)=0.5*(xu(1:nx+1,2:ny+2)+xu(1:nx+1,1:ny+1));
    yy(1:nx+1,1:ny+1)=0.5*(yv(2:nx+2,1:ny+1)+yv(1:nx+1,1:ny+1));
    w(1:nx+1,1:ny+1)=(u(1:nx+1,2:ny+2)-u(1:nx+1,1:ny+1)-...
        v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1))/(2*h);
    hold off;  subplot(3,3,[1 2 4 5 7 8]);% subplot(3,3,[1 2 4 5 7 8]);
    quiver(xi,yj,flipud(rot90(uu)),flipud(rot90(vv)),'r');
    axis equal; hold on;  contour(xP,yP,ro,5);
    plot(xFront(1:NFront+1),yFront(1:NFront+1),'k','linewidth',1.5);pause(0.01)
    % contour(flipud(rot90(xx)),flipud(rot90(yy)),flipud(rot90(w)),100);
    hold off;
    % --------------------------------------------------------------------
    subplot(3,3,3);
    plot(Nit(1:it)*dt,XcBubble(1:it,2),'k',Nit(1:it)*dt,XcBubble(1:it,2),'ok');
    xlabel('Tiempo [-]'); ylabel('yc [-]');
    subplot(3,3,6);
    plot(Nit(1:it)*dt,VelBubble(1:it,2),'b',Nit(1:it)*dt,VelBubble(1:it,2),'ob');
    xlabel('Tiempo [-]'); ylabel('Vely [-]');
    subplot(3,3,9);
    plot(Nit(1:it)*dt,CirBubble(1:it),'r',Nit(1:it)*dt,CirBubble(1:it),'or');
    xlabel('Tiempo [-]'); ylabel('Cir [-]');
    pause(0.005);
end
% ========================================================================