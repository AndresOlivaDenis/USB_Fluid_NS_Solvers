% =======================================================================
% Plot Results
% =======================================================================

if (Plot1==1)
    nx=Nx; ny=Ny; h=(X2-X1)/Nx;
    uu(1:nx+1,1:ny+1)=0.5*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
    vv(1:nx+1,1:ny+1)=0.5*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
    xx(1:nx+1,1:ny+1)=0.5*(xu(1:nx+1,2:ny+2)+xu(1:nx+1,1:ny+1));
    yy(1:nx+1,1:ny+1)=0.5*(yv(2:nx+2,1:ny+1)+yv(1:nx+1,1:ny+1));
    w(1:nx+1,1:ny+1)=(u(1:nx+1,2:ny+2)-u(1:nx+1,1:ny+1)-...
        v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1))/(2*h);
    hold off;  subplot(4,4,[1 2 5 6 9 10 13 14]);% subplot(3,3,[1 2 4 5 7 8]);
    quiver(xi,yj,flipud(rot90(uu)),flipud(rot90(vv)),'r');
    axis equal; hold on;
    contour(flipud(rot90(xx)),flipud(rot90(yy)),flipud(rot90(w)),100);
    hold off;
    % ------------------------------ %
    ErrorNit(it)=Error; niteraNit(it)=nitera; subplot(4,4,[3 4]);
    plot(Nit(1:it),ErrorNit(1:it),'ok',Nit(1:it),ones(it,1)*Tol,'r',...
        'linewidth',2.5);
    xlabel('# Iteracion'); ylabel('Error'); subplot(4,4,[7 8]);
    plot(Nit(1:it),niteraNit(1:it),'ok',...
        Nit(1:it),ones(it,1)*nmax,'r','linewidth',2.5);
    xlabel('# Iteracion'); ylabel('#SubInt');
    pause(0.01);
    % ------------------------------ %
    %         ErrorNitPE(it)=ErrorPE; niteraNitPE(it)=niteraPE;
    %         ErrorNitPC(it)=ErrorPC; niteraNitPC(it)=niteraPC; subplot(4,4,[11 12]);
    %         plot(Nit(1:it),ErrorNitPE(1:it),'og',Nit(1:it),ErrorNitPC(1:it),'ob',...
    %             Nit(1:it),ones(it,1)*TolP,'r','linewidth',2.5);
    % ------------------------------ %
    xlabel('# Iteracion'); ylabel('Error'); subplot(4,4,[15 16]);
    plot(Nit(1:it),niteraNitPE(1:it),'og',Nit(1:it),niteraNitPC(1:it),'ob',...
        Nit(1:it),ones(it,1)*nmaxP,'r','linewidth',2.5);
    xlabel('# Iteracion'); ylabel('#SubInt');
    % ------------------------------ %
    %         uvels1(it)=uvels(2); vvels1(it)=vvels(2);
    %         subplot(4,4,[11 12]);
    %         plot(Nit(1:it),uvels1(1:it),'og',Nit(1:it),vvels1(1:it),'ob');
    %         pause(0.01);
    %
elseif (Plot1==2)
    nx=Nx; ny=Ny; h=(X2-X1)/Nx;
    uu(1:nx+1,1:ny+1)=0.5*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
    vv(1:nx+1,1:ny+1)=0.5*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
    xx(1:nx+1,1:ny+1)=0.5*(xu(1:nx+1,2:ny+2)+xu(1:nx+1,1:ny+1));
    yy(1:nx+1,1:ny+1)=0.5*(yv(2:nx+2,1:ny+1)+yv(1:nx+1,1:ny+1));
    w(1:nx+1,1:ny+1)=(u(1:nx+1,2:ny+2)-u(1:nx+1,1:ny+1)-...
        v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1))/(2*h);
    hold off;
    contour(xP,yP,ro); axis equal; hold on;
    quiver(xi,yj,flipud(rot90(uu)),flipud(rot90(vv)),'r');
    contour(flipud(rot90(xx)),flipud(rot90(yy)),flipud(rot90(w)),100);
    pause(0.01);
end