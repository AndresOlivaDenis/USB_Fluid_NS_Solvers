% ======================================================================= %
% Advection of The Front
% To interpolate Velocities
%     A -> Bilinear Interpolation (uses the nearest 4points)
%     B -> Peskin (1977) (uses the nearest 4points + 12)
%     C -> Peskin & McQueen(1989) (uses the nearest 4points + 12)
%     D -> Brackbill and Ruppel (1986) (uses the nearest 4points + 12)
% ======================================================================= %

function [xFront,yFront,uFront,vFront,...
    xmapFront,ymapFront,umapFront,vmapFront]...
    =AdvectionoftheFront03_B(NFront0,xFront0,yFront0,...
    X1,X2,Y1,Y2,Nx,Ny,u,v,u0,v0,xu,yu,xv,yv,dt,...
    xmapFront0,ymapFront0,FMx,FMy)

% Interpolation function -------------------------------------------------
FIx=@(r) 0.25*(1+cos(pi*r*0.5)); % Peskin (1977) Interpolation 
FIy=@(r) 0.25*(1+cos(pi*r*0.5)); % Peskin (1977) Interpolation
% ------------------------------------------------------------------------

% Initialize -------------------------------------------------------------
uFront=zeros(1,NFront0+2); vFront=zeros(1,NFront0+2);
umapFront=zeros(1,NFront0+2); vmapFront=zeros(1,NFront0+2);
xFront=zeros(1,NFront0+2); yFront=zeros(1,NFront0+2);
xmapFront=zeros(1,NFront0+2); ymapFront=zeros(1,NFront0+2);
uFront0=zeros(1,NFront0+2); vFront0=zeros(1,NFront0+2);
umapFront0=zeros(1,NFront0+2); vmapFront0=zeros(1,NFront0+2);
% weigth=zeros(2,NFront0+2);
% ------------------------------------------------------------------------

% ========================================================================
for il=1:NFront0+2
    %    Velocities of Front Points =====================================
    
    %  u Velocity -------------------------------------------------------
    ixF=floor(xmapFront0(il)*Nx+1); jyF=floor(ymapFront0(il)*Ny+1+0.5);
    dxu=xu(ixF+1,jyF)-xu(ixF,jyF); dyu=yu(ixF,jyF+1)-yu(ixF,jyF);
    
    rx1=(xFront0(il)-xu(ixF,jyF))/dxu; ry1=(yFront0(il)-yu(ixF,jyF))/dyu;
    rx2=(xu(ixF+1,jyF)-xFront0(il))/dxu; ry2=(yFront0(il)-yu(ixF+1,jyF))/dyu;
    rx3=(xFront0(il)-xu(ixF,jyF+1))/dxu; ry3=(yu(ixF,jyF+1)-yFront0(il))/dyu;
    rx4=(xu(ixF+1,jyF+1)-xFront0(il))/dxu ;ry4=(yu(ixF+1,jyF+1)-yFront0(il))/dyu;
    
    rAx1=rx1; rAy1=ry1+1; rAx2=rx2; rAy2=ry2+1;
    rAx3=rx1+1; rAy3=ry1; rAx4=rx2+1; rAy4=ry2;
    rAx5=rx3+1; rAy5=ry3; rAx6=rx4+1; rAy6=ry4;
    rAx7=rx3; rAy7=ry3+1; rAx8=rx4; rAy8=ry4+1;
    rAx9=rx1+1; rAy9=ry1+1; rAx10=rx2+1; rAy10=ry2+1;
    rAx11=rx3+1; rAy11=ry3+1; rAx12=rx4+1; rAy12=ry4+1;
    
    uFront(il)=u(ixF,jyF)*FIx(rx1)*FIy(ry1)+u(ixF+1,jyF)*FIx(rx2)*FIy(ry2)+...
        u(ixF,jyF+1)*FIx(rx3)*FIy(ry3)+u(ixF+1,jyF+1)*FIx(rx4)*FIy(ry4)+...
        u(ixF,jyF-1)*FIx(rAx1)*FIy(rAy1)+u(ixF+1,jyF-1)*FIx(rAx2)*FIy(rAy2)+...
        u(ixF-1,jyF)*FIx(rAx3)*FIy(rAy3)+u(ixF+2,jyF)*FIx(rAx4)*FIy(rAy4)+...
        u(ixF-1,jyF+1)*FIx(rAx5)*FIy(rAy5)+u(ixF+2,jyF+1)*FIx(rAx6)*FIy(rAy6)+...
        u(ixF,jyF+2)*FIx(rAx7)*FIy(rAy7)+u(ixF+1,jyF+2)*FIx(rAx8)*FIy(rAy8)+...
        u(ixF-1,jyF-1)*FIx(rAx9)*FIy(rAy9)+u(ixF+2,jyF-1)*FIx(rAx10)*FIy(rAy10)+...
        u(ixF-1,jyF+2)*FIx(rAx11)*FIy(rAy11)+u(ixF+2,jyF+2)*FIx(rAx12)*FIy(rAy12);
    umapFront(il)=uFront(il)/dxu*(1/Nx);
    
    uFront0(il)=u0(ixF,jyF)*FIx(rx1)*FIy(ry1)+u0(ixF+1,jyF)*FIx(rx2)*FIy(ry2)+...
        u0(ixF,jyF+1)*FIx(rx3)*FIy(ry3)+u0(ixF+1,jyF+1)*FIx(rx4)*FIy(ry4)+...
        u0(ixF,jyF-1)*FIx(rAx1)*FIy(rAy1)+u0(ixF+1,jyF-1)*FIx(rAx2)*FIy(rAy2)+...
        u0(ixF-1,jyF)*FIx(rAx3)*FIy(rAy3)+u0(ixF+2,jyF)*FIx(rAx4)*FIy(rAy4)+...
        u0(ixF-1,jyF+1)*FIx(rAx5)*FIy(rAy5)+u0(ixF+2,jyF+1)*FIx(rAx6)*FIy(rAy6)+...
        u0(ixF,jyF+2)*FIx(rAx7)*FIy(rAy7)+u0(ixF+1,jyF+2)*FIx(rAx8)*FIy(rAy8)+...
        u0(ixF-1,jyF-1)*FIx(rAx9)*FIy(rAy9)+u0(ixF+2,jyF-1)*FIx(rAx10)*FIy(rAy10)+...
        u0(ixF-1,jyF+2)*FIx(rAx11)*FIy(rAy11)+u0(ixF+2,jyF+2)*FIx(rAx12)*FIy(rAy12);
    umapFront0(il)=uFront0(il)/dxu*(1/Nx);
    % --------------------------------------------------------------------
    
    % v Velocity ---------------------------------------------------------
    ixF=floor(xmapFront0(il)*Nx+1+0.5); jyF=floor(ymapFront0(il)*Ny+1);
    dxv=xv(ixF+1,jyF)-xv(ixF,jyF); dyv=yv(ixF,jyF+1)-yv(ixF,jyF);
    
    rx1=(xFront0(il)-xv(ixF,jyF))/dxv; ry1=(yFront0(il)-yv(ixF,jyF))/dyv;
    rx2=(xv(ixF+1,jyF)-xFront0(il))/dxv; ry2=(yFront0(il)-yv(ixF+1,jyF))/dyv;
    rx3=(xFront0(il)-xv(ixF,jyF+1))/dxv; ry3=(yv(ixF,jyF+1)-yFront0(il))/dyv;
    rx4=(xv(ixF+1,jyF+1)-xFront0(il))/dxv; ry4=(yv(ixF+1,jyF+1)-yFront0(il))/dyv;
        
    rAx1=rx1; rAy1=ry1+1; rAx2=rx2; rAy2=ry2+1; 
    rAx3=rx1+1; rAy3=ry1; rAx4=rx2+1; rAy4=ry2;
    rAx5=rx3+1; rAy5=ry3; rAx6=rx4+1; rAy6=ry4;
    rAx7=rx3; rAy7=ry3+1; rAx8=rx4; rAy8=ry4+1;
    rAx9=rx1+1; rAy9=ry1+1; rAx10=rx2+1; rAy10=ry2+1;
    rAx11=rx3+1; rAy11=ry3+1; rAx12=rx4+1; rAy12=ry4+1;
    
    vFront(il)=v(ixF,jyF)*FIx(rx1)*FIy(ry1)+v(ixF+1,jyF)*FIx(rx2)*FIy(ry2)+...
        v(ixF,jyF+1)*FIx(rx3)*FIy(ry3)+v(ixF+1,jyF+1)*FIx(rx4)*FIy(ry4)+...
        v(ixF,jyF-1)*FIx(rAx1)*FIy(rAy1)+v(ixF+1,jyF-1)*FIx(rAx2)*FIy(rAy2)+...
        v(ixF-1,jyF)*FIx(rAx3)*FIy(rAy3)+v(ixF+2,jyF)*FIx(rAx4)*FIy(rAy4)+...
        v(ixF-1,jyF+1)*FIx(rAx5)*FIy(rAy5)+v(ixF+2,jyF+1)*FIx(rAx6)*FIy(rAy6)+...
        v(ixF,jyF+2)*FIx(rAx7)*FIy(rAy7)+v(ixF+1,jyF+2)*FIx(rAx8)*FIy(rAy8)+...
        v(ixF-1,jyF-1)*FIx(rAx9)*FIy(rAy9)+v(ixF+2,jyF-1)*FIx(rAx10)*FIy(rAy10)+...
        v(ixF-1,jyF+2)*FIx(rAx11)*FIy(rAy11)+v(ixF+2,jyF+2)*FIx(rAx12)*FIy(rAy12);
    vmapFront(il)=vFront(il)/dyv*(1/Ny);
    
    vFront0(il)=v0(ixF,jyF)*FIx(rx1)*FIy(ry1)+v0(ixF+1,jyF)*FIx(rx2)*FIy(ry2)+...
        v0(ixF,jyF+1)*FIx(rx3)*FIy(ry3)+v0(ixF+1,jyF+1)*FIx(rx4)*FIy(ry4)+...
        v0(ixF,jyF-1)*FIx(rAx1)*FIy(rAy1)+v0(ixF+1,jyF-1)*FIx(rAx2)*FIy(rAy2)+...
        v0(ixF-1,jyF)*FIx(rAx3)*FIy(rAy3)+v0(ixF+2,jyF)*FIx(rAx4)*FIy(rAy4)+...
        v0(ixF-1,jyF+1)*FIx(rAx5)*FIy(rAy5)+v0(ixF+2,jyF+1)*FIx(rAx6)*FIy(rAy6)+...
        v0(ixF,jyF+2)*FIx(rAx7)*FIy(rAy7)+v0(ixF+1,jyF+2)*FIx(rAx8)*FIy(rAy8)+...
        v0(ixF-1,jyF-1)*FIx(rAx9)*FIy(rAy9)+v0(ixF+2,jyF-1)*FIx(rAx10)*FIy(rAy10)+...
        v0(ixF-1,jyF+2)*FIx(rAx11)*FIy(rAy11)+v0(ixF+2,jyF+2)*FIx(rAx12)*FIy(rAy12);
    vmapFront0(il)=vFront0(il)/dyv*(1/Ny);
    % -------------------------------------------------------------------
    
%     rAx1=(xFront0(il)-xu(ixF,jyF-1))/dxu; rAy1=(yFront0(il)-yu(ixF,jyF-1))/dyu; 
%     rAx2=(xFront0(il)-xu(ixF+1,jyF-1))/dxu; rAy2=(yFront0(il)-yu(ixF+1,jyF-1))/dyu; 
%     rAx3=(xFront0(il)-xu(ixF-1,jyF))/dxu; rAy3=(yFront0(il)-yu(ixF-1,jyF))/dyu; 
%     rAx4=(xFront0(il)-xu(ixF+2,jyF))/dxu; rAy4=(yFront0(il)-yu(ixF+2,jyF))/dyu;
%     rAx5=(xFront0(il)-xu(ixF-1,jyF+1))/dxu; rAy5=(yFront0(il)-yu(ixF-1,jyF+1))/dyu; 
%     rAx6=(xFront0(il)-xu(ixF+2,jyF+1))/dxu; rAy6=(yFront0(il)-yu(ixF+2,jyF+1))/dyu;
%     rAx7=(xFront0(il)-xu(ixF,jyF+2))/dxu; rAy7=(yFront0(il)-yu(ixF,jyF+2))/dyu;
%     rAx8=(xFront0(il)-xu(ixF+1,jyF+2))/dxu; rAy8=(yFront0(il)-yu(ixF+1,jyF+2))/dyu;
%     
end

% --------------------------------------------------------------------
% WeigthSum=FIx(rx1)*FIy(ry1)+FIx(rx2)*FIy(ry2)+FIx(rx3)*FIy(ry3)+FIx(rx4)*FIy(ry4)+...
%     FIx(rAx1)*FIy(rAy1)+FIx(rAx2)*FIy(rAy2)+FIx(rAx3)*FIy(rAy3)+FIx(rAx4)*FIy(rAy4)+...
%     FIx(rAx5)*FIy(rAy5)+FIx(rAx6)*FIy(rAy6)+FIx(rAx7)*FIy(rAy7)+FIx(rAx8)*FIy(rAy8)+...
%     FIx(rAx9)*FIy(rAy9)+FIx(rAx10)*FIy(rAy10)+FIx(rAx11)*FIy(rAy11)+FIx(rAx12)*FIy(rAy12)

% --------------------------------------------------------------------

% Move of Front Points ===================================================
for il=1:NFront0+2
    xmapFront(il)=xmapFront0(il)+dt*(umapFront(il));
    ymapFront(il)=ymapFront0(il)+dt*(vmapFront(il));
    % Coordinates in Real Space
    xFront(il)=X1+(X2-X1)*FMx(xmapFront(il)); yFront(il)=Y1+(Y2-Y1)*FMy(ymapFront(il));
    
    
%     ixF=floor(xmapFront(il)*Nx+1); jyF=floor(ymapFront(il)*Ny+1);
%     xFront3(il)=(ixF+1-(xmapFront(il)*Nx+1))*xi(ixF)+...
%         ((xmapFront(il)*Nx+1)-ixF)*xi(ixF);
%     
%     xFront2(il)=xFront0(il)+dt*uFront(il); yFront2(il)=yFront0(il)+dt*vFront(il);
end
xFront(1)=xFront(NFront0+1); xFront(NFront0+2)=xFront(2);
yFront(1)=yFront(NFront0+1); yFront(NFront0+2)=yFront(2);

% ========================================================================

