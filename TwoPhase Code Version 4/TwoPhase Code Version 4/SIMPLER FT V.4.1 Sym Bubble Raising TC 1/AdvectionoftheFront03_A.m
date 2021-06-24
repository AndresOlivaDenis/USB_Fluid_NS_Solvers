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
    =AdvectionoftheFront03_A(NFront0,xFront0,yFront0,...
    X1,X2,Y1,Y2,Nx,Ny,u,v,u0,v0,xu,yu,xv,yv,dt,...
    xmapFront0,ymapFront0,FMx,FMy)

% Interpolation function -------------------------------------------------
FIx=@(r) 1-abs(r); % lineal Interpolation 
FIy=@(r) 1-abs(r); % lineal Interpolation
% ------------------------------------------------------------------------

% Initialize -------------------------------------------------------------
uFront=zeros(1,NFront0+2); vFront=zeros(1,NFront0+2);
umapFront=zeros(1,NFront0+2); vmapFront=zeros(1,NFront0+2);
xFront=zeros(1,NFront0+2); yFront=zeros(1,NFront0+2);
xmapFront=zeros(1,NFront0+2); ymapFront=zeros(1,NFront0+2);
uFront0=zeros(1,NFront0+2); vFront0=zeros(1,NFront0+2);
umapFront0=zeros(1,NFront0+2); vmapFront0=zeros(1,NFront0+2);
% ------------------------------------------------------------------------

% weigth=zeros(2,NFront0+2);
 
% ========================================================================
for il=2:NFront0+1
    %    Velocities of Front Points =====================================
    
    %  u Velocity -------------------------------------------------------
    ixF=floor(xmapFront0(il)*Nx+1); jyF=floor(ymapFront0(il)*Ny+1+0.5);
    dxu=xu(ixF+1,jyF)-xu(ixF,jyF); dyu=yu(ixF,jyF+1)-yu(ixF,jyF);
    
    rx1=(xFront0(il)-xu(ixF,jyF))/dxu; ry1=(yFront0(il)-yu(ixF,jyF))/dyu;
    rx2=(xu(ixF+1,jyF)-xFront0(il))/dxu; ry2=(yFront0(il)-yu(ixF+1,jyF))/dyu;
    rx3=(xFront0(il)-xu(ixF,jyF+1))/dxu; ry3=(yu(ixF,jyF+1)-yFront0(il))/dyu;
    rx4=(xu(ixF+1,jyF+1)-xFront0(il))/dxu ;ry4=(yu(ixF+1,jyF+1)-yFront0(il))/dyu;
    
    uFront(il)=u(ixF,jyF)*FIx(rx1)*FIy(ry1)+u(ixF+1,jyF)*FIx(rx2)*FIy(ry2)+...
        u(ixF,jyF+1)*FIx(rx3)*FIy(ry3)+u(ixF+1,jyF+1)*FIx(rx4)*FIy(ry4);
    umapFront(il)=uFront(il)/dxu*(1/Nx);
    
    uFront0(il)=u0(ixF,jyF)*FIx(rx1)*FIy(ry1)+u0(ixF+1,jyF)*FIx(rx2)*FIy(ry2)+...
        u0(ixF,jyF+1)*FIx(rx3)*FIy(ry3)+u0(ixF+1,jyF+1)*FIx(rx4)*FIy(ry4);
    umapFront0(il)=uFront0(il)/dxu*(1/Nx);
    % --------------------------------------------------------------------
    
    % v Velocity ---------------------------------------------------------
    ixF=floor(xmapFront0(il)*Nx+1+0.5); jyF=floor(ymapFront0(il)*Ny+1);
    dxv=xv(ixF+1,jyF)-xv(ixF,jyF); dyv=yv(ixF,jyF+1)-yv(ixF,jyF);
    
    rx1=(xFront0(il)-xv(ixF,jyF))/dxv; ry1=(yFront0(il)-yv(ixF,jyF))/dyv;
    rx2=(xv(ixF+1,jyF)-xFront0(il))/dxv; ry2=(yFront0(il)-yv(ixF+1,jyF))/dyv;
    rx3=(xFront0(il)-xv(ixF,jyF+1))/dxv; ry3=(yv(ixF,jyF+1)-yFront0(il))/dyv;
    rx4=(xv(ixF+1,jyF+1)-xFront0(il))/dxv; ry4=(yv(ixF+1,jyF+1)-yFront0(il))/dyv;
    
    vFront(il)=v(ixF,jyF)*FIx(rx1)*FIy(ry1)+v(ixF+1,jyF)*FIx(rx2)*FIy(ry2)+...
        v(ixF,jyF+1)*FIx(rx3)*FIy(ry3)+v(ixF+1,jyF+1)*FIx(rx4)*FIy(ry4);
    vmapFront(il)=vFront(il)/dyv*(1/Ny);
    
    vFront0(il)=v0(ixF,jyF)*FIx(rx1)*FIy(ry1)+v0(ixF+1,jyF)*FIx(rx2)*FIy(ry2)+...
        v0(ixF,jyF+1)*FIx(rx3)*FIy(ry3)+v0(ixF+1,jyF+1)*FIx(rx4)*FIy(ry4);
    vmapFront0(il)=vFront0(il)/dyv*(1/Ny);
    % --------------------------------------------------------------------
    
end
% WeigthSum=....
%     FIx(rx1)*FIy(ry1)+FIx(rx2)*FIy(ry2)+FIx(rx3)*FIy(ry3)+FIx(rx4)*FIy(ry4)

% Move of Front Points ===================================================
for il=2:NFront0+1
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

% ========================================================================
xFront(1)=-xFront(3);       xFront(NFront0+2)=-xFront(NFront0);
yFront(1)=yFront(3);        yFront(NFront0+2)=yFront(NFront0);

uFront(1)=-uFront(3);        uFront(NFront0+2)=-uFront(NFront0);
vFront(1)=vFront(3);        vFront(NFront0+2)=vFront(NFront0);

xmapFront(1)=-xmapFront(3);       xmapFront(NFront0+2)=-xmapFront(NFront0);
ymapFront(1)=ymapFront(3);        ymapFront(NFront0+2)=ymapFront(NFront0);

umapFront(1)=-umapFront(3);        umapFront(NFront0+2)=-umapFront(NFront0);
vmapFront(1)=vmapFront(3);        vmapFront(NFront0+2)=vmapFront(NFront0);
% ========================================================================

