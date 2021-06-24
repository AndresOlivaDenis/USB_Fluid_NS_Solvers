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
    xmapFront0,ymapFront0,FMx,FMy,Lx,Ly)

% Interpolation function -------------------------------------------------
FIx=@(r) 1-abs(r); % lineal Interpolation 
FIy=@(r) 1-abs(r); % lineal Interpolation
% ------------------------------------------------------------------------

% Initialize -------------------------------------------------------------
uFront{1}=zeros(1,NFront0{1}+2);    vFront{1}=zeros(1,NFront0{1}+2);
umapFront{1}=zeros(1,NFront0{1}+2); vmapFront{1}=zeros(1,NFront0{1}+2);
xFront{1}=zeros(1,NFront0{1}+2);    yFront{1}=zeros(1,NFront0{1}+2);
xmapFront{1}=zeros(1,NFront0{1}+2); ymapFront{1}=zeros(1,NFront0{1}+2);
uFront0{1}=zeros(1,NFront0{1}+2);   vFront0{1}=zeros(1,NFront0{1}+2);
umapFront0{1}=zeros(1,NFront0{1}+2); vmapFront0{1}=zeros(1,NFront0{1}+2);

uFront{2}=zeros(1,NFront0{2}+2);    vFront{2}=zeros(1,NFront0{2}+2);
umapFront{2}=zeros(1,NFront0{2}+2); vmapFront{2}=zeros(1,NFront0{2}+2);
xFront{2}=zeros(1,NFront0{2}+2);    yFront{2}=zeros(1,NFront0{2}+2);
xmapFront{2}=zeros(1,NFront0{2}+2); ymapFront{2}=zeros(1,NFront0{2}+2);
uFront0{2}=zeros(1,NFront0{2}+2);   vFront0{2}=zeros(1,NFront0{2}+2);
umapFront0{2}=zeros(1,NFront0{2}+2); vmapFront0{2}=zeros(1,NFront0{2}+2);

% weigth=zeros(2,NFront0+2);
 
% ------------------------------------------------------------------------

% ========================================================================
for i =1 : 2
    for il=2:NFront0{i}+1
        %    Velocities of Front Points =====================================
        
        %  u Velocity -------------------------------------------------------
        ixF=floor(xmapFront0{i}(il)*Nx+1);   ixF = min(ixF,Nx); 
        jyF=floor(ymapFront0{i}(il)*Ny+1+0.5);
        dxu=xu(ixF+1,jyF)-xu(ixF,jyF); dyu=yu(ixF,jyF+1)-yu(ixF,jyF);
        
        rx1=(xFront0{i}(il)-xu(ixF,jyF))/dxu; ry1=(yFront0{i}(il)-yu(ixF,jyF))/dyu;
        rx2=(xu(ixF+1,jyF)-xFront0{i}(il))/dxu; ry2=(yFront0{i}(il)-yu(ixF+1,jyF))/dyu;
        rx3=(xFront0{i}(il)-xu(ixF,jyF+1))/dxu; ry3=(yu(ixF,jyF+1)-yFront0{i}(il))/dyu;
        rx4=(xu(ixF+1,jyF+1)-xFront0{i}(il))/dxu ;ry4=(yu(ixF+1,jyF+1)-yFront0{i}(il))/dyu;
        
        uFront{i}(il)=u(ixF,jyF)*FIx(rx1)*FIy(ry1)+u(ixF+1,jyF)*FIx(rx2)*FIy(ry2)+...
            u(ixF,jyF+1)*FIx(rx3)*FIy(ry3)+u(ixF+1,jyF+1)*FIx(rx4)*FIy(ry4);
        umapFront{i}(il)=uFront{i}(il)/dxu*(1/Nx);
        
        uFront0{i}(il)=u0(ixF,jyF)*FIx(rx1)*FIy(ry1)+u0(ixF+1,jyF)*FIx(rx2)*FIy(ry2)+...
            u0(ixF,jyF+1)*FIx(rx3)*FIy(ry3)+u0(ixF+1,jyF+1)*FIx(rx4)*FIy(ry4);
        umapFront0{i}(il)=uFront0{i}(il)/dxu*(1/Nx);
        % --------------------------------------------------------------------
        
        % v Velocity ---------------------------------------------------------
        ixF=floor(xmapFront0{i}(il)*Nx+1+0.5);  jyF=floor(ymapFront0{i}(il)*Ny+1);
        dxv=xv(ixF+1,jyF)-xv(ixF,jyF); dyv=yv(ixF,jyF+1)-yv(ixF,jyF);
        
        rx1=(xFront0{i}(il)-xv(ixF,jyF))/dxv; ry1=(yFront0{i}(il)-yv(ixF,jyF))/dyv;
        rx2=(xv(ixF+1,jyF)-xFront0{i}(il))/dxv; ry2=(yFront0{i}(il)-yv(ixF+1,jyF))/dyv;
        rx3=(xFront0{i}(il)-xv(ixF,jyF+1))/dxv; ry3=(yv(ixF,jyF+1)-yFront0{i}(il))/dyv;
        rx4=(xv(ixF+1,jyF+1)-xFront0{i}(il))/dxv; ry4=(yv(ixF+1,jyF+1)-yFront0{i}(il))/dyv;
        
        vFront{i}(il)=v(ixF,jyF)*FIx(rx1)*FIy(ry1)+v(ixF+1,jyF)*FIx(rx2)*FIy(ry2)+...
            v(ixF,jyF+1)*FIx(rx3)*FIy(ry3)+v(ixF+1,jyF+1)*FIx(rx4)*FIy(ry4);
        vmapFront{i}(il)=vFront{i}(il)/dyv*(1/Ny);
        
        vFront0{i}(il)=v0(ixF,jyF)*FIx(rx1)*FIy(ry1)+v0(ixF+1,jyF)*FIx(rx2)*FIy(ry2)+...
            v0(ixF,jyF+1)*FIx(rx3)*FIy(ry3)+v0(ixF+1,jyF+1)*FIx(rx4)*FIy(ry4);
        vmapFront0{i}(il)=vFront0{i}(il)/dyv*(1/Ny);
        
    end
end
% ------------------------------------------------------------------------
% WeigthSum=....
%     FIx(rx1)*FIy(ry1)+FIx(rx2)*FIy(ry2)+FIx(rx3)*FIy(ry3)+FIx(rx4)*FIy(ry4)
% ------------------------------------------------------------------------

% Move of Front Points ===================================================
for i = 1:2
    for il=2:NFront0{i}+1
    xmapFront{i}(il)=xmapFront0{i}(il)+dt*(umapFront{i}(il));
    ymapFront{i}(il)=ymapFront0{i}(il)+dt*(vmapFront{i}(il));
    % Coordinates in Real Space
    xFront{i}(il)=X1+(X2-X1)*FMx(xmapFront{i}(il)); 
    yFront{i}(il)=Y1+(Y2-Y1)*FMy(ymapFront{i}(il));        
        %     ixF=floor(xmapFront(il)*Nx+1); jyF=floor(ymapFront(il)*Ny+1);
        %     xFront3(il)=(ixF+1-(xmapFront(il)*Nx+1))*xi(ixF)+...
        %         ((xmapFront(il)*Nx+1)-ixF)*xi(ixF);
        %
        %     xFront2(il)=xFront0(il)+dt*uFront(il); yFront2(il)=yFront0(il)+dt*vFront(il);
    end
end
% ========================================================================
xFront{1}(1)=-xFront{1}(3);       yFront{1}(1)=yFront{1}(3);      
xFront{1}(NFront0{1}+2) = xFront{1}(NFront0{1}+1) + ...
    ( xFront{1}(NFront0{1}+1) - xFront{1}(NFront0{1}) );        
yFront{1}(NFront0{1}+2) = yFront{1}(NFront0{1});

uFront{1}(1)=-uFront{1}(3);       vFront{1}(1)=vFront{1}(3);              
uFront{1}(NFront0{1}+2)=uFront{1}(NFront0{1});
vFront{1}(NFront0{1}+2)=vFront{1}(NFront0{1});

xmapFront{1}(1)=-xmapFront{1}(3);       ymapFront{1}(1)=ymapFront{1}(3);      
xmapFront{1}(NFront0{1}+2) = xmapFront{1}(NFront0{1}+1) + ...
    ( xmapFront{1}(NFront0{1}+1) - xmapFront{1}(NFront0{1}) );        
ymapFront{1}(NFront0{1}+2) = ymapFront{1}(NFront0{1});

umapFront{1}(1)=-umapFront{1}(3);       vmapFront{1}(1)=vmapFront{1}(3);              
umapFront{1}(NFront0{1}+2)=umapFront{1}(NFront0{1});
vmapFront{1}(NFront0{1}+2)=vmapFront{1}(NFront0{1});

% ------------------------------------------------------------------------

xFront{2}(1)=-xFront{2}(3);       xFront{2}(NFront0{2}+2)=-xFront{2}(NFront0{2});
yFront{2}(1)=yFront{2}(3);        yFront{2}(NFront0{2}+2)=yFront{2}(NFront0{2});

uFront{2}(1)=-uFront{2}(3);       uFront{2}(NFront0{2}+2)=-uFront{2}(NFront0{2});
vFront{2}(1)=vFront{2}(3);        vFront{2}(NFront0{2}+2)=vFront{2}(NFront0{2});

xmapFront{2}(1)=-xmapFront{2}(3); xmapFront{2}(NFront0{2}+2)=-xmapFront{2}(NFront0{2});
ymapFront{2}(1)=ymapFront{2}(3);  ymapFront{2}(NFront0{2}+2)=ymapFront{2}(NFront0{2});

umapFront{2}(1)=-umapFront{2}(3); umapFront{2}(NFront0{2}+2)=-umapFront{2}(NFront0{2});
vmapFront{2}(1)=vmapFront{2}(3);  vmapFront{2}(NFront0{2}+2)=vmapFront{2}(NFront0{2});
% ========================================================================