% ========================================================================
% Surface Tension Force
% Weigth Functions For interpolation
%     A -> Bilinear Interpolation (uses the nearest 4points)
%     B -> Peskin (1977) (uses the nearest 4points + 12)
%     C -> Peskin & McQueen(1989) (uses the nearest 4points + 12)
%     D -> Brackbill and Ruppel (1986) (uses the nearest 4points + 12)
% ========================================================================

function [sftx,sfty] = SurfaceTensionForce04_A(Nx,Ny,sigma,...
    xu,yu,xv,yv,NFront,xFront,yFront,xmapFront,ymapFront)

% Interpolation function ------------------------------------------------
FIx=@(r) 1-abs(r); % lineal Interpolation 
FIy=@(r) 1-abs(r); % lineal Interpolation
% -----------------------------------------------------------------------

% Initialization of Surface Tension Matriz ------------------------------
sftx=zeros(Nx+1,Ny+2); sfty=zeros(Nx+2,Ny+1);
% -----------------------------------------------------------------------

% Tangent Vectors to the Interface (eval in l+1/2) ----------------------
ty=zeros(NFront+2,1); tx=zeros(NFront+2,1);
for il=2:NFront+1
    ds=((yFront(il+1)-yFront(il))^2+(xFront(il+1)-xFront(il))^2)^0.5;
    tx(il)=(xFront(il+1)-xFront(il))/ds; 
    ty(il)=(yFront(il+1)-yFront(il))/ds;
end
tx(1)=tx(NFront+1); ty(1)=ty(NFront+1); tx(NFront+2)=tx(2); ty(NFront+2)=ty(2);
% -----------------------------------------------------------------------

% -----------------------------------------------------------------------
STFx=zeros(NFront+2,1); STFy=zeros(NFront+2,1);
for il=2:NFront+1
    % Values in the Front of Surface Tension (Its Contain Delta segment)
    STFx(il)=sigma*(tx(il)-tx(il-1)); STFy(il)=sigma*(ty(il)-ty(il-1));
    % -------------------------------------------------------------------
    
    % Transfer of values to the Fixed Grid (a Stagered One!!)
    % Surface Tension in x ----------------------------------------------
    ixF=floor(xmapFront(il)*Nx+1); jyF=floor(ymapFront(il)*Ny+1+0.5);
    dxu=xu(ixF+1,jyF)-xu(ixF,jyF); dyu=yu(ixF,jyF+1)-yu(ixF,jyF);
    
    rx1=(xFront(il)-xu(ixF,jyF))/dxu; ry1=(yFront(il)-yu(ixF,jyF))/dyu;
    rx2=(xu(ixF+1,jyF)-xFront(il))/dxu; ry2=(yFront(il)-yu(ixF+1,jyF))/dyu;
    rx3=(xFront(il)-xu(ixF,jyF+1))/dxu; ry3=(yu(ixF,jyF+1)-yFront(il))/dyu;
    rx4=(xu(ixF+1,jyF+1)-xFront(il))/dxu ;ry4=(yu(ixF+1,jyF+1)-yFront(il))/dyu;
    
    sftx(ixF,jyF)=sftx(ixF,jyF)+STFx(il)/(dxu*dyu)*FIx(rx1)*FIy(ry1);
    sftx(ixF+1,jyF)=sftx(ixF+1,jyF)+STFx(il)/(dxu*dyu)*FIx(rx2)*FIy(ry2);
    sftx(ixF,jyF+1)=sftx(ixF,jyF+1)+STFx(il)/(dxu*dyu)*FIx(rx3)*FIy(ry3);
    sftx(ixF+1,jyF+1)=sftx(ixF+1,jyF+1)+STFx(il)/(dxu*dyu)*FIx(rx4)*FIy(ry4);
    % -------------------------------------------------------------------
    
    % Surface Tension in y ----------------------------------------------
    ixF=floor(xmapFront(il)*Nx+1+0.5); jyF=floor(ymapFront(il)*Ny+1);
    dxv=xv(ixF+1,jyF)-xv(ixF,jyF); dyv=yv(ixF,jyF+1)-yv(ixF,jyF);
    
    rx1=(xFront(il)-xv(ixF,jyF))/dxv; ry1=(yFront(il)-yv(ixF,jyF))/dyv;
    rx2=(xv(ixF+1,jyF)-xFront(il))/dxv; ry2=(yFront(il)-yv(ixF+1,jyF))/dyv;
    rx3=(xFront(il)-xv(ixF,jyF+1))/dxv; ry3=(yv(ixF,jyF+1)-yFront(il))/dyv;
    rx4=(xv(ixF+1,jyF+1)-xFront(il))/dxv; ry4=(yv(ixF+1,jyF+1)-yFront(il))/dyv;
    
    sfty(ixF,jyF)=sfty(ixF,jyF)+STFy(il)/(dxv*dyv)*FIx(rx1)*FIy(ry1);
    sfty(ixF+1,jyF)=sfty(ixF+1,jyF)+STFy(il)/(dxv*dyv)*FIx(rx2)*FIy(ry2);
    sfty(ixF,jyF+1)=sfty(ixF,jyF+1)+STFy(il)/(dxv*dyv)*FIx(rx3)*FIy(ry3);
    sfty(ixF+1,jyF+1)=sfty(ixF+1,jyF+1)+STFy(il)/(dxv*dyv)*FIx(rx4)*FIy(ry4); 
    % -------------------------------------------------------------------
end
% =======================================================================

% %%%%%%%%%%%%%%%%%%%%%%%%%%%% % Plots
% ftxx(1:Nx+1,1:Ny+1)=0.5*(sftx(1:Nx+1,2:Ny+2)+sftx(1:Nx+1,1:Ny+1));
% ftyy(1:Nx+1,1:Ny+1)=0.5*(sfty(2:Nx+2,1:Ny+1)+sfty(1:Nx+1,1:Ny+1));
% figure; quiver(xi,yj,flipud(rot90(ftxx)),flipud(rot90(ftyy))); axis equal