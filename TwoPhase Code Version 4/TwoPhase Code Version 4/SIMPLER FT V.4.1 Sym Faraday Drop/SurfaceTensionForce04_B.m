% ========================================================================
% Surface Tension Force
% Weigth Functions For interpolation
%     A -> Bilinear Interpolation (uses the nearest 4points)
%     B -> Peskin (1977) (uses the nearest 4points + 12)
%     C -> Peskin & McQueen(1989) (uses the nearest 4points + 12)
%     D -> Brackbill and Ruppel (1986) (uses the nearest 4points + 12)
% ========================================================================

function [sftx,sfty] = SurfaceTensionForce04_B(Nx,Ny,sigma,...
    xu,yu,xv,yv,NFront,xFront,yFront,xmapFront,ymapFront,Lx,Ly)

% Interpolation function ------------------------------------------------
FIx=@(r) 0.25*(1+cos(pi*r*0.5)); % Peskin (1977) Interpolation 
FIy=@(r) 0.25*(1+cos(pi*r*0.5)); % Peskin (1977) Interpolation
% -----------------------------------------------------------------------

% Initialization of Surface Tension Matriz ------------------------------
sftx=zeros(Nx+1,Ny+2); sfty=zeros(Nx+2,Ny+1);
% -----------------------------------------------------------------------

%  Tangent Vectors to the Interface (eval in l+1/2) ---------------------
ty=zeros(NFront+2,1); tx=zeros(NFront+2,1);
for il=2:NFront+1
    [dx,Ix]=min([abs(xFront(il+1)-xFront(il)),...
        abs(xFront(il+1)-Lx-xFront(il)),...
        abs(xFront(il+1)+Lx-xFront(il))]);
    signx=sign([(xFront(il+1)-xFront(il)),...
        (xFront(il+1)-Lx-xFront(il)),...
        (xFront(il+1)+Lx-xFront(il))]);
    
    [dy,Iy]=min([abs(yFront(il+1)-yFront(il)),...
        abs(yFront(il+1)-Ly-yFront(il)),...
        abs(yFront(il+1)+Ly-yFront(il))]);
    signy=sign([(yFront(il+1)-yFront(il)),...
        (yFront(il+1)-Ly-yFront(il)),...
        (yFront(il+1)+Ly-yFront(il))]);

    ds=(dy^2+dx^2)^0.5;
    tx(il)=signx(Ix)*dx/ds; ty(il)=signy(Iy)*dy/ds;
end
tx(1)=tx(NFront+1); ty(1)=ty(NFront+1); tx(NFront+2)=tx(2); ty(NFront+2)=ty(2);
% -----------------------------------------------------------------------
% 
% % Tangent Vectors to the Interface (eval in l+1/2) --------------------
% ty=zeros(NFront+2,1); tx=zeros(NFront+2,1);
% for il=2:NFront+1
%     ds=((yFront(il+1)-yFront(il))^2+(xFront(il+1)-xFront(il))^2)^0.5;
%     tx(il)=(xFront(il+1)-xFront(il))/ds; 
%     ty(il)=(yFront(il+1)-yFront(il))/ds;
% end
% tx(1)=tx(NFront+1); ty(1)=ty(NFront+1); tx(NFront+2)=tx(2); ty(NFront+2)=ty(2);
% -----------------------------------------------------------------------

% -----------------------------------------------------------------------
STFx=zeros(NFront+2,1); STFy=zeros(NFront+2,1);
for il=2:NFront+1
    % Values in the Front of Surface Tension (Its Contain Delta segment)
    STFx(il)=sigma*(tx(il)-tx(il-1)); STFy(il)=sigma*(ty(il)-ty(il-1));
    
    % Transfer of values to the Fixed Grid (a Stagered One!!) -----------
    % Surface Tension in x ----------------------------------------------
%     ixF=floor(xmapFront(il)*Nx+1); jyF=floor(ymapFront(il)*Ny+1+0.5);
    ixF=floor(xmapFront(il)*Nx+1); ixF=ixF-floor(ixF/(Nx+1))*Nx;
    jyF=floor(ymapFront(il)*Ny+1+0.5); 
    dxu=xu(ixF+1,jyF)-xu(ixF,jyF); dyu=yu(ixF,jyF+1)-yu(ixF,jyF);
    
    rx1=min([abs(xFront(il)-xu(ixF,jyF)),abs(xFront(il)-xu(ixF,jyF)-Lx),...
        abs(xFront(il)-xu(ixF,jyF)+Lx)])/dxu;
    rx2=min([abs(xu(ixF+1,jyF)-xFront(il)),abs(xu(ixF+1,jyF)-xFront(il)-Lx),...
        abs(xu(ixF+1,jyF)-xFront(il)+Lx)])/dxu;
    ry1=(yFront(il)-yu(ixF,jyF))/dyu; ry2=(yFront(il)-yu(ixF+1,jyF))/dyu;
    rx3=min([abs(xFront(il)-xu(ixF,jyF+1)),abs(xFront(il)-xu(ixF,jyF+1)-Lx),...
        abs(xFront(il)-xu(ixF,jyF+1)+Lx)])/dxu;
    rx4=min([abs(xu(ixF+1,jyF+1)-xFront(il)),abs(xu(ixF+1,jyF+1)-xFront(il)-Lx),...
        abs(xu(ixF+1,jyF+1)-xFront(il)+Lx)])/dxu;
    ry3=(yu(ixF,jyF+1)-yFront(il))/dyu; ry4=(yu(ixF+1,jyF+1)-yFront(il))/dyu;
    
    rAx1=rx1; rAy1=ry1+1; rAx2=rx2; rAy2=ry2+1;
    rAx3=rx1+1; rAy3=ry1; rAx4=rx2+1; rAy4=ry2;
    rAx5=rx3+1; rAy5=ry3; rAx6=rx4+1; rAy6=ry4;
    rAx7=rx3; rAy7=ry3+1; rAx8=rx4; rAy8=ry4+1;
    rAx9=rx1+1; rAy9=ry1+1; rAx10=rx2+1; rAy10=ry2+1;
    rAx11=rx3+1; rAy11=ry3+1; rAx12=rx4+1; rAy12=ry4+1;
    
    sftx(ixF,jyF)=sftx(ixF,jyF)+STFx(il)/(dxu*dyu)*FIx(rx1)*FIy(ry1);
    sftx(ixF+1,jyF)=sftx(ixF+1,jyF)+STFx(il)/(dxu*dyu)*FIx(rx2)*FIy(ry2);
    sftx(ixF,jyF+1)=sftx(ixF,jyF+1)+STFx(il)/(dxu*dyu)*FIx(rx3)*FIy(ry3);
    sftx(ixF+1,jyF+1)=sftx(ixF+1,jyF+1)+STFx(il)/(dxu*dyu)*FIx(rx4)*FIy(ry4);
    
    ixFA1=max( ixF-2-floor((ixF-2)/(Nx))*(Nx+1), ixF-1 ); % for ixF-1
    ixFA2=ixF+2-floor((ixF+2)/(Nx+1))*Nx; % for ixF+2
    sftx(ixF,jyF-1)=sftx(ixF,jyF-1)+STFx(il)/(dxu*dyu)*FIx(rAx1)*FIy(rAy1);
    sftx(ixF+1,jyF-1)=sftx(ixF+1,jyF-1)+STFx(il)/(dxu*dyu)*FIx(rAx2)*FIy(rAy2);
    sftx(ixFA1,jyF)=sftx(ixFA1,jyF)+STFx(il)/(dxu*dyu)*FIx(rAx3)*FIy(rAy3);
    sftx(ixFA2,jyF)=sftx(ixFA2,jyF)+STFx(il)/(dxu*dyu)*FIx(rAx4)*FIy(rAy4);
    sftx(ixFA1,jyF+1)=sftx(ixFA1,jyF+1)+STFx(il)/(dxu*dyu)*FIx(rAx5)*FIy(rAy5);
    sftx(ixFA2,jyF+1)=sftx(ixFA2,jyF+1)+STFx(il)/(dxu*dyu)*FIx(rAx6)*FIy(rAy6);
    sftx(ixF,jyF+2)=sftx(ixF,jyF+2)+STFx(il)/(dxu*dyu)*FIx(rAx7)*FIy(rAy7);
    sftx(ixF+1,jyF+2)=sftx(ixF+1,jyF+2)+STFx(il)/(dxu*dyu)*FIx(rAx8)*FIy(rAy8);
    sftx(ixFA1,jyF-1)=sftx(ixFA1,jyF-1)+STFx(il)/(dxu*dyu)*FIx(rAx9)*FIy(rAy9);
    sftx(ixFA2,jyF-1)=sftx(ixFA2,jyF-1)+STFx(il)/(dxu*dyu)*FIx(rAx10)*FIy(rAy10);
    sftx(ixFA1,jyF+2)=sftx(ixFA1,jyF+2)+STFx(il)/(dxu*dyu)*FIx(rAx11)*FIy(rAy11);
    sftx(ixFA2,jyF+2)=sftx(ixFA2,jyF+2)+STFx(il)/(dxu*dyu)*FIx(rAx12)*FIy(rAy12);
    % --------------------------------------------------------------------
    
    % Surface Tension in y -----------------------------------------------
%     ixF=floor(xmapFront(il)*Nx+1+0.5); jyF=floor(ymapFront(il)*Ny+1);
    ixF=floor(xmapFront(il)*Nx+1+0.5); 
    jyF=floor(ymapFront(il)*Ny+1);  jyF=jyF-floor(jyF/(Ny+1))*Ny;
    dxv=xv(ixF+1,jyF)-xv(ixF,jyF); dyv=yv(ixF,jyF+1)-yv(ixF,jyF);
    
    rx1=(xFront(il)-xv(ixF,jyF))/dxv; ry1=(yFront(il)-yv(ixF,jyF))/dyv;
    rx2=(xv(ixF+1,jyF)-xFront(il))/dxv; ry2=(yFront(il)-yv(ixF+1,jyF))/dyv;
    rx3=(xFront(il)-xv(ixF,jyF+1))/dxv; ry3=(yv(ixF,jyF+1)-yFront(il))/dyv;
    rx4=(xv(ixF+1,jyF+1)-xFront(il))/dxv; ry4=(yv(ixF+1,jyF+1)-yFront(il))/dyv;
    rAx1=rx1; rAy1=ry1+1; rAx2=rx2; rAy2=ry2+1;
    rAx3=rx1+1; rAy3=ry1; rAx4=rx2+1; rAy4=ry2;
    rAx5=rx3+1; rAy5=ry3; rAx6=rx4+1; rAy6=ry4;
    rAx7=rx3; rAy7=ry3+1; rAx8=rx4; rAy8=ry4+1;
    rAx9=rx1+1; rAy9=ry1+1; rAx10=rx2+1; rAy10=ry2+1;
    rAx11=rx3+1; rAy11=ry3+1; rAx12=rx4+1; rAy12=ry4+1;
    
    sfty(ixF,jyF)=sfty(ixF,jyF)+STFy(il)/(dxv*dyv)*FIx(rx1)*FIy(ry1);
    sfty(ixF+1,jyF)=sfty(ixF+1,jyF)+STFy(il)/(dxv*dyv)*FIx(rx2)*FIy(ry2);
    sfty(ixF,jyF+1)=sfty(ixF,jyF+1)+STFy(il)/(dxv*dyv)*FIx(rx3)*FIy(ry3);
    sfty(ixF+1,jyF+1)=sfty(ixF+1,jyF+1)+STFy(il)/(dxv*dyv)*FIx(rx4)*FIy(ry4);
    
    ixFA1=max( ixF-2-floor((ixF-2)/(Nx+1))*Nx+1, ixF-1 ); % for ixF-1
    ixFA2=ixF+2-floor((ixF+2)/(Nx+3))*(Nx); % for ixF+2
    sfty(ixF,jyF-1)=sfty(ixF,jyF-1)+STFy(il)/(dxv*dyv)*FIx(rAx1)*FIy(rAy1);
    sfty(ixF+1,jyF-1)=sfty(ixF+1,jyF-1)+STFy(il)/(dxv*dyv)*FIx(rAx2)*FIy(rAy2);
    sfty(ixFA1,jyF)=sfty(ixFA1,jyF)+STFy(il)/(dxv*dyv)*FIx(rAx3)*FIy(rAy3);
    sfty(ixFA2,jyF)=sfty(ixFA2,jyF)+STFy(il)/(dxv*dyv)*FIx(rAx4)*FIy(rAy4);
    sfty(ixFA1,jyF+1)=sfty(ixFA1,jyF+1)+STFy(il)/(dxv*dyv)*FIx(rAx5)*FIy(rAy5);
    sfty(ixFA2,jyF+1)=sfty(ixFA2,jyF+1)+STFy(il)/(dxv*dyv)*FIx(rAx6)*FIy(rAy6);
    sfty(ixF,jyF+2)=sfty(ixF,jyF+2)+STFy(il)/(dxv*dyv)*FIx(rAx7)*FIy(rAy7);
    sfty(ixF+1,jyF+2)=sfty(ixF+1,jyF+2)+STFy(il)/(dxv*dyv)*FIx(rAx8)*FIy(rAy8);
    sfty(ixFA1,jyF-1)=sfty(ixFA1,jyF-1)+STFy(il)/(dxv*dyv)*FIx(rAx9)*FIy(rAy9);
    sfty(ixFA2,jyF-1)=sfty(ixFA2,jyF-1)+STFy(il)/(dxv*dyv)*FIx(rAx10)*FIy(rAy10);
    sfty(ixFA1,jyF+2)=sfty(ixFA1,jyF+2)+STFy(il)/(dxv*dyv)*FIx(rAx11)*FIy(rAy11);
    sfty(ixFA2,jyF+2)=sfty(ixFA2,jyF+2)+STFy(il)/(dxv*dyv)*FIx(rAx12)*FIy(rAy12);
    % --------------------------------------------------------------------
    
end
% -----------------------------------------------------------------------

% Simetria ---------------------------------------------------------------
% for j=1:Ny+2
%     [sftxSime,ISime]=max([abs(sftx(1,j)),abs(sftx(Nx+1,j))]);
%     signSime=sign([sftx(1,j),sftx(Nx+1,j)]);
%     sftx(1,j)=signSime(ISime)*sftxSime;
%     sftx(Nx+1,j)=signSime(ISime)*sftxSime;
% end
% for j=1:Ny+1
%     [sfty1Sime,ISime1]=max([abs(sfty(Nx+1,j)),abs(sfty(1,j))]);
%     signSime1=sign([sfty(Nx+1,j),sfty(1,j)]);
%     sfty(1,j)=signSime1(ISime1)*sfty1Sime;
%     sfty(Nx+1,j)=signSime1(ISime1)*sfty1Sime;
%     
%     [sfty2Sime,ISime2]=max([abs(sfty(Nx+2,j)),abs(sfty(2,j))]);
%     signSime2=sign([sfty(Nx+2,j),sfty(2,j)]);
%     sfty(Nx+2,j)=signSime2(ISime2)*sfty2Sime;
%     sfty(2,j)=signSime2(ISime2)*sfty2Sime;
% end
sftx(1,:)=sftx(1,:)+sftx(Nx+1,:); sftx(Nx+1,:)=sftx(1,:);
sfty(1,:)=sfty(Nx+1,:)+sfty(1,:); sfty(Nx+1,:)=sfty(1,:);
sfty(Nx+2,:)=sfty(Nx+2,:)+sfty(2,:); sfty(2,:)=sfty(Nx+2,:);
% -----------------------------------------------------------------------

% %%%%%%%%%%%%%%%%%%%%%%%%%%%% % Plots
% ftxx(1:Nx+1,1:Ny+1)=0.5*(sftx(1:Nx+1,2:Ny+2)+sftx(1:Nx+1,1:Ny+1));
% ftyy(1:Nx+1,1:Ny+1)=0.5*(sfty(2:Nx+2,1:Ny+1)+sfty(1:Nx+1,1:Ny+1));
% figure; quiver(xi,yj,flipud(rot90(ftxx)),flipud(rot90(ftyy))); axis equal