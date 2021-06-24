% ========================================================================
% Constructing the density field
% Weigth Functions
%     A -> Bilinear Interpolation (uses the nearest 4points)
%     B -> Peskin (1977) (uses the nearest 4points + 12)
%     C -> Peskin & McQueen(1989) (uses the nearest 4points + 12)
%     D -> Brackbill and Ruppel (1986) (uses the nearest 4points + 12)
% ========================================================================

function [MarkF,MarkFdx,MarkFdy]...
    =ConstructMarkerFunction03_A(MarkF0,MF01,MF02,Nx,Ny,...
    xu,yu,xv,yv,dx_PE,dy_PN,dx_WP,dy_SP,...
    NFront,xFront,yFront,xmapFront,ymapFront)

% Interpolation function -------------------------------------------------
FIx=@(r) 1-abs(r); % lineal Interpolation 
FIy=@(r) 1-abs(r); % lineal Interpolation
% ------------------------------------------------------------------------

% =======================================================================
% Gradient Distribution -------------------------------------------------
MarkFdx=zeros(Nx+1,Ny+2); MarkFdy=zeros(Nx+2,Ny+1);
Flagx=zeros(Nx+2,Ny+2); Flagy=zeros(Nx+2,Ny+2);
for il=2:NFront+1
    % Values in the Front (with Normal Vectors) -------------------------
    nfx=-0.5*(yFront(il+1)-yFront(il-1))*(MF01-MF02);
    nfy=0.5*(xFront(il+1)-xFront(il-1))*(MF01-MF02);
    % -------------------------------------------------------------------
    
    % Gradients in the x direction (Located in laterals C.V Faces)
    ixF=floor(xmapFront(il)*Nx+1); jyF=floor(ymapFront(il)*Ny+1+0.5);
    dxu=xu(ixF+1,jyF)-xu(ixF,jyF); dyu=yu(ixF,jyF+1)-yu(ixF,jyF);
    
    rx1=(xFront(il)-xu(ixF,jyF))/dxu; ry1=(yFront(il)-yu(ixF,jyF))/dyu;
    rx2=(xu(ixF+1,jyF)-xFront(il))/dxu; ry2=(yFront(il)-yu(ixF+1,jyF))/dyu;
    rx3=(xFront(il)-xu(ixF,jyF+1))/dxu; ry3=(yu(ixF,jyF+1)-yFront(il))/dyu;
    rx4=(xu(ixF+1,jyF+1)-xFront(il))/dxu ;ry4=(yu(ixF+1,jyF+1)-yFront(il))/dyu;
    
    MarkFdx(ixF,jyF)=MarkFdx(ixF,jyF)+nfx/(dxu*dyu)*FIx(rx1)*FIy(ry1);
    MarkFdx(ixF+1,jyF)=MarkFdx(ixF+1,jyF)+nfx/(dxu*dyu)*FIx(rx2)*FIy(ry2);
    MarkFdx(ixF,jyF+1)=MarkFdx(ixF,jyF+1)+nfx/(dxu*dyu)*FIx(rx3)*FIy(ry3);
    MarkFdx(ixF+1,jyF+1)=MarkFdx(ixF+1,jyF+1)+nfx/(dxu*dyu)*FIx(rx4)*FIy(ry4);
    % -------------------------------------------------------------------
    
    % Flags x -----------------------------------------------------------
    Flagx(ixF,jyF)=1; Flagx(ixF+1,jyF)=1; Flagx(ixF+2,jyF)=1; 
    Flagx(ixF,jyF+1)=1; Flagx(ixF+1,jyF+1)=1; Flagx(ixF+2,jyF+1)=1; 
    % -------------------------------------------------------------------
    
    % Gradients in the y direction (Located in top and bottom C.V Faces)
    ixF=floor(xmapFront(il)*Nx+1+0.5); jyF=floor(ymapFront(il)*Ny+1);
    dxv=xv(ixF+1,jyF)-xv(ixF,jyF); dyv=yv(ixF,jyF+1)-yv(ixF,jyF);
    
    rx1=(xFront(il)-xv(ixF,jyF))/dxv; ry1=(yFront(il)-yv(ixF,jyF))/dyv;
    rx2=(xv(ixF+1,jyF)-xFront(il))/dxv; ry2=(yFront(il)-yv(ixF+1,jyF))/dyv;
    rx3=(xFront(il)-xv(ixF,jyF+1))/dxv; ry3=(yv(ixF,jyF+1)-yFront(il))/dyv;
    rx4=(xv(ixF+1,jyF+1)-xFront(il))/dxv; ry4=(yv(ixF+1,jyF+1)-yFront(il))/dyv;
    
    MarkFdy(ixF,jyF)=MarkFdy(ixF,jyF)+nfy/(dxv*dyv)*FIx(rx1)*FIy(ry1);
    MarkFdy(ixF+1,jyF)=MarkFdy(ixF+1,jyF)+nfy/(dxv*dyv)*FIx(rx2)*FIy(ry2);
    MarkFdy(ixF,jyF+1)=MarkFdy(ixF,jyF+1)+nfy/(dxv*dyv)*FIx(rx3)*FIy(ry3);
    MarkFdy(ixF+1,jyF+1)=MarkFdy(ixF+1,jyF+1)+nfy/(dxv*dyv)*FIx(rx4)*FIy(ry4); 
    % -------------------------------------------------------------------
    
    % Flags y -----------------------------------------------------------
    Flagy(ixF,jyF)=1; Flagy(ixF,jyF+1)=1; Flagy(ixF,jyF+2)=1;
    Flagy(ixF+1,jyF)=1; Flagy(ixF+1,jyF+1)=1; Flagy(ixF+1,jyF+2)=1;
    % -------------------------------------------------------------------

end
% -------------------------------------------------------------------
% WeigthSum=....
%     FIx(rx1)*FIy(ry1)+FIx(rx2)*FIy(ry2)+FIx(rx3)*FIy(ry3)+FIx(rx4)*FIy(ry4)
% -------------------------------------------------------------------

% ========================================================================
% Construction of Density -----------------------------------------------
%  Location of The nearest grid Points ------------ 
% ixmin=min(floor(xmapFront*Nx+1+0.5))-3;  jymin=min(floor(ymapFront*Ny+1+0.5))-3;
% ixmax=max(floor(xmapFront*Nx+1+0.5)+1)+3; jymax=max(floor(ymapFront*Ny+1+0.5)+1)+3;

ixmin=2;  jymin=2; ixmax=Nx+1; jymax=Ny+1;
% Iteration for the Marker Function --------------------------------------
Tol=1e-6; Error=1; nitera=0; nmax=500; % ome=1.0;
MarkF=round(MarkF0); 
% MarkF(ixmin-4,:)=0; MarkF(ixmax+4,:)=0; MarkF(:,jymin-4)=0; MarkF(:,jymax+4)=0;
while (Error>=Tol)&&(nitera<=nmax)
    roLast=MarkF; nitera=nitera+1;
    for i=ixmin:ixmax
        for j=jymin:jymax
            MarkF(i,j)=max(Flagx(i,j),Flagy(i,j))*...
                (0.25*(MarkF(i+1,j)+MarkF(i-1,j)+MarkF(i,j+1)+MarkF(i,j-1)+...
                dx_WP(i,j)*MarkFdx(i-1,j)+dy_SP(i,j)*MarkFdy(i,j-1)+...
                -dx_PE(i,j)*MarkFdx(i,j)-dy_PN(i,j)*MarkFdy(i,j)))+...
                (1-max(Flagx(i,j),Flagy(i,j)))*MarkF(i,j);
        end
    end
    Error=max(max(abs((MarkF-roLast)./MarkF)));
end
% disp('Construccion de la densidad')
info=strcat('#It. en la construccion de la MarkerF = ',num2str(nitera),'_ Error = ',num2str(Error));
disp(info)
% -----------------------------------------------------------------------
% ========================================================================
% figure; surf(xP,yP,ro); figure; surf(xP,yP,ro0);
