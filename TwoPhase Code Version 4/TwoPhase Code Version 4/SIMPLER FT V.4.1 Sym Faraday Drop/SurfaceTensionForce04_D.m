% ========================================================================
% Surface Tension Force
% Weigth Functions For interpolation
%     A -> Bilinear Interpolation (uses the nearest 4points)
%     B -> Peskin (1977) (uses the nearest 4points + 12)
%     C -> Peskin & McQueen(1989) (uses the nearest 4points + 12)
%     D -> Brackbill and Ruppel (1986) (uses the nearest 4points + 12)
% ========================================================================

function [sftx,sfty] = SurfaceTensionForce04_D(Nx,Ny,sigma,...
    xu,yu,xv,yv,NFront,xFront,yFront,xmapFront,ymapFront,Lx,Ly)

% Interpolation function -------------------------------------------------
FIx=@(r) 2/3-(abs(r))^2+1/2*(abs(r))^3; % BrackBill and Ruppel (1986) 
FIy=@(r) 2/3-(abs(r))^2+1/2*(abs(r))^3; % BrackBill and Ruppel (1986) 
FIx2=@(r) 1/6*(2-abs(r))^3; % BrackBill and Ruppel (1986) 
FIy2=@(r) 1/6*(2-abs(r))^3; % BrackBill and Ruppel (1986) 
% -----------------------------------------------------------------------

% Initialization of Surface Tension Matriz ------------------------------
sftx=zeros(Nx+1,Ny+2); sfty=zeros(Nx+2,Ny+1);
% -----------------------------------------------------------------------

% Tangent Vectors to the Interface (eval in l+1/2) ----------------------
ty{1}=zeros(NFront{1}+2,1); tx{1}=zeros(NFront{1}+2,1);
ty{2}=zeros(NFront{2}+2,1); tx{2}=zeros(NFront{2}+2,1);
for i = 1 : 2
    for il=2:NFront{i}+1
        ds=((yFront{i}(il+1)-yFront{i}(il))^2+(xFront{i}(il+1)-xFront{i}(il))^2)^0.5;
        tx{i}(il)=(xFront{i}(il+1)-xFront{i}(il))/ds;
        ty{i}(il)=(yFront{i}(il+1)-yFront{i}(il))/ds;
    end
end
tx{1}(1)=tx{1}(2); ty{1}(1)=-ty{1}(2);
tx{1}(NFront{1}+2)=tx{1}(NFront{1}+1); ty{1}(NFront{1}+2)=-ty{1}(NFront{1}+1);

tx{2}(1)=tx{2}(2); ty{2}(1)=-ty{2}(2);
tx{2}(NFront{2}+2)=tx{2}(NFront{2}+1); ty{2}(NFront{2}+2)=-ty{2}(NFront{2}+1);
% -----------------------------------------------------------------------

% % --- Tangent Vectors to the Interface (eval in l+1/2) ----------------
% ty=zeros(NFront+2,1); tx=zeros(NFront+2,1);
% for il=2:NFront+1
%     ds=((yFront(il+1)-yFront(il))^2+(xFront(il+1)-xFront(il))^2)^0.5;
%     tx(il)=(xFront(il+1)-xFront(il))/ds; 
%     ty(il)=(yFront(il+1)-yFront(il))/ds;
% end
% tx(1)=tx(NFront+1); ty(1)=ty(NFront+1); tx(NFront+2)=tx(2); ty(NFront+2)=ty(2);
% -----------------------------------------------------------------------

% -----------------------------------------------------------------------
STFx{1}=zeros(NFront{1}+2,1); STFy{1}=zeros(NFront{1}+2,1);
STFx{2}=zeros(NFront{2}+2,1); STFy{2}=zeros(NFront{2}+2,1);

for i=1:2
    for il=2:NFront{i}+1
        % Values in the Front of Surface Tension (Its Contain Delta segment)
        STFx{i}(il)=sigma*(tx{i}(il)-tx{i}(il-1)); 
        STFy{i}(il)=sigma*(ty{i}(il)-ty{i}(il-1));
        
        % Transfer of values to the Fixed Grid (a Stagered One!!) ===========
        % Surface Tension in x ----------------------------------------------
        %     ixF=floor(xmapFront(il)*Nx+1); jyF=floor(ymapFront(il)*Ny+1+0.5);
        ixF=floor(xmapFront{i}(il)*Nx+1);   ixF = min(ixF,Nx); 
        jyF=floor(ymapFront{i}(il)*Ny+1+0.5);
        dxu=xu(ixF+1,jyF)-xu(ixF,jyF); dyu=yu(ixF,jyF+1)-yu(ixF,jyF);
        
        rx1=(xFront{i}(il)-xu(ixF,jyF))/dxu; ry1=(yFront{i}(il)-yu(ixF,jyF))/dyu;
        rx2=(xu(ixF+1,jyF)-xFront{i}(il))/dxu; ry2=(yFront{i}(il)-yu(ixF+1,jyF))/dyu;
        rx3=(xFront{i}(il)-xu(ixF,jyF+1))/dxu; ry3=(yu(ixF,jyF+1)-yFront{i}(il))/dyu;
        rx4=(xu(ixF+1,jyF+1)-xFront{i}(il))/dxu ;ry4=(yu(ixF+1,jyF+1)-yFront{i}(il))/dyu;
        
        rAx1=rx1; rAy1=ry1+1; rAx2=rx2; rAy2=ry2+1;
        rAx3=rx1+1; rAy3=ry1; rAx4=rx2+1; rAy4=ry2;
        rAx5=rx3+1; rAy5=ry3; rAx6=rx4+1; rAy6=ry4;
        rAx7=rx3; rAy7=ry3+1; rAx8=rx4; rAy8=ry4+1;
        rAx9=rx1+1; rAy9=ry1+1; rAx10=rx2+1; rAy10=ry2+1;
        rAx11=rx3+1; rAy11=ry3+1; rAx12=rx4+1; rAy12=ry4+1;
        
        sftx(ixF,jyF)=sftx(ixF,jyF)+STFx{i}(il)/(dxu*dyu)*FIx(rx1)*FIy(ry1);
        sftx(ixF+1,jyF)=sftx(ixF+1,jyF)+STFx{i}(il)/(dxu*dyu)*FIx(rx2)*FIy(ry2);
        sftx(ixF,jyF+1)=sftx(ixF,jyF+1)+STFx{i}(il)/(dxu*dyu)*FIx(rx3)*FIy(ry3);
        sftx(ixF+1,jyF+1)=sftx(ixF+1,jyF+1)+STFx{i}(il)/(dxu*dyu)*FIx(rx4)*FIy(ry4);
        
        if ixF == 1
            sftx(ixF,jyF-1)=sftx(ixF,jyF-1)+STFx{i}(il)/(dxu*dyu)*FIx(rAx1)*FIy2(rAy1);
            sftx(ixF+1,jyF-1)=sftx(ixF+1,jyF-1)+STFx{i}(il)/(dxu*dyu)*FIx(rAx2)*FIy2(rAy2);
%             sftx(ixF-1,jyF)=sftx(ixF-1,jyF)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx3)*FIy(rAy3);
            sftx(ixF+2,jyF)=sftx(ixF+2,jyF)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx4)*FIy(rAy4);
%             sftx(ixF-1,jyF+1)=sftx(ixF-1,jyF+1)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx5)*FIy(rAy5);
            sftx(ixF+2,jyF+1)=sftx(ixF+2,jyF+1)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx6)*FIy(rAy6);
            sftx(ixF,jyF+2)=sftx(ixF,jyF+2)+STFx{i}(il)/(dxu*dyu)*FIx(rAx7)*FIy2(rAy7);
            sftx(ixF+1,jyF+2)=sftx(ixF+1,jyF+2)+STFx{i}(il)/(dxu*dyu)*FIx(rAx8)*FIy2(rAy8);
%             sftx(ixF-1,jyF-1)=sftx(ixF-1,jyF-1)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx9)*FIy2(rAy9);
            sftx(ixF+2,jyF-1)=sftx(ixF+2,jyF-1)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx10)*FIy2(rAy10);
%             sftx(ixF-1,jyF+2)=sftx(ixF-1,jyF+2)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx11)*FIy2(rAy11);
            sftx(ixF+2,jyF+2)=sftx(ixF+2,jyF+2)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx12)*FIy2(rAy12);
            
        elseif ixF == Nx
            sftx(ixF,jyF-1)=sftx(ixF,jyF-1)+STFx{i}(il)/(dxu*dyu)*FIx(rAx1)*FIy2(rAy1);
            sftx(ixF+1,jyF-1)=sftx(ixF+1,jyF-1)+STFx{i}(il)/(dxu*dyu)*FIx(rAx2)*FIy2(rAy2);
            sftx(ixF-1,jyF)=sftx(ixF-1,jyF)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx3)*FIy(rAy3);
%             sftx(ixF+2,jyF)=sftx(ixF+2,jyF)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx4)*FIy(rAy4);
            sftx(ixF-1,jyF+1)=sftx(ixF-1,jyF+1)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx5)*FIy(rAy5);
%             sftx(ixF+2,jyF+1)=sftx(ixF+2,jyF+1)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx6)*FIy(rAy6);
            sftx(ixF,jyF+2)=sftx(ixF,jyF+2)+STFx{i}(il)/(dxu*dyu)*FIx(rAx7)*FIy2(rAy7);
            sftx(ixF+1,jyF+2)=sftx(ixF+1,jyF+2)+STFx{i}(il)/(dxu*dyu)*FIx(rAx8)*FIy2(rAy8);
            sftx(ixF-1,jyF-1)=sftx(ixF-1,jyF-1)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx9)*FIy2(rAy9);
%             sftx(ixF+2,jyF-1)=sftx(ixF+2,jyF-1)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx10)*FIy2(rAy10);
            sftx(ixF-1,jyF+2)=sftx(ixF-1,jyF+2)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx11)*FIy2(rAy11);
%             sftx(ixF+2,jyF+2)=sftx(ixF+2,jyF+2)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx12)*FIy2(rAy12);
        else
            sftx(ixF,jyF-1)=sftx(ixF,jyF-1)+STFx{i}(il)/(dxu*dyu)*FIx(rAx1)*FIy2(rAy1);
            sftx(ixF+1,jyF-1)=sftx(ixF+1,jyF-1)+STFx{i}(il)/(dxu*dyu)*FIx(rAx2)*FIy2(rAy2);
            sftx(ixF-1,jyF)=sftx(ixF-1,jyF)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx3)*FIy(rAy3);
            sftx(ixF+2,jyF)=sftx(ixF+2,jyF)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx4)*FIy(rAy4);
            sftx(ixF-1,jyF+1)=sftx(ixF-1,jyF+1)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx5)*FIy(rAy5);
            sftx(ixF+2,jyF+1)=sftx(ixF+2,jyF+1)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx6)*FIy(rAy6);
            sftx(ixF,jyF+2)=sftx(ixF,jyF+2)+STFx{i}(il)/(dxu*dyu)*FIx(rAx7)*FIy2(rAy7);
            sftx(ixF+1,jyF+2)=sftx(ixF+1,jyF+2)+STFx{i}(il)/(dxu*dyu)*FIx(rAx8)*FIy2(rAy8);
            sftx(ixF-1,jyF-1)=sftx(ixF-1,jyF-1)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx9)*FIy2(rAy9);
            sftx(ixF+2,jyF-1)=sftx(ixF+2,jyF-1)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx10)*FIy2(rAy10);
            sftx(ixF-1,jyF+2)=sftx(ixF-1,jyF+2)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx11)*FIy2(rAy11);
            sftx(ixF+2,jyF+2)=sftx(ixF+2,jyF+2)+STFx{i}(il)/(dxu*dyu)*FIx2(rAx12)*FIy2(rAy12);
        end
        % --------------------------------------------------------------------
        
        % Surface Tension in y -----------------------------------------------
        %     ixF=floor(xmapFront(il)*Nx+1+0.5); jyF=floor(ymapFront(il)*Ny+1);
        ixF=floor(xmapFront{i}(il)*Nx+1+0.5);
        jyF=floor(ymapFront{i}(il)*Ny+1);  jyF=jyF-floor(jyF/(Ny+1))*Ny;
        dxv=xv(ixF+1,jyF)-xv(ixF,jyF); dyv=yv(ixF,jyF+1)-yv(ixF,jyF);
        
        rx1=(xFront{i}(il)-xv(ixF,jyF))/dxv; ry1=(yFront{i}(il)-yv(ixF,jyF))/dyv;
        rx2=(xv(ixF+1,jyF)-xFront{i}(il))/dxv; ry2=(yFront{i}(il)-yv(ixF+1,jyF))/dyv;
        rx3=(xFront{i}(il)-xv(ixF,jyF+1))/dxv; ry3=(yv(ixF,jyF+1)-yFront{i}(il))/dyv;
        rx4=(xv(ixF+1,jyF+1)-xFront{i}(il))/dxv; ry4=(yv(ixF+1,jyF+1)-yFront{i}(il))/dyv;
        
        rAx1=rx1; rAy1=ry1+1; rAx2=rx2; rAy2=ry2+1;
        rAx3=rx1+1; rAy3=ry1; rAx4=rx2+1; rAy4=ry2;
        rAx5=rx3+1; rAy5=ry3; rAx6=rx4+1; rAy6=ry4;
        rAx7=rx3; rAy7=ry3+1; rAx8=rx4; rAy8=ry4+1;
        rAx9=rx1+1; rAy9=ry1+1; rAx10=rx2+1; rAy10=ry2+1;
        rAx11=rx3+1; rAy11=ry3+1; rAx12=rx4+1; rAy12=ry4+1;
        
        sfty(ixF,jyF)=sfty(ixF,jyF)+STFy{i}(il)/(dxv*dyv)*FIx(rx1)*FIy(ry1);
        sfty(ixF+1,jyF)=sfty(ixF+1,jyF)+STFy{i}(il)/(dxv*dyv)*FIx(rx2)*FIy(ry2);
        sfty(ixF,jyF+1)=sfty(ixF,jyF+1)+STFy{i}(il)/(dxv*dyv)*FIx(rx3)*FIy(ry3);
        sfty(ixF+1,jyF+1)=sfty(ixF+1,jyF+1)+STFy{i}(il)/(dxv*dyv)*FIx(rx4)*FIy(ry4);
        
        if ixF == 1
            sfty(ixF,jyF-1)=sfty(ixF,jyF-1)+STFy{i}(il)/(dxv*dyv)*FIx(rAx1)*FIy2(rAy1);
            sfty(ixF+1,jyF-1)=sfty(ixF+1,jyF-1)+STFy{i}(il)/(dxv*dyv)*FIx(rAx2)*FIy2(rAy2);
%             sfty(ixF-1,jyF)=sfty(ixF-1,jyF)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx3)*FIy(rAy3);
            sfty(ixF+2,jyF)=sfty(ixF+2,jyF)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx4)*FIy(rAy4);
%             sfty(ixF-1,jyF+1)=sfty(ixF-1,jyF+1)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx5)*FIy(rAy5);
            sfty(ixF+2,jyF+1)=sfty(ixF+2,jyF+1)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx6)*FIy(rAy6);
            sfty(ixF,jyF+2)=sfty(ixF,jyF+2)+STFy{i}(il)/(dxv*dyv)*FIx(rAx7)*FIy2(rAy7);
            sfty(ixF+1,jyF+2)=sfty(ixF+1,jyF+2)+STFy{i}(il)/(dxv*dyv)*FIx(rAx8)*FIy2(rAy8);
%             sfty(ixF-1,jyF-1)=sfty(ixF-1,jyF-1)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx9)*FIy2(rAy9);
            sfty(ixF+2,jyF-1)=sfty(ixF+2,jyF-1)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx10)*FIy2(rAy10);
%             sfty(ixF-1,jyF+2)=sfty(ixF-1,jyF+2)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx11)*FIy2(rAy11);
            sfty(ixF+2,jyF+2)=sfty(ixF+2,jyF+2)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx12)*FIy2(rAy12);
        elseif ixF == Nx + 1
            sfty(ixF,jyF-1)=sfty(ixF,jyF-1)+STFy{i}(il)/(dxv*dyv)*FIx(rAx1)*FIy2(rAy1);
            sfty(ixF+1,jyF-1)=sfty(ixF+1,jyF-1)+STFy{i}(il)/(dxv*dyv)*FIx(rAx2)*FIy2(rAy2);
            sfty(ixF-1,jyF)=sfty(ixF-1,jyF)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx3)*FIy(rAy3);
%             sfty(ixF+2,jyF)=sfty(ixF+2,jyF)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx4)*FIy(rAy4);
            sfty(ixF-1,jyF+1)=sfty(ixF-1,jyF+1)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx5)*FIy(rAy5);
%             sfty(ixF+2,jyF+1)=sfty(ixF+2,jyF+1)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx6)*FIy(rAy6);
            sfty(ixF,jyF+2)=sfty(ixF,jyF+2)+STFy{i}(il)/(dxv*dyv)*FIx(rAx7)*FIy2(rAy7);
            sfty(ixF+1,jyF+2)=sfty(ixF+1,jyF+2)+STFy{i}(il)/(dxv*dyv)*FIx(rAx8)*FIy2(rAy8);
            sfty(ixF-1,jyF-1)=sfty(ixF-1,jyF-1)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx9)*FIy2(rAy9);
%             sfty(ixF+2,jyF-1)=sfty(ixF+2,jyF-1)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx10)*FIy2(rAy10);
            sfty(ixF-1,jyF+2)=sfty(ixF-1,jyF+2)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx11)*FIy2(rAy11);
%             sfty(ixF+2,jyF+2)=sfty(ixF+2,jyF+2)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx12)*FIy2(rAy12);
        else
            sfty(ixF,jyF-1)=sfty(ixF,jyF-1)+STFy{i}(il)/(dxv*dyv)*FIx(rAx1)*FIy2(rAy1);
            sfty(ixF+1,jyF-1)=sfty(ixF+1,jyF-1)+STFy{i}(il)/(dxv*dyv)*FIx(rAx2)*FIy2(rAy2);
            sfty(ixF-1,jyF)=sfty(ixF-1,jyF)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx3)*FIy(rAy3);
            sfty(ixF+2,jyF)=sfty(ixF+2,jyF)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx4)*FIy(rAy4);
            sfty(ixF-1,jyF+1)=sfty(ixF-1,jyF+1)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx5)*FIy(rAy5);
            sfty(ixF+2,jyF+1)=sfty(ixF+2,jyF+1)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx6)*FIy(rAy6);
            sfty(ixF,jyF+2)=sfty(ixF,jyF+2)+STFy{i}(il)/(dxv*dyv)*FIx(rAx7)*FIy2(rAy7);
            sfty(ixF+1,jyF+2)=sfty(ixF+1,jyF+2)+STFy{i}(il)/(dxv*dyv)*FIx(rAx8)*FIy2(rAy8);
            sfty(ixF-1,jyF-1)=sfty(ixF-1,jyF-1)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx9)*FIy2(rAy9);
            sfty(ixF+2,jyF-1)=sfty(ixF+2,jyF-1)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx10)*FIy2(rAy10);
            sfty(ixF-1,jyF+2)=sfty(ixF-1,jyF+2)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx11)*FIy2(rAy11);
            sfty(ixF+2,jyF+2)=sfty(ixF+2,jyF+2)+STFy{i}(il)/(dxv*dyv)*FIx2(rAx12)*FIy2(rAy12);
        end
        % --------------------------------------------------------------------
        
    end
end
% -----------------------------------------------------------------------
% =======================================================================

% %%%%%%%%%%%%%%%%%%%%%%%%%%%% % Plots
% ftxx(1:Nx+1,1:Ny+1)=0.5*(sftx(1:Nx+1,2:Ny+2)+sftx(1:Nx+1,1:Ny+1));
% ftyy(1:Nx+1,1:Ny+1)=0.5*(sfty(2:Nx+2,1:Ny+1)+sfty(1:Nx+1,1:Ny+1));
% figure; quiver(xi,yj,flipud(rot90(ftxx)),flipud(rot90(ftyy))); axis equal