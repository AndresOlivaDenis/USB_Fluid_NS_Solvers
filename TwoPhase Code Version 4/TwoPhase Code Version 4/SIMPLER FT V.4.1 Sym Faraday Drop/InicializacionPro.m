% ======================================================================= %
% Inicializacion de propiedades y variables
% =======================================================================

% Initial Marker Function -----------------------------------------------
% Pool Surface -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
MF01=1; MF02=0; MarkF=MF02*ones(Nx+2,Ny+2);
for i=2:Nx+1,
    for j=1:Ny+2
        if yP(i,j) >= ysurf+amp*cos(k*xP(i,j) + Lx*k)
            MarkF(i,j)=MF01; end
    end,
end
% -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
% Bubble  .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
for i=2:Nx+1,
    for j=2:Ny+1
        if (((xP(i,j)-xdrop)^2+(yP(i,j)-ydrop)^2)<=rdrop^2)
            MarkF(i,j)=MF02; end
    end,
end
% -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
MarkF(1,:)=MarkF(2,:);      MarkF(Nx+2,:)=MarkF(Nx+1,:);
% -----------------------------------------------------------------------

% Imposicion de las propiedades Fisicas ---------------------------------
% Propiedades y terminos fuentes constantes:
ro=RO2.*MarkF+(1-MarkF).*RO1; % miu=MIU2*ones(Nx+2,Ny+2);
miu=MIU2.*MarkF+(1-MarkF).*MIU1;
b=B*ones(Nx+2,Ny+2);

ro0=ro; miu0=miu;
% -----------------------------------------------------------------------

% Condicion Inicial -----------------------------------------------------
% Velocidades y Presion
u0=U0*ones(Nx+1,Ny+2); v0=V0*ones(Nx+2,Ny+1);
u=u0; v=v0; p=zeros(Nx+2,Ny+2); pc=zeros(Nx+2,Ny+2);
us=u0; vs=v0; ut=u0; vt=v0;
% -----------------------------------------------------------------------

% Termino fuente ----------------------------------------------------
Su=gx*ones(Nx+1,Ny+2); Sv=gy*ones(Nx+2,Ny+1);
% -------------------------------------------------------------------

% =======================================================================
% Set and Initialization of the Front --------------------------------
% Pool Front -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
NFront{1}=Nx*4; xFront{1}=zeros(1,NFront{1}+2); yFront{1}=zeros(1,NFront{1}+2);
uFront{1}=zeros(1,NFront{1}+2); vFront{1}=zeros(1,NFront{1}+2);
rodxFront{1}=zeros(1,NFront{1}+2); rodyFront{1}=zeros(1,NFront{1}+2);

% for il=2:NFront+1
%     xFront(il)=(il-2)*(Lx)/(NFront);
%     yFront(il)=ysurf+amp*cos(ome*xFront(il));
% end
for il=2:NFront{1}+1
    xFront{1}(il)=(il-2)*(Lx)/(NFront{1}-1);
    yFront{1}(il)=ysurf+amp*cos(k*xFront{1}(il)+ Lx*k);
end
xFront{1}(1)=-xFront{1}(3);  
xFront{1}(NFront{1}+2)=xFront{1}(NFront{1}+1) + ...
    ( xFront{1}(NFront{1}+1) - xFront{1}(NFront{1}) ) ;
yFront{1}(1)=yFront{1}(3);  yFront{1}(NFront{1}+2)=yFront{1}(NFront{1}+1);

% For de Reestructure of the Front
FrontDistMax{1}=2.50*((xFront{1}(2)-xFront{1}(1))^2+(yFront{1}(2)-yFront{1}(1))^2)^0.5;
FrontDistMin{1}=0.75*((xFront{1}(2)-xFront{1}(1))^2+(yFront{1}(2)-yFront{1}(1))^2)^0.5;

% -- Mapped Front Coordinates -- %
xmapFront{1}=zeros(1,NFront{1}+2); ymapFront{1}=zeros(1,NFront{1}+2);
umapFront{1}=zeros(1,NFront{1}+2); vmapFront{1}=zeros(1,NFront{1}+2);

for il=2:NFront{1}+1
    for i=1:Nx
        if xFront{1}(il)>xi(i+1), else
            xmapFront{1}(il)=(i-1)/Nx*(xi(i+1)-xFront{1}(il))/(xi(i+1)-xi(i))+...
                (i)/Nx*(xFront{1}(il)-xi(i))/(xi(i+1)-xi(i)); break
        end
    end
    for j=1:Ny
        if yFront{1}(il)>yj(j+1), else
            ymapFront{1}(il)=(j-1)/Ny*(yj(j+1)-yFront{1}(il))/(yj(j+1)-yj(j))+...
                (j)/Ny*(yFront{1}(il)-yj(j))/(yj(j+1)-yj(j)); break
        end
    end
end
xmapFront{1}(1)=-xmapFront{1}(3);  
xmapFront{1}(NFront{1}+2)=xmapFront{1}(NFront{1}+1) + ...
    ( xmapFront{1}(NFront{1}+1) - xmapFront{1}(NFront{1}) ) ;
ymapFront{1}(1)=ymapFront{1}(3);  ymapFront{1}(NFront{1}+2)=ymapFront{1}(NFront{1}+1);
% -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

% Bubble Front -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
NFront{2}=floor(450*Nx/80); xFront{2}=zeros(1,NFront{2}+2); yFront{2}=zeros(1,NFront{2}+2);
uFront{2}=zeros(1,NFront{2}+2); vFront{2}=zeros(1,NFront{2}+2);
rodxFront{2}=zeros(1,NFront{2}+2); rodyFront{2}=zeros(1,NFront{2}+2);
for il=2:NFront{2}+1
    xFront{2}(il)=xdrop-rdrop*sin(-pi + pi*(il-2)/(NFront{2}-1));
    yFront{2}(il)=ydrop+rdrop*cos(-pi + pi*(il-2)/(NFront{2}-1));
end
xFront{2}(1) = -xFront{2}(3); yFront{2}(1) = yFront{2}(3);
xFront{2}(NFront{2}+2) = -xFront{2}(NFront{2}); yFront{2}(NFront{2}+2) = yFront{2}(NFront{2});
% For de Reestructure of the Front 
FrontDistMax{2}=1.95*((xFront{2}(2)-xFront{2}(1))^2+(yFront{2}(2)-yFront{2}(1))^2)^0.5;
FrontDistMin{2}=0.5*((xFront{2}(2)-xFront{2}(1))^2+(yFront{2}(2)-yFront{2}(1))^2)^0.5;

% -- Mapped Front Coordinates -- %
xmapFront{2}=zeros(1,NFront{2}+2); ymapFront{2}=zeros(1,NFront{2}+2);
umapFront{2}=zeros(1,NFront{2}+2); vmapFront{2}=zeros(1,NFront{2}+2);

for il=2:NFront{2}+1
    for i=1:Nx
        if xFront{2}(il)>xi(i+1), else
            xmapFront{2}(il)=(i-1)/Nx*(xi(i+1)-xFront{2}(il))/(xi(i+1)-xi(i))+...
                (i)/Nx*(xFront{2}(il)-xi(i))/(xi(i+1)-xi(i)); break
        end
    end
    for j=1:Ny
        if yFront{2}(il)>yj(j+1), else
            ymapFront{2}(il)=(j-1)/Ny*(yj(j+1)-yFront{2}(il))/(yj(j+1)-yj(j))+...
                (j)/Ny*(yFront{2}(il)-yj(j))/(yj(j+1)-yj(j)); break
        end
    end
end
xmapFront{2}(1) = -xmapFront{2}(3); ymapFront{2}(1) = ymapFront{2}(3);
xmapFront{2}(NFront{2}+2) = -xmapFront{2}(NFront{2}); 
ymapFront{2}(NFront{2}+2) = ymapFront{2}(NFront{2});

xFront{2}(2)=0; xFront{2}(NFront{2}+1)=0;
xmapFront{2}(2)=0; xmapFront{2}(NFront{2}+1)=0;

% -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
% -------------------------------------------------------------------

% Construction of The Marker Function  ------------------------------
MarkF0=MarkF;

[MarkF,MarkFdx,MarkFdy]...
    =ConstructMarkerFunction03_A(MarkF0,MF01,MF02,Nx,Ny,...
    xu,yu,xv,yv,dx_PE,dy_PN,dx_WP,dy_SP,...
    NFront,xFront,yFront,xmapFront,ymapFront,Lx);

for i=2:Nx+1
    for j=2:Ny+1
        MarkF(i,j)=max(MarkF(i,j),0);   MarkF(i,j)=min(MarkF(i,j),1);
    end
end
MarkF(1,:)=MarkF(2,:); MarkF(Nx+2,:)=MarkF(Nx+1,:);

% Viscosity & Density
ro=RO2.*MarkF+(1-MarkF).*RO1;       miu=MIU2.*MarkF+(1-MarkF).*MIU1;

% Surface Tension Force
[sftx,sfty] = SurfaceTensionForce04_D(Nx,Ny,sigma,...
    xu,yu,xv,yv,NFront,xFront,yFront,xmapFront,ymapFront,Lx,Ly);
% -------------------------------------------------------------------
% =======================================================================
