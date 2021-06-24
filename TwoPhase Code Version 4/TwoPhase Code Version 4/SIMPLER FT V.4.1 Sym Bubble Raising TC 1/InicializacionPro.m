% ======================================================================= %
% Inicializacion de propiedades y variables
% =======================================================================

% Initial Marker Function -----------------------------------------------
MF01=1; MF02=0; MarkF=MF02*ones(Nx+2,Ny+2);
for i=2:Nx+1,
    for j=2:Ny+1
        if (((xP(i,j)-xdrop)^2+(yP(i,j)-ydrop)^2)<=rdrop^2)
            MarkF(i,j)=MF01; end
    end,
end
MarkF(1,:)=MarkF(2,:);
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
NFront=floor(450*h/80)*2; xFront=zeros(1,NFront+2); yFront=zeros(1,NFront+2);
uFront=zeros(1,NFront+2); vFront=zeros(1,NFront+2);
rodxFront=zeros(1,NFront+2); rodyFront=zeros(1,NFront+2);
for il=2:NFront+1
    xFront(il)=xdrop-rdrop*sin(-pi + pi*(il-2)/(NFront-1));
    yFront(il)=ydrop+rdrop*cos(-pi + pi*(il-2)/(NFront-1));
end
xFront(1) = -xFront(3); yFront(1) = yFront(3);
xFront(NFront+2) = -xFront(NFront); yFront(NFront+2) = yFront(NFront);
% For de Reestructure of the Front 
FrontDistMax=1.95*((xFront(2)-xFront(1))^2+(yFront(2)-yFront(1))^2)^0.5;
FrontDistMin=0.5*((xFront(2)-xFront(1))^2+(yFront(2)-yFront(1))^2)^0.5;

% -- Mapped Front Coordinates -- %
xmapFront=zeros(1,NFront+2); ymapFront=zeros(1,NFront+2);
umapFront=zeros(1,NFront+2); vmapFront=zeros(1,NFront+2);

for il=2:NFront+1
    for i=1:Nx
        if xFront(il)>xi(i+1), else
            xmapFront(il)=(i-1)/Nx*(xi(i+1)-xFront(il))/(xi(i+1)-xi(i))+...
                (i)/Nx*(xFront(il)-xi(i))/(xi(i+1)-xi(i)); break
        end
    end
    for j=1:Ny
        if yFront(il)>yj(j+1), else
            ymapFront(il)=(j-1)/Ny*(yj(j+1)-yFront(il))/(yj(j+1)-yj(j))+...
                (j)/Ny*(yFront(il)-yj(j))/(yj(j+1)-yj(j)); break
        end
    end
end
xmapFront(1) = -xmapFront(3); ymapFront(1) = ymapFront(3);
xmapFront(NFront+2) = -xmapFront(NFront); ymapFront(NFront+2) = ymapFront(NFront);

% -------------------------------------------------------------------
xFront(2)=0; xFront(NFront+1)=0;
xmapFront(2)=0; xmapFront(NFront+1)=0;

% Construction of The Marker Function  ------------------------------
MarkF0=MarkF;
[MarkF,MarkFdx,MarkFdy]...
    =ConstructMarkerFunction03_A(MarkF0,MF01,MF02,Nx,Ny,...
    xu,yu,xv,yv,dx_PE,dy_PN,dx_WP,dy_SP,...
    NFront,xFront,yFront,xmapFront,ymapFront);
for i=2:Nx+1
    for j=2:Ny+1
        MarkF(i,j)=max(MarkF(i,j),0); MarkF(i,j)=min(MarkF(i,j),1);
    end
end
MarkF(1,:)=MarkF(2,:); MarkF(Nx+2,:)=MarkF(Nx+1,:);
% Viscosity & Density
ro=RO2.*MarkF+(1-MarkF).*RO1; miu=MIU2.*MarkF+(1-MarkF).*MIU1;
% Surface Tension Force
[sftx,sfty] = SurfaceTensionForce04_D(Nx,Ny,sigma,...
    xu,yu,xv,yv,NFront,xFront,yFront,xmapFront,ymapFront);
% -------------------------------------------------------------------
% =======================================================================
