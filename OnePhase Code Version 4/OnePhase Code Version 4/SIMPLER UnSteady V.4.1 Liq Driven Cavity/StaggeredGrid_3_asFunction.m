% =======================================================================
% A Staggered Grid.
% Control Volume for u & v -> Forward U,V Volume Control
% =======================================================================
function [xP,yP,dx_PE,dx_Pe,dy_PN,dy_Pn,dx_WP,dx_wP,dy_SP,dy_sP,...
    xu,yu,dx_uE,dx_ue,dy_uN,dy_un,dx_Wu,dx_wu,dy_Su,dy_su,...
    xv,yv,dx_vE,dx_ve,dy_vN,dy_vn,dx_Wv,dx_wv,dy_Sv,dy_sv,...
    FMx,FMy,xi,yj,psi,psj]...
    =StaggeredGrid_3_asFunction(X1,X2,Y1,Y2,Nx,Ny)


% =======================================================================
% Transformacion escojida para mallar -----------------------------------
% Tener cuidado! considerar que existe un pto antes y uno despues
% La funcion evaluada en i=1 =0; y i=Nx+1=1

FMx=@(psi) psi; % Nodos espaciados linealmente! 
FMy=@(psi) psi; % Nodos espaciados linealmente!
 
% FMx=@(psi) (exp(psi)-1)/(exp(1)-1); % Nodos espaciados Crecientemente 1! 
% FMy=@(psi) (exp(psi)-1)/(exp(1)-1); % Nodos espaciados Crecientemente 1!
 
% FMx=@(psi) (exp(psi^2)-1)/(exp(1)-1); % Nodos espaciados Crecientemente 2! 
% FMy=@(psi) (exp(psi^2)-1)/(exp(1)-1); % Nodos espaciados Crecientemente 2!

% FMx=@(psi) 0.5*(1+tan(-pi/4+(psi)*pi/2)); % Nodos espaciados haca la mitad! 
% FMy=@(psi) 0.5*(1+tan(-pi/4+(psi)*pi/2)); % Nodos espaciados Crecientemente 3!

% Polinomio Cubico!
% A=-1.25; FMx=@(psi) psi+A*psi.*(1-psi)*(0.5-psi); % Nodos espaciados hacia la mitad! 
% A=-1.25; FMy=@(psi) psi+A*psi.*(1-psi)*(0.5-psi); % Nodos espaciados hacia la mitad! 

% psj=(j-1)/Ny; psi=(i-1)/Nx;
% =======================================================================

% =======================================================================
% % Initialization of variables ----------------------------------------
xP=zeros(Nx+2,Ny+2); yP=zeros(Nx+2,Ny+2);
dx_PE=zeros(Nx+2,Ny+2); dx_Pe=zeros(Nx+2,Ny+2); dy_PN=zeros(Nx+2,Ny+2);
dy_Pn=zeros(Nx+2,Ny+2); dx_WP=zeros(Nx+2,Ny+2); dx_wP=zeros(Nx+2,Ny+2); 
dy_SP=zeros(Nx+2,Ny+2); dy_sP=zeros(Nx+2,Ny+2);

xu=zeros(Nx+1,Ny+2); % yu=zeros(Nx+1,Ny+2);
dx_uE=zeros(Nx+1,Ny+2); dx_ue=zeros(Nx+1,Ny+2); dy_uN=zeros(Nx+1,Ny+2); 
dy_un=zeros(Nx+1,Ny+2); dx_Wu=zeros(Nx+1,Ny+2); dx_wu=zeros(Nx+1,Ny+2); 
dy_Su=zeros(Nx+1,Ny+2); dy_su=zeros(Nx+1,Ny+2);

yv=zeros(Nx+2,Ny+1); % xv=zeros(Nx+2,Ny+1);
dx_vE=zeros(Nx+2,Ny+1); dx_ve=zeros(Nx+2,Ny+1); dy_vN=zeros(Nx+2,Ny+1); 
dy_vn=zeros(Nx+2,Ny+1); dx_Wv=zeros(Nx+2,Ny+1); dx_wv=zeros(Nx+2,Ny+1); 
dy_Sv=zeros(Nx+2,Ny+1); dy_sv=zeros(Nx+2,Ny+1);

psi=zeros(Nx+2,1); psj=zeros(1,Ny+2); xi=zeros(Nx+1,1); yj=zeros(1,Ny+1);
% -----------------------------------------------------------------------

% Control Volume for Displaced Locations variables ----------------------
for i=1:Nx+1
    psi(i)=(i-1)/(Nx);
    xi(i)=X1+(X2-X1)*FMx(psi(i));
    for j=1:Ny+2
        xu(i,j)=X1+(X2-X1)*FMx(psi(i));
    end
end
for j=1:Ny+1
    psj(j)=(j-1)/(Ny);
    yj(j)=Y1+(Y2-Y1)*FMy(psj(j));
    for  i=1:Nx+2
        yv(i,j)=Y1+(Y2-Y1)*FMy(psj(j));
    end
end
% -----------------------------------------------------------------------

% Control Volume for Non-Displaced Locations variables ------------------
for j=1:Ny+2
    xP(1,j)=xu(1,j)-0.5*(xu(2,j)-xu(1,j));
    % xv(1,j)=xu(1,j)-0.5*(xu(2,j)-xu(1,j));
    for i=2:Nx+1
        xP(i,j)=0.5*(xu(i,j)+xu(i-1,j));
        % xv(i,j)=0.5*(xu(i,j)+xu(i-1,j));
    end
    xP(Nx+2,j)=xu(Nx+1,j)+0.5*(xu(Nx+1,j)-xu(Nx,j));
    % xv(Nx+2,j)=xu(Nx+1,j)+0.5*(xu(Nx+1,j)-xu(Nx,j));
end
xv=xP(1:Nx+2,1:Ny+1);

for i=1:Nx+2
    yP(i,1)=yv(i,1)-0.5*(yv(i,2)-yv(i,1));
    % yu(i,1)=yv(i,1)-0.5*(yv(i,2)-yv(i,1));
    for j=2:Ny+1
        yP(i,j)=0.5*(yv(i,j)+yv(i,j-1));
        % yu(i,j)=0.5*(yv(i,j)+yv(i,j-1));
    end
    yP(i,Ny+2)=yv(i,Ny+1)+0.5*(yv(i,Ny+1)-yv(i,Ny));
    % yu(i,Ny+2)=yv(i,Ny+1)+0.5*(yv(i,Ny+1)-yv(i,Ny));
end
yu=yP(1:Nx+1,1:Ny+2);

% =======================================================================
% Calculos de los delta de espaciamiento ---------------------------------
% Se supone que las caras del VC estan a la mitad entre los nodos
   
for i=2:Nx+1 % Volumenes Internos ----------------------------------------
    for j=2:Ny+1
        dx_PE(i,j)=xP(i+1,j)-xP(i,j); dx_WP(i,j)=xP(i,j)-xP(i-1,j);
        dx_Pe(i,j)=(xP(i+1,j)-xP(i,j))/2; dx_wP(i,j)=(xP(i,j)-xP(i-1,j))/2;
        dy_PN(i,j)=yP(i,j+1)-yP(i,j); dy_SP(i,j)=yP(i,j)-yP(i,j-1);
        dy_Pn(i,j)=(yP(i,j+1)-yP(i,j))/2; dy_sP(i,j)=(yP(i,j)-yP(i,j-1))/2;
    end
end
% ------------------------------------------------------------------------
for i=2:Nx % Volumenes Internos - Velocidad u ----------------------------
    for j=2:Ny+1
        dx_uE(i,j)=xu(i+1,j)-xu(i,j); dx_Wu(i,j)=xu(i,j)-xu(i-1,j);
        dx_ue(i,j)=0.5*(xu(i+1,j)-xu(i,j)); dx_wu(i,j)=0.5*(xu(i,j)-xu(i-1,j));
        dy_uN(i,j)=yu(i,j+1)-yu(i,j); dy_Su(i,j)=yu(i,j)-yu(i,j-1);
        dy_un(i,j)=0.5*(yu(i,j+1)-yu(i,j)); dy_su(i,j)=0.5*(yu(i,j)-yu(i,j-1));
    end
end
% ------------------------------------------------------------------------
for i=2:Nx+1 % Volumenes Internos - Velocidad v --------------------------
    for j=2:Ny
        dx_vE(i,j)=xv(i+1,j)-xv(i,j); dx_Wv(i,j)=xv(i,j)-xv(i-1,j);
        dx_ve(i,j)=0.5*(xv(i+1,j)-xv(i,j)); dx_wv(i,j)=0.5*(xv(i,j)-xv(i-1,j));
        dy_vN(i,j)=yv(i,j+1)-yv(i,j); dy_Sv(i,j)=yv(i,j)-yv(i,j-1);
        dy_vn(i,j)=0.5*(yv(i,j+1)-yv(i,j)); dy_sv(i,j)=0.5*(yv(i,j)-yv(i,j-1));
    end
end
% ------------------------------------------------------------------------
for j=2:Ny+1 % Volumenes Laterales  --------------------------------------
    i=1;
    dx_PE(i,j)=xP(i+1,j)-xP(i,j); dx_WP(i,j)=0;
    dx_Pe(i,j)=(xP(i+1,j)-xP(i,j))/2; dx_wP(i,j)=0;
    dy_PN(i,j)=yP(i,j+1)-yP(i,j); dy_SP(i,j)=yP(i,j)-yP(i,j-1);
    dy_Pn(i,j)=(yP(i,j+1)-yP(i,j))/2; dy_sP(i,j)=(yP(i,j)-yP(i,j-1))/2;
    i=Nx+2;
    dx_PE(i,j)=0; dx_WP(i,j)=xP(i,j)-xP(i-1,j);
    dx_Pe(i,j)=0; dx_wP(i,j)=(xP(i,j)-xP(i-1,j))/2;
    dy_PN(i,j)=yP(i,j+1)-yP(i,j); dy_SP(i,j)=yP(i,j)-yP(i,j-1);
    dy_Pn(i,j)=(yP(i,j+1)-yP(i,j))/2; dy_sP(i,j)=(yP(i,j)-yP(i,j-1))/2;
end
% ------------------------------------------------------------------------
for j=2:Ny+1 % Volumenes Laterales - Velocidad u -------------------------
    i=1;
    dx_uE(i,j)=xu(i+1,j)-xu(i,j); dx_Wu(i,j)=0;
    dx_ue(i,j)=(xu(i+1,j)-xu(i,j))/2; dx_wu(i,j)=0;
    dy_uN(i,j)=yu(i,j+1)-yu(i,j); dy_Su(i,j)=yu(i,j)-yu(i,j-1);
    dy_un(i,j)=(yu(i,j+1)-yu(i,j))/2; dy_su(i,j)=(yu(i,j)-yu(i,j-1))/2;
    i=Nx+1;
    dx_uE(i,j)=0; dx_Wu(i,j)=xu(i,j)-xu(i-1,j);
    dx_ue(i,j)=0; dx_wu(i,j)=(xu(i,j)-xu(i-1,j))/2;
    dy_uN(i,j)=yu(i,j+1)-yu(i,j); dy_Su(i,j)=yu(i,j)-yu(i,j-1);
    dy_un(i,j)=(yu(i,j+1)-yu(i,j))/2; dy_su(i,j)=(yu(i,j)-yu(i,j-1))/2;
end
% ------------------------------------------------------------------------
for j=2:Ny % Volumenes Laterales - Velocidad v ---------------------------
    i=1;
    dx_vE(i,j)=xv(i+1,j)-xv(i,j); dx_Wv(i,j)=0;
    dx_ve(i,j)=(xv(i+1,j)-xv(i,j))/2; dx_wv(i,j)=0;
    dy_vN(i,j)=yv(i,j+1)-yv(i,j); dy_Sv(i,j)=yv(i,j)-yv(i,j-1);
    dy_vn(i,j)=(yv(i,j+1)-yv(i,j))/2; dy_sv(i,j)=(yv(i,j)-yv(i,j-1))/2;
    i=Nx+2;
    dx_vE(i,j)=0; dx_Wv(i,j)=xv(i,j)-xv(i-1,j);
    dx_ve(i,j)=0; dx_wv(i,j)=(xv(i,j)-xv(i-1,j))/2;
    dy_vN(i,j)=yv(i,j+1)-yv(i,j); dy_Sv(i,j)=yv(i,j)-yv(i,j-1);
    dy_vn(i,j)=(yv(i,j+1)-yv(i,j))/2; dy_sv(i,j)=(yv(i,j)-yv(i,j-1))/2;
end
% ------------------------------------------------------------------------
for i=2:Nx+1 % Volumenes Superior e inferior -----------------------------
    j=1;
    dx_PE(i,j)=xP(i+1,j)-xP(i,j); dx_WP(i,j)=xP(i,j)-xP(i-1,j);
    dx_Pe(i,j)=(xP(i+1,j)-xP(i,j))/2; dx_wP(i,j)=(xP(i,j)-xP(i-1,j))/2;
    dy_PN(i,j)=yP(i,j+1)-yP(i,j); dy_SP(i,j)=0;
    dy_Pn(i,j)=(yP(i,j+1)-yP(i,j))/2; dy_sP(i,j)=0;
    j=Ny+2;
    dx_PE(i,j)=xP(i+1,j)-xP(i,j); dx_WP(i,j)=xP(i,j)-xP(i-1,j);
    dx_Pe(i,j)=(xP(i+1,j)-xP(i,j))/2; dx_wP(i,j)=(xP(i,j)-xP(i-1,j))/2;
    dy_SP(i,j)=yP(i,j)-yP(i,j-1); dy_PN(i,j)=0; 
    dy_sP(i,j)=(yP(i,j)-yP(i,j-1))/2; dy_Pn(i,j)=0;
end
% ------------------------------------------------------------------------
for i=2:Nx % Volumenes Superior e inferior - Velocidad u -----------------
    j=1;
    dx_uE(i,j)=xu(i+1,j)-xu(i,j); dx_Wu(i,j)=xu(i,j)-xu(i-1,j);
    dx_ue(i,j)=(xu(i+1,j)-xu(i,j))/2; dx_wu(i,j)=(xu(i,j)-xu(i-1,j))/2;
    dy_uN(i,j)=yu(i,j+1)-yu(i,j); dy_Su(i,j)=0;
    dy_un(i,j)=(yu(i,j+1)-yu(i,j))/2; dy_su(i,j)=0;
    j=Ny+2;
    dx_uE(i,j)=xu(i+1,j)-xu(i,j); dx_Wu(i,j)=xu(i,j)-xu(i-1,j);
    dx_ue(i,j)=(xu(i+1,j)-xu(i,j))/2; dx_wu(i,j)=(xu(i,j)-xu(i-1,j))/2;
    dy_Su(i,j)=yu(i,j)-yu(i,j-1); dy_uN(i,j)=0; 
    dy_su(i,j)=(yu(i,j)-yu(i,j-1))/2; dy_un(i,j)=0;
end
% ------------------------------------------------------------------------
for i=2:Nx+1 % Volumenes Superior e inferior - Velocidad v ---------------
    j=1;
    dx_vE(i,j)=xv(i+1,j)-xv(i,j); dx_Wv(i,j)=xv(i,j)-xv(i-1,j);
    dx_ve(i,j)=(xv(i+1,j)-xv(i,j))/2; dx_wv(i,j)=(xv(i,j)-xv(i-1,j))/2;
    dy_vN(i,j)=yv(i,j+1)-yv(i,j); dy_Sv(i,j)=0;
    dy_vn(i,j)=(yv(i,j+1)-yv(i,j))/2; dy_sv(i,j)=0;
    j=Ny+1;
    dx_vE(i,j)=xv(i+1,j)-xv(i,j); dx_Wv(i,j)=xv(i,j)-xv(i-1,j);
    dx_ve(i,j)=(xv(i+1,j)-xv(i,j))/2; dx_wv(i,j)=(xv(i,j)-xv(i-1,j))/2;
    dy_Sv(i,j)=yv(i,j)-yv(i,j-1); dy_vN(i,j)=0; 
    dy_sv(i,j)=(yv(i,j)-yv(i,j-1))/2; dy_vn(i,j)=0;
end
% ------------------------------------------------------------------------
i=1; j=1; % Volumenes de las esquinas ------------------------------------
dx_PE(i,j)=xP(i+1,j)-xP(i,j); dx_WP(i,j)=0;
dx_Pe(i,j)=(xP(i+1,j)-xP(i,j))/2; dx_wP(i,j)=0;
dy_PN(i,j)=yP(i,j+1)-yP(i,j); dy_SP(i,j)=0;
dy_Pn(i,j)=(yP(i,j+1)-yP(i,j))/2; dy_sP(i,j)=0;
i=1; j=Ny+2;
dx_PE(i,j)=xP(i+1,j)-xP(i,j); dx_WP(i,j)=0;
dx_Pe(i,j)=(xP(i+1,j)-xP(i,j))/2; dx_wP(i,j)=0;
dy_SP(i,j)=yP(i,j)-yP(i,j-1); dy_PN(i,j)=0;
dy_sP(i,j)=(yP(i,j)-yP(i,j-1))/2; dy_Pn(i,j)=0;
i=Nx+2; j=1;
dx_PE(i,j)=0; dx_WP(i,j)=xP(i,j)-xP(i-1,j);
dx_Pe(i,j)=0; dx_wP(i,j)=(xP(i,j)-xP(i-1,j))/2;
dy_PN(i,j)=yP(i,j+1)-yP(i,j); dy_SP(i,j)=0;
dy_Pn(i,j)=(yP(i,j+1)-yP(i,j))/2; dy_sP(i,j)=0;
i=Nx+2; j=Ny+2;
dx_PE(i,j)=0; dx_WP(i,j)=xP(i,j)-xP(i-1,j);
dx_Pe(i,j)=0; dx_wP(i,j)=(xP(i,j)-xP(i-1,j))/2;
dy_SP(i,j)=yP(i,j)-yP(i,j-1); dy_PN(i,j)=0;
dy_sP(i,j)=(yP(i,j)-yP(i,j-1))/2; dy_Pn(i,j)=0;
% ------------------------------------------------------------------------

i=1; j=1; % Volumenes de las esquinas - Velocidad u ----------------------
dx_uE(i,j)=xu(i+1,j)-xu(i,j); dx_Wu(i,j)=0;
dx_ue(i,j)=(xu(i+1,j)-xu(i,j))/2; dx_wu(i,j)=0;
dy_uN(i,j)=yu(i,j+1)-yu(i,j); dy_Su(i,j)=0;
dy_un(i,j)=(yu(i,j+1)-yu(i,j))/2; dy_su(i,j)=0;
i=1; j=Ny+2;
dx_uE(i,j)=xu(i+1,j)-xu(i,j); dx_Wu(i,j)=0;
dx_ue(i,j)=(xu(i+1,j)-xu(i,j))/2; dx_wu(i,j)=0;
dy_Su(i,j)=yu(i,j)-yu(i,j-1); dy_uN(i,j)=0;
dy_su(i,j)=(yu(i,j)-yu(i,j-1))/2; dy_un(i,j)=0;
i=Nx+1; j=1;
dx_uE(i,j)=0; dx_Wu(i,j)=xu(i,j)-xu(i-1,j);
dx_ue(i,j)=0; dx_wu(i,j)=(xu(i,j)-xu(i-1,j))/2;
dy_uN(i,j)=yu(i,j+1)-yu(i,j); dy_Su(i,j)=0;
dy_un(i,j)=(yu(i,j+1)-yu(i,j))/2; dy_su(i,j)=0;
i=Nx+1; j=Ny+2;
dx_uE(i,j)=0; dx_Wu(i,j)=xu(i,j)-xu(i-1,j);
dx_ue(i,j)=0; dx_wu(i,j)=(xu(i,j)-xu(i-1,j))/2;
dy_Su(i,j)=yu(i,j)-yu(i,j-1); dy_uN(i,j)=0;
dy_su(i,j)=(yu(i,j)-yu(i,j-1))/2; dy_un(i,j)=0;
% ------------------------------------------------------------------------

i=1; j=1; % Volumenes de las esquinas - Velocidad v ----------------------
dx_vE(i,j)=xv(i+1,j)-xv(i,j); dx_Wv(i,j)=0;
dx_ve(i,j)=(xv(i+1,j)-xv(i,j))/2; dx_wv(i,j)=0;
dy_vN(i,j)=yv(i,j+1)-yv(i,j); dy_Sv(i,j)=0;
dy_vn(i,j)=(yv(i,j+1)-yv(i,j))/2; dy_sv(i,j)=0;
i=1; j=Ny+1;
dx_vE(i,j)=xv(i+1,j)-xv(i,j); dx_Wv(i,j)=0;
dx_ve(i,j)=(xv(i+1,j)-xv(i,j))/2; dx_wv(i,j)=0;
dy_Sv(i,j)=yv(i,j)-yv(i,j-1); dy_vN(i,j)=0;
dy_sv(i,j)=(yv(i,j)-yv(i,j-1))/2; dy_vn(i,j)=0;
i=Nx+2; j=1;
dx_vE(i,j)=0; dx_Wv(i,j)=xv(i,j)-xv(i-1,j);
dx_ve(i,j)=0; dx_wv(i,j)=(xv(i,j)-xv(i-1,j))/2;
dy_vN(i,j)=yv(i,j+1)-yv(i,j); dy_Sv(i,j)=0;
dy_vn(i,j)=(yv(i,j+1)-yv(i,j))/2; dy_sv(i,j)=0;
i=Nx+2; j=Ny+1;
dx_vE(i,j)=0; dx_Wv(i,j)=xv(i,j)-xv(i-1,j);
dx_ve(i,j)=0; dx_wv(i,j)=(xv(i,j)-xv(i-1,j))/2;
dy_Sv(i,j)=yv(i,j)-yv(i,j-1); dy_vN(i,j)=0;
dy_sv(i,j)=(yv(i,j)-yv(i,j-1))/2; dy_vn(i,j)=0;
% ------------------------------------------------------------------------
display('Mallado->Listo')
% =======================================================================
