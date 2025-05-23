%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution of 2D steady momentum eq.
% Time -> General CN
% Space -> Centered Diference 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ut,vt,aPu,aPv] = MomentumMiu2_STFg_asFunction(Nx,Ny,...
    dx_uE,dx_ue,dy_uN,dy_un,dx_Wu,dx_wu,dy_Su,dy_su,...
    dx_vE,dx_ve,dy_vN,dy_vn,dx_Wv,dx_wv,dy_Sv,dy_sv,...
    ro,miu,b,Su,Sv,roref,...
    p,u0,v0,ro0,miu0,dt,...
    uCBnorth,uCBsouth,uCBwest,uCBeast,vCBnorth,vCBsouth,vCBwest,vCBeast,...
    CuCB,CvCB,fc,fd,fSu,fSv,...
    sftx,sfty,sftx0,sfty0,Tol)

Error=100; nitera=0; nmax=1000; % omeu=1.5; omev=1.5;
% fc=0.5; fd=0.5; fSu=0.5; fSv=0.5; % Nsim=1; % Time Eval Weight,
%                 0->Full explicit, 1-> Full implicit, 0.5-> ClankNicolson
%                 fc -> For Convective terms, fd -> For difusion terms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Condiciones de borde
% unorth=1; usouth=0; uwest=0; ueast=0; 
% vnorth=0; vsouth=0; vwest=0; veast=0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ut=u0; vt=v0; % Para inicio del Proceso iterativo de solucion
%%% Condiciones de borde %%%
% Normal Velocities
ut(1,:)=uCBwest*(1-CuCB(1))  +(ut(2,:)-dx_uE(1,:)*uCBwest)*CuCB(1);
ut(Nx+1,:)=uCBeast*(1-CuCB(3))  +(ut(Nx,:)+dx_Wu(Nx+1,:)*uCBeast)*CuCB(3);
vt(:,1)=vCBsouth*(1-CvCB(4))  +(vt(:,2)-dy_vN(:,1)*vCBsouth)*CvCB(4);
vt(:,Ny+1)=vCBnorth*(1-CvCB(2))  +(vt(:,Ny)+dy_Sv(:,Ny+1)*vCBnorth)*CvCB(2);
% Tangential Velocities
ut(:,1)=(2*uCBsouth-ut(:,2))*(1-CuCB(4))+...
    (ut(:,2)-dy_uN(:,2)*uCBsouth)*CuCB(4);
ut(:,Ny+2)=(2*uCBnorth-ut(:,Ny+1))*(1-CuCB(2))+...
    (ut(:,Ny+1)+dy_Su(:,Ny+2)*uCBnorth)*CuCB(2);
vt(1,:)=(2*vCBwest-vt(2,:))*(1-CvCB(1))+...
    (vt(2,:)-dx_vE(1,:)*vCBwest)*CvCB(1);
vt(Nx+2,:)=(2*vCBeast-vt(Nx+1,:))*(1-CvCB(3))+...
    (vt(Nx+1,:)+dx_Wv(Nx+2,:)*vCBeast)*CvCB(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% aP0 %%%
aP0u=zeros(Nx+1,Ny+2); aP0v=zeros(Nx+2,Ny+1);

aP0u(2:Nx,2:Ny+1)=0.5*(ro0(2:Nx,2:Ny+1)+ro0(3:Nx+1,2:Ny+1))/dt.*...
    b(2:Nx,2:Ny+1).*(dy_un(2:Nx,2:Ny+1)+dy_su(2:Nx,2:Ny+1)).*...
    (dx_ue(2:Nx,2:Ny+1)+dx_wu(2:Nx,2:Ny+1));

aP0v(2:Nx+1,2:Ny)=0.5*(ro0(2:Nx+1,2:Ny)+ro0(2:Nx+1,3:Ny+1))/dt.*...
    b(2:Nx+1,2:Ny).*(dy_vn(2:Nx+1,2:Ny)+dy_sv(2:Nx+1,2:Ny)).*...
    (dx_ve(2:Nx+1,2:Ny)+dx_wv(2:Nx+1,2:Ny));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluacion de los terminos de la iteracion anterior y del termino
% independiente
%%%% Coeficientes para u %%%%
% aEuCero=zeros(Nx+1,Ny+2); aWuCero=zeros(Nx+1,Ny+2); aNuCero=zeros(Nx+1,Ny+2);
% aSuCero=zeros(Nx+1,Ny+2);
bu=zeros(Nx+1,Ny+2); % aPu=zeros(Nx+1,Ny+2); 

for i=2:Nx
    for j=2:Ny+1
        Fe=0.5*(b(i+1,j)*(dy_un(i,j)+dy_su(i,j)))*...
            ro0(i+1,j)*(u0(i,j)+u0(i+1,j));
        De=2*miu0(i+1,j)*(b(i+1,j)*(dy_un(i,j)+dy_su(i,j)))/dx_uE(i,j);
        %aEuCero(i,j)=-0.5*Fe+De;
        
        Fw=0.5*(b(i,j)*(dy_un(i,j)+dy_su(i,j)))*...
            ro0(i,j)*(u0(i-1,j)+u0(i,j));
        Dw=2*miu0(i,j)*(b(i,j)*(dy_un(i,j)+dy_su(i,j)))/dx_Wu(i,j);
        %aWuCero(i,j)=0.5*Fw+Dw;
        
        Fn=0.5*(0.25*(b(i,j)+b(i+1,j)+b(i,j+1)+b(i+1,j+1))*...
            (dx_ue(i,j)+dx_wu(i,j)))*...
            (0.25*(ro0(i,j)+ro0(i+1,j)+ro0(i,j+1)+ro0(i+1,j+1)))*...
            (v0(i,j)+v0(i+1,j));
        Dn=0.25*(miu0(i,j)*b(i,j)+miu0(i+1,j)*b(i+1,j)+...
            miu0(i,j+1)*b(i,j+1)+miu0(i+1,j+1)*b(i+1,j+1))*...
            (dx_wu(i,j)+dx_ue(i,j))/dy_uN(i,j);
        %aNuCero(i,j)=-0.5*Fn+Dn;
        
        Fs=0.5*(0.25*(b(i,j-1)+b(i+1,j-1)+b(i,j)+b(i+1,j))*...
            (dx_ue(i,j)+dx_wu(i,j)))*...
            (0.25*(ro0(i,j)+ro0(i+1,j)+ro0(i,j-1)+ro0(i+1,j-1)))*...
            (v0(i,j-1)+v0(i+1,j-1));
        Ds=0.25*(miu0(i,j-1)*b(i,j-1)+miu0(i+1,j-1)*b(i+1,j-1)+...
            miu0(i,j)*b(i,j)+miu0(i+1,j)*b(i+1,j))*...
            (dx_wu(i,j)+dx_ue(i,j))/dy_Su(i,j);
        %aSuCero(i,j)=0.5*Fs+Ds;
        
        Dnv=0.25*(miu0(i,j)*b(i,j)+miu0(i+1,j)*b(i+1,j)+...
            miu0(i,j+1)*b(i,j+1)+miu0(i+1,j+1)*b(i+1,j+1))*...
            (dx_wu(i,j)+dx_ue(i,j))/(dx_ue(i,j)+dx_wu(i,j));
        
        Dsv=0.25*(miu0(i,j-1)*b(i,j-1)+miu0(i+1,j-1)*b(i+1,j-1)+...
            miu0(i,j)*b(i,j)+miu0(i+1,j)*b(i+1,j))*...
            (dx_wu(i,j)+dx_ue(i,j))/(dx_ue(i,j)+dx_wu(i,j));
        
        bu(i,j)=-p(i+1,j)*(b(i+1,j)*(dy_un(i,j)+dy_su(i,j)))+... % Pressure Terms
            p(i,j)*(b(i,j)*(dy_un(i,j)+dy_su(i,j)))+... 
            (1-fSu)*Su(i,j)*b(i,j)*(0.5*(ro0(i,j)+ro0(i+1,j))-roref)*... % Source Terms
            (dy_un(i,j)+dy_su(i,j))*(dx_ue(i,j)+dx_wu(i,j))+...
            (fSu)*Su(i,j)*b(i,j)*(0.5*(ro(i,j)+ro(i+1,j))-roref)*...
            (dy_un(i,j)+dy_su(i,j))*(dx_ue(i,j)+dx_wu(i,j))+...
            (1-fc)*0.5*(-Fe*(u0(i+1,j)+u0(i,j))+Fw*(u0(i-1,j)+u0(i,j))+... % Convective Terms
            -Fn*(u0(i,j+1)+u0(i,j))+Fs*(u0(i,j-1)+u0(i,j)))+...
            (1-fd)*(De*(u0(i+1,j)-u0(i,j))+Dw*(u0(i-1,j)-u0(i,j))+... % Difusion Terms
            Dn*(u0(i,j+1)-u0(i,j))+Ds*(u0(i,j-1)-u0(i,j)))+...
            (1-fSu)*sftx0(i,j)*... % Surface Tension Terms
            b(i,j)*(dy_un(i,j)+dy_su(i,j))*(dx_ue(i,j)+dx_wu(i,j))+...
            (fSu)*sftx(i,j)*... 
            b(i,j)*(dy_un(i,j)+dy_su(i,j))*(dx_ue(i,j)+dx_wu(i,j))+...
            (1-fd)*(Dnv*(v0(i+1,j)-v0(i,j))-Dsv*(v0(i+1,j-1)-v0(i,j-1))); % Difusion 2 
    end
end
%%%% Coeficientes para v %%%%
% aEvCero=zeros(Nx+1,Ny+2); aWvCero=zeros(Nx+1,Ny+2); aNvCero=zeros(Nx+1,Ny+2);
% aSvCero=zeros(Nx+1,Ny+2);
bv=zeros(Nx+1,Ny+2); % aPu=zeros(Nx+1,Ny+2);

for i=2:Nx+1
    for j=2:Ny
        Fe=0.5*(0.25*(b(i,j)+b(i+1,j)+b(i,j+1)+b(i+1,j+1))*...
            (dy_vn(i,j)+dy_sv(i,j)))*...
            (0.25*(ro0(i,j)+ro0(i+1,j)+ro0(i,j+1)+ro0(i+1,j+1)))*...
            (u0(i,j)+u0(i,j+1));
        De=0.25*(miu0(i,j)*b(i,j)+miu0(i+1,j)*b(i+1,j)+...
            miu0(i,j+1)*b(i,j+1)+miu0(i+1,j+1)*b(i+1,j+1))*...
            (dy_vn(i,j)+dy_sv(i,j))/dx_vE(i,j);
        %aEvCero(i,j)=-0.5*Fe+De;
        
        Fw=0.5*(0.25*(b(i-1,j)+b(i,j)+b(i-1,j+1)+b(i,j+1))*...
            (dy_vn(i,j)+dy_sv(i,j)))*...
            (0.25*(ro0(i,j)+ro0(i-1,j)+ro0(i,j+1)+ro0(i-1,j+1)))*...
            (u0(i-1,j)+u0(i-1,j+1));
        Dw=0.25*(miu0(i-1,j)*b(i-1,j)+miu0(i,j)*b(i,j)+...
            miu0(i-1,j+1)*b(i-1,j+1)+miu0(i,j+1)*b(i,j+1))*...
            (dy_vn(i,j)+dy_sv(i,j))/dx_Wv(i,j);
        %aWvCero(i,j)=0.5*Fw+Dw;
        
        Fn=0.5*(b(i,j+1)*(dx_wv(i,j)+dx_ve(i,j)))*...
            ro0(i,j+1)*(v0(i,j+1)+v0(i,j));
        Dn=2*miu0(i,j+1)*(b(i,j+1)*(dx_wv(i,j)+dx_ve(i,j)))/dy_vN(i,j);
        %aNvCero(i,j)=-0.5*Fn+Dn;
        
        Fs=0.5*(b(i,j)*(dx_wv(i,j)+dx_ve(i,j)))*...
            ro0(i,j)*(v0(i,j-1)+v0(i,j));
        Ds=2*miu0(i,j)*(b(i,j)*(dx_wv(i,j)+dx_ve(i,j)))/dy_Sv(i,j);
        %aSvCero(i,j)=0.5*Fs+Ds;
        
        Deu=0.25*(miu0(i,j)*b(i,j)+miu0(i+1,j)*b(i+1,j)+...
            miu0(i,j+1)*b(i,j+1)+miu0(i+1,j+1)*b(i+1,j+1))*...
            (dy_vn(i,j)+dy_sv(i,j))/(dy_vn(i,j)+dy_sv(i,j));
        
        Dwu=0.25*(miu0(i-1,j)*b(i-1,j)+miu0(i,j)*b(i,j)+...
            miu0(i-1,j+1)*b(i-1,j+1)+miu0(i,j+1)*b(i,j+1))*...
            (dy_vn(i,j)+dy_sv(i,j))/(dy_vn(i,j)+dy_sv(i,j));
        
        bv(i,j)=-p(i,j+1)*(b(i,j+1)*(dx_wv(i,j)+dx_ve(i,j)))+... % Pressure 
            p(i,j)*(b(i,j)*(dx_wv(i,j)+dx_ve(i,j)))+...
            (1-fSv)*Sv(i,j)*b(i,j)*(0.5*(ro0(i,j)+ro0(i,j+1))-roref)*... % Source
            (dy_vn(i,j)+dy_sv(i,j))*(dx_ve(i,j)+dx_wv(i,j))+...
            (fSv)*Sv(i,j)*b(i,j)*(0.5*(ro(i,j)+ro(i,j+1))-roref)*...
            (dy_vn(i,j)+dy_sv(i,j))*(dx_ve(i,j)+dx_wv(i,j))+...
            (1-fc)*0.5*(-Fe*(v0(i+1,j)+v0(i,j))+Fw*(v0(i-1,j)+v0(i,j))+... % Convective
            -Fn*(v0(i,j+1)+v0(i,j))+Fs*(v0(i,j-1)+v0(i,j)))+...
            (1-fd)*(De*(v0(i+1,j)-v0(i,j))+Dw*(v0(i-1,j)-v0(i,j))+... % Difusion
            Dn*(v0(i,j+1)-v0(i,j))+Ds*(v0(i,j-1)-v0(i,j)))+...
            (1-fSv)*sfty0(i,j)*... % Surface Tension Terms
            b(i,j)*(dy_vn(i,j)+dy_sv(i,j))*(dx_ve(i,j)+dx_wv(i,j))+...
            (fSv)*sfty(i,j)*... 
            b(i,j)*(dy_vn(i,j)+dy_sv(i,j))*(dx_ve(i,j)+dx_wv(i,j))+...
            (1-fd)*(Deu*(u0(i,j+1)-u0(i,j))-Dwu*(u0(i-1,j+1)-u0(i-1,j))); % Difusion 2
    end
end

while (Error>=Tol)&&(nitera<=nmax) % -> Esto es si se usa metodo iterativo!
    nitera=nitera+1; ulast=ut; vlast=vt;
  
    %info=strcat('#Iteracion = ',num2str(nitera),'_ Error = ',num2str(Error));
    %disp(info)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Coeficientes aW,aE,aN,aS,..., y Termino Independiente 
    % Central Diference Scheme & Linear Interpolation for unstored Values
    
    %%%% Coeficientes para u %%%%
    aEu=zeros(Nx+1,Ny+2); aWu=zeros(Nx+1,Ny+2); aNu=zeros(Nx+1,Ny+2);
    aSu=zeros(Nx+1,Ny+2); aPu=zeros(Nx+1,Ny+2); % bu=zeros(Nx+1,Ny+2);
    Dnv=zeros(Nx+1,Ny+2); Dsv=zeros(Nx+1,Ny+2);
    for i=2:Nx
        for j=2:Ny+1
            Fe=0.5*(b(i+1,j)*(dy_un(i,j)+dy_su(i,j)))*...
                ro(i+1,j)*(ut(i,j)+ut(i+1,j));
            De=2*miu(i+1,j)*(b(i+1,j)*(dy_un(i,j)+dy_su(i,j)))/dx_uE(i,j);
            aEu(i,j)=-fc*0.5*Fe+fd*De;
            
            Fw=0.5*(b(i,j)*(dy_un(i,j)+dy_su(i,j)))*...
                ro(i,j)*(ut(i-1,j)+ut(i,j));
            Dw=2*miu(i,j)*(b(i,j)*(dy_un(i,j)+dy_su(i,j)))/dx_Wu(i,j);
            aWu(i,j)=fc*0.5*Fw+fd*Dw;
            
            Fn=0.5*(0.25*(b(i,j)+b(i+1,j)+b(i,j+1)+b(i+1,j+1))*...
                (dx_ue(i,j)+dx_wu(i,j)))*...
                (0.25*(ro(i,j)+ro(i+1,j)+ro(i,j+1)+ro(i+1,j+1)))*...
                (vt(i,j)+vt(i+1,j));
            Dn=0.25*(miu(i,j)*b(i,j)+miu(i+1,j)*b(i+1,j)+...
                miu(i,j+1)*b(i,j+1)+miu(i+1,j+1)*b(i+1,j+1))*...
                (dx_wu(i,j)+dx_ue(i,j))/dy_uN(i,j);
            aNu(i,j)=-fc*0.5*Fn+fd*Dn;
            
            Fs=0.5*(0.25*(b(i,j-1)+b(i+1,j-1)+b(i,j)+b(i+1,j))*...
                (dx_ue(i,j)+dx_wu(i,j)))*...
                (0.25*(ro(i,j)+ro(i+1,j)+ro(i,j-1)+ro(i+1,j-1)))*...
                (vt(i,j-1)+vt(i+1,j-1));
            Ds=0.25*(miu(i,j-1)*b(i,j-1)+miu(i+1,j-1)*b(i+1,j-1)+...
                miu(i,j)*b(i,j)+miu(i+1,j)*b(i+1,j))*...
                (dx_wu(i,j)+dx_ue(i,j))/dy_Su(i,j);
            aSu(i,j)=fc*0.5*Fs+fd*Ds;
            
            aPu(i,j)=aEu(i,j)+aWu(i,j)+aNu(i,j)+aSu(i,j)+...
                fc*((Fe-Fw)+(Fn-Fs))+...
                0.5*(ro(i,j)+ro(i+1,j))/dt*...
                b(i,j)*(dy_un(i,j)+dy_su(i,j))*(dx_ue(i,j)+dx_wu(i,j));
            
            Dnv(i,j)=0.25*(miu0(i,j)*b(i,j)+miu0(i+1,j)*b(i+1,j)+...
                miu0(i,j+1)*b(i,j+1)+miu0(i+1,j+1)*b(i+1,j+1))*...
                (dx_wu(i,j)+dx_ue(i,j))/(dx_ue(i,j)+dx_wu(i,j));
            
            Dsv(i,j)=0.25*(miu0(i,j-1)*b(i,j-1)+miu0(i+1,j-1)*b(i+1,j-1)+...
                miu0(i,j)*b(i,j)+miu0(i+1,j)*b(i+1,j))*...
                (dx_wu(i,j)+dx_ue(i,j))/(dx_ue(i,j)+dx_wu(i,j));
        end
    end
    
    %%%% Coeficientes para v %%%%
    aEv=zeros(Nx+2,Ny+1); aWv=zeros(Nx+2,Ny+1); aNv=zeros(Nx+2,Ny+1);
    aSv=zeros(Nx+2,Ny+1); aPv=zeros(Nx+2,Ny+1); % bv=zeros(Nx+2,Ny+1);
    Deu=zeros(Nx+2,Ny+1); Dwu=zeros(Nx+2,Ny+1);
    for i=2:Nx+1
        for j=2:Ny
            Fe=0.5*(0.25*(b(i,j)+b(i+1,j)+b(i,j+1)+b(i+1,j+1))*...
                (dy_vn(i,j)+dy_sv(i,j)))*...
                (0.25*(ro(i,j)+ro(i+1,j)+ro(i,j+1)+ro(i+1,j+1)))*...
                (ut(i,j)+ut(i,j+1));
            De=0.25*(miu(i,j)*b(i,j)+miu(i+1,j)*b(i+1,j)+...
                miu(i,j+1)*b(i,j+1)+miu(i+1,j+1)*b(i+1,j+1))*...
                (dy_vn(i,j)+dy_sv(i,j))/dx_vE(i,j);
            aEv(i,j)=-fc*0.5*Fe+fd*De;
            
            Fw=0.5*(0.25*(b(i-1,j)+b(i,j)+b(i-1,j+1)+b(i,j+1))*...
                (dy_vn(i,j)+dy_sv(i,j)))*...
                (0.25*(ro(i,j)+ro(i-1,j)+ro(i,j+1)+ro(i-1,j+1)))*...
                (ut(i-1,j)+ut(i-1,j+1));
            Dw=0.25*(miu(i-1,j)*b(i-1,j)+miu(i,j)*b(i,j)+...
                miu(i-1,j+1)*b(i-1,j+1)+miu(i,j+1)*b(i,j+1))*...
                (dy_vn(i,j)+dy_sv(i,j))/dx_Wv(i,j);
            aWv(i,j)=fc*0.5*Fw+fd*Dw;
            
            Fn=0.5*(b(i,j+1)*(dx_wv(i,j)+dx_ve(i,j)))*...
                ro(i,j+1)*(vt(i,j+1)+vt(i,j));
            Dn=2*miu(i,j+1)*(b(i,j+1)*(dx_wv(i,j)+dx_ve(i,j)))/dy_vN(i,j);
            aNv(i,j)=-fc*0.5*Fn+fd*Dn;
            
            Fs=0.5*(b(i,j)*(dx_wv(i,j)+dx_ve(i,j)))*...
                ro(i,j)*(vt(i,j-1)+vt(i,j));
            Ds=2*miu(i,j)*(b(i,j)*(dx_wv(i,j)+dx_ve(i,j)))/dy_Sv(i,j);
            aSv(i,j)=fc*0.5*Fs+fd*Ds;
            
            aPv(i,j)=aEv(i,j)+aWv(i,j)+aNv(i,j)+aSv(i,j)+...
                fc*((Fe-Fw)+(Fn-Fs))+...
                0.5*(ro(i,j)+ro(i,j+1))/dt*...
                b(i,j)*(dy_vn(i,j)+dy_sv(i,j))*(dx_ve(i,j)+dx_wv(i,j));
            
            Deu(i,j)=0.25*(miu0(i,j)*b(i,j)+miu0(i+1,j)*b(i+1,j)+...
                miu0(i,j+1)*b(i,j+1)+miu0(i+1,j+1)*b(i+1,j+1))*...
                (dy_vn(i,j)+dy_sv(i,j))/(dy_vn(i,j)+dy_sv(i,j));
            
            Dwu(i,j)=0.25*(miu0(i-1,j)*b(i-1,j)+miu0(i,j)*b(i,j)+...
                miu0(i-1,j+1)*b(i-1,j+1)+miu0(i,j+1)*b(i,j+1))*...
                (dy_vn(i,j)+dy_sv(i,j))/(dy_vn(i,j)+dy_sv(i,j));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Solution of discrete momentun eq using the guessed pressure field
    % The term of pressure its already putted in the independent term bu & bv
    
    % Normal Velocities
    ut(1,:)=uCBwest*(1-CuCB(1))  +(ut(2,:)-dx_uE(1,:)*uCBwest)*CuCB(1);
    ut(Nx+1,:)=uCBeast*(1-CuCB(3))  +(ut(Nx,:)+dx_Wu(Nx+1,:)*uCBeast)*CuCB(3);
    vt(:,1)=vCBsouth*(1-CvCB(4))  +(vt(:,2)-dy_vN(:,1)*vCBsouth)*CvCB(4);
    vt(:,Ny+1)=vCBnorth*(1-CvCB(2))  +(vt(:,Ny)+dy_Sv(:,Ny+1)*vCBnorth)*CvCB(2);
    % Tangential Velocities
    ut(:,1)=(2*uCBsouth-ut(:,2))*(1-CuCB(4))+...
        (ut(:,2)-dy_uN(:,2)*uCBsouth)*CuCB(4);
    ut(:,Ny+2)=(2*uCBnorth-ut(:,Ny+1))*(1-CuCB(2))+...
        (ut(:,Ny+1)+dy_Su(:,Ny+2)*uCBnorth)*CuCB(2);
    vt(1,:)=(2*vCBwest-vt(2,:))*(1-CvCB(1))+...
        (vt(2,:)-dx_vE(1,:)*vCBwest)*CvCB(1);
    vt(Nx+2,:)=(2*vCBeast-vt(Nx+1,:))*(1-CvCB(3))+...
        (vt(Nx+1,:)+dx_Wv(Nx+2,:)*vCBeast)*CvCB(3);
    
    for i=2:Nx
        for j=2:Ny+1
            ut(i,j)=1/aPu(i,j)*...
                (aWu(i,j)*ut(i-1,j)+aEu(i,j)*ut(i+1,j)+...
                aNu(i,j)*ut(i,j+1)+aSu(i,j)*ut(i,j-1)+bu(i,j)+...
                aP0u(i,j)*u0(i,j)+...
                fd*(Dnv(i,j)*(vt(i+1,j)-vt(i,j))-...
                Dsv(i,j)*(vt(i+1,j-1)-vt(i,j-1))) );
        end
    end
    for i=2:Nx+1
        for j=2:Ny
            vt(i,j)=1/aPv(i,j)*...
                (aWv(i,j)*vt(i-1,j)+aEv(i,j)*vt(i+1,j)+...
                aNv(i,j)*vt(i,j+1)+aSv(i,j)*vt(i,j-1)+bv(i,j)+...
                aP0v(i,j)*v0(i,j)+...
                fd*(Deu(i,j)*(ut(i,j+1)-ut(i,j))-...
                Dwu(i,j)*(ut(i-1,j+1)-ut(i-1,j))) );
        end
    end
%     Erroru=max(max(abs((ut-ulast)./ulast)));
%     Errorv=max(max(abs((vt-vlast)./vlast)));
%     Error=max(Erroru,Errorv);
    
    Erroru1=max(max(abs((ut-ulast)./ulast)));
    Erroru2=max(max(abs((ut-ulast)./mean(mean(abs(ulast))))));
    Errorv1=max(max(abs((vt-vlast)./vlast)));
    Errorv2=max(max(abs((vt-vlast)./mean(mean(abs(vlast))))));
    Erroru=min(Erroru1,Erroru2*10); Errorv=min(Errorv1,Errorv2*10);
    Error=max(Erroru,Errorv)*100;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Notas:
    %%%  El termino independiente puede sacarse del ciclo iterativo para reducir
    %    costos compu (el desolucion de momentum!)
end

disp('Sol. Ec. Momentun en 2D')
info=strcat('#Iteraciones realizadas = ',num2str(nitera),'_ Error = ',num2str(Error));
disp(info)
    