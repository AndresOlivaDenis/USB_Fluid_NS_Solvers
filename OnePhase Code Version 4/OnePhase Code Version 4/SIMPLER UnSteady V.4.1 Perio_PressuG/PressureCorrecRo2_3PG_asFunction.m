% ========================================================================
% Pressure Correction Eq.
% ========================================================================

function [pc,niteraPC,ErrorPC] = PressureCorrecRo2_3PG_asFunction(Nx,Ny,...
    dx_Pe,dy_Pn,dx_wP,dy_sP,dy_un,dy_su,dx_ve,dx_wv,...
    ro,ut,vt,aPu,aPv,pc,dt,ro0,nmaxP,TolP,...
    dpdy,yP)

% The Pressure Correction Equation
% ========================================================================
ErrorPC=100; niteraPC=0; omepc=1.6;


% ========================================================================
% Coeficientes aW,aE,aN,aS,..., y Termino Independiente
aEpc=zeros(Nx+2,Ny+2); aWpc=zeros(Nx+2,Ny+2); aNpc=zeros(Nx+2,Ny+2);
aSpc=zeros(Nx+2,Ny+2); aPpc=zeros(Nx+2,Ny+2); bpc=zeros(Nx+2,Ny+2);
for i=2:Nx+1
    for j=2:Ny+1
        aEpc(i,j)=(dy_Pn(i,j)+dy_sP(i,j))*(dy_un(i,j)+dy_su(i,j))/aPu(i,j);
        aWpc(i,j)=(dy_Pn(i,j)+dy_sP(i,j))*(dy_un(i-1,j)+dy_su(i-1,j))/aPu(i-1,j);
        aNpc(i,j)=(dx_Pe(i,j)+dx_wP(i,j))*(dx_ve(i,j)+dx_wv(i,j))/aPv(i,j);
        aSpc(i,j)=(dx_Pe(i,j)+dx_wP(i,j))*(dx_ve(i,j-1)+dx_wv(i,j-1))/aPv(i,j-1);
        aPpc(i,j)=aEpc(i,j)+aWpc(i,j)+aNpc(i,j)+aSpc(i,j);
        bpc(i,j)=-(dy_Pn(i,j)+dy_sP(i,j))*ut(i,j)+...
            (dy_Pn(i,j)+dy_sP(i,j))*ut(i-1,j)+...
            -(dx_Pe(i,j)+dx_wP(i,j))*vt(i,j)+...
            (dx_Pe(i,j)+dx_wP(i,j))*vt(i,j-1)+...
            (dy_Pn(i,j)+dy_sP(i,j))*(dx_Pe(i,j)+dx_wP(i,j))*...
            (-ro(i,j)+ro0(i,j))/dt;
    end
end

i=2; j=2; aPpc(i,j)=aEpc(i,j)+aNpc(i,j); % Esquina 1

i=2; % Nodos del tope izquierdo
for j=3:Ny,aPpc(i,j)=aEpc(i,j)+aNpc(i,j)+aSpc(i,j);end

i=2; j=Ny+1; aPpc(i,j)=aEpc(i,j)+aSpc(i,j); % Esquina  3

for i=3:Nx
    j=2; aPpc(i,j)=aEpc(i,j)+aWpc(i,j)+aNpc(i,j); % Nodos del tope inferior
    j=Ny+1; aPpc(i,j)=aEpc(i,j)+aWpc(i,j)+aSpc(i,j); % Nodos del tope superior
end

i=Nx+1; j=2; aPpc(i,j)=aWpc(i,j)+aNpc(i,j); % Esquina 2

i=Nx+1; % Nodos del tope derecho
for j=3:Ny, aPpc(i,j)=aWpc(i,j)+aNpc(i,j)+aSpc(i,j); end

i=Nx+1; j=Ny+1; aPpc(i,j)=aWpc(i,j)+aSpc(i,j); % Esquina  4

aWpc(2,:)=zeros(1,Ny+2); aEpc(Nx+1,:)=zeros(1,Ny+2);
aNpc(:,Ny+1)=zeros(Nx+2,1); aSpc(:,2)=zeros(Nx+2,1);

% ========================================================================

while (ErrorPC>=TolP)&&(niteraPC<=nmaxP) % -> Esto es si se usa metodo iterativo!
    niteraPC=niteraPC+1; pclast=pc;
   % Gradiente de Presion en y ----------------------------
    pc(2:Nx+1,2)= pc(2:Nx+1,Ny+1) - ...
        dpdy*(yP(2:Nx+1,Ny+1)-yP(2:Nx+1,2));
    pc(2:Nx+1,Ny+1)= pc(2:Nx+1,3) - ...
        dpdy*(yP(2:Nx+1,3)-yP(2:Nx+1,Ny+1));
    % -----------------------------------------------------
    for i=2:Nx+1
        for j=3:Ny
            pc(i,j)=omepc*1/aPpc(i,j)*...
                (aEpc(i,j)*pc(i+1,j)+aWpc(i,j)*pc(i-1,j)+...
                aNpc(i,j)*pc(i,j+1)+aSpc(i,j)*pc(i,j-1)+...
                bpc(i,j))+...
                (1-omepc)*pc(i,j);
        end
    end
%     ErrorPC=max(max(abs((pc-pclast)./...
%         max(max(abs(pc/100))))))*100;

%     ErrorPC=max(max(abs((pc-pclast)./pc)))*100;
%     ErrorPC1=max(max(abs((pc-pclast)./pc)))*100;
%     ErrorPC2=max(max(abs((pc-pclast)./mean(mean(abs(pc))))))*100;
% %     ErrorPC=max(ErrorPC1,mean([ErrorPC2,ErrorPC1])); 
%     ErrorPC=min(ErrorPC1,ErrorPC2);

    ErrorPC=(sum((pc-pclast).^2)/sum(pc.^2))^0.5;
    ErrorPC=ErrorPC*100;
end

disp('Sol. Ec. Presion de Correccion')
info=strcat('#Iteraciones realizadas = ',num2str(niteraPC),'_ Error = ',num2str(ErrorPC));
disp(info)

% ========================================================================

% i=5; j=5; disp('Coeficientes Nodos Internos')
% aEpc(i,j)/aPpc(i,j)
% 0.5*(ro(i,j)*b(i,j)+ro(i+1,j)*b(i+1,j))*(dy_Pn(i,j)+dy_sP(i,j))/aPpc(i,j)
% 
% disp('Coeficientes Nodos izq'); i=2; j=5;
% aEpc(i,j)/(aEpc(i,j)+aNpc(i,j)+aSpc(i,j))
% 0.5*(ro(i,j)*b(i,j)+ro(i+1,j)*b(i+1,j))*(dy_Pn(i,j)+dy_sP(i,j))/(aEpc(i,j)+aNpc(i,j)+aSpc(i,j))
% 
% disp('Coeficientes Nodos der'); i=Nx+1; j=5;
% aWpc(i,j)/(aWpc(i,j)+aNpc(i,j)+aSpc(i,j))
% 0.5*(ro(i,j)*b(i,j)+ro(i+1,j)*b(i+1,j))*(dy_Pn(i,j)+dy_sP(i,j))/(aWpc(i,j)+aNpc(i,j)+aSpc(i,j))
%         
% disp('Coeficientes Nodos sup'); i=5; j=Ny+1;
% aEpc(i,j)/(aEpc(i,j)+aWpc(i,j)+aSpc(i,j))
% 0.5*(ro(i,j)*b(i,j)+ro(i+1,j)*b(i+1,j))*(dy_Pn(i,j)+dy_sP(i,j))/(aEpc(i,j)+aWpc(i,j)+aSpc(i,j))
%         
% disp('Coeficientes Nodos infe'); i=5; j=2;
% aEpc(i,j)/(aEpc(i,j)+aWpc(i,j)+aNpc(i,j))
% 0.5*(ro(i,j)*b(i,j)+ro(i+1,j)*b(i+1,j))*(dy_Pn(i,j)+dy_sP(i,j))/(aEpc(i,j)+aWpc(i,j)+aNpc(i,j))

% disp('Coeficientes Nodos esquina1'); i=2; j=2;
% aEpc(i,j)/(aEpc(i,j)+aNpc(i,j))
% 0.5*(ro(i,j)*b(i,j)+ro(i+1,j)*b(i+1,j))*(dy_Pn(i,j)+dy_sP(i,j))/(aEpc(i,j)+aNpc(i,j))
% 
% disp('Coeficientes Nodos esquina2'); i=Nx+1; j=2;
% aWpc(i,j)/(aWpc(i,j)+aNpc(i,j))
% 0.5*(ro(i,j)*b(i,j)+ro(i+1,j)*b(i+1,j))*(dy_Pn(i,j)+dy_sP(i,j))/(aWpc(i,j)+aNpc(i,j))
% 
% disp('Coeficientes Nodos esquina3'); i=2; j=Ny+1;
% aSpc(i,j)/(aEpc(i,j)+aSpc(i,j))
% 0.5*(ro(i,j)*b(i,j)+ro(i+1,j)*b(i+1,j))*(dy_Pn(i,j)+dy_sP(i,j))/(aEpc(i,j)+aSpc(i,j))
%     
% disp('Coeficientes Nodos esquina4'); i=Nx+1; j=Ny+1;
% aWpc(i,j)/(aWpc(i,j)+aSpc(i,j))
% 0.5*(ro(i,j)*b(i,j)+ro(i+1,j)*b(i+1,j))*(dy_Pn(i,j)+dy_sP(i,j))/(aWpc(i,j)+aSpc(i,j))
    
% ========================================================================