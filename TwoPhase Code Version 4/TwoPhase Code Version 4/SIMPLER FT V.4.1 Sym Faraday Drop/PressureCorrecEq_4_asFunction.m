%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pressure Correction Eq.
% 2D Finite Volumes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pc,niteraPC,ErrorPC] = PressureCorrecEq_4_asFunction(Nx,Ny,...
    dx_Pe,dy_Pn,dx_wP,dy_sP,dy_un,dy_su,dx_ve,dx_wv,...
    ro,b,ut,vt,aPu,aPv,pc,dt,ro0,nmaxP,TolP)

% The Pressure Correction Equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ErrorPC=100; niteraPC=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coeficientes aW,aE,aN,aS,..., y Termino Independiente
aEpc=zeros(Nx+2,Ny+2); aWpc=zeros(Nx+2,Ny+2); aNpc=zeros(Nx+2,Ny+2);
aSpc=zeros(Nx+2,Ny+2); aPpc=zeros(Nx+2,Ny+2); bpc=zeros(Nx+2,Ny+2);
for i=2:Nx+1
    for j=2:Ny+1
        aEpc(i,j)=0.5*(ro(i,j)*b(i,j)+ro(i+1,j)*b(i+1,j))*(dy_Pn(i,j)+dy_sP(i,j))*...
            0.5*(b(i,j)+b(i+1,j))*(dy_un(i,j)+dy_su(i,j))/aPu(i,j);
        aWpc(i,j)=0.5*(ro(i-1,j)*b(i-1,j)+ro(i,j)*b(i,j))*(dy_Pn(i,j)+dy_sP(i,j))*...
            0.5*(b(i-1,j)+b(i,j))*(dy_un(i-1,j)+dy_su(i-1,j))/aPu(i-1,j);
        aNpc(i,j)=0.5*(ro(i,j)*b(i,j)+ro(i,j+1)*b(i,j+1))*(dx_Pe(i,j)+dx_wP(i,j))*...
            0.5*(b(i,j)+b(i,j+1))*(dx_ve(i,j)+dx_wv(i,j))/aPv(i,j);
        aSpc(i,j)=0.5*(ro(i,j-1)*b(i,j-1)+ro(i,j)*b(i,j))*(dx_Pe(i,j)+dx_wP(i,j))*...
            0.5*(b(i,j-1)+b(i,j))*(dx_ve(i,j-1)+dx_wv(i,j-1))/aPv(i,j-1);
        aPpc(i,j)=aEpc(i,j)+aWpc(i,j)+aNpc(i,j)+aSpc(i,j);
        bpc(i,j)=-0.5*(ro(i,j)*b(i,j)+ro(i+1,j)*b(i+1,j))*(dy_Pn(i,j)+dy_sP(i,j))*...
            ut(i,j)+...
            0.5*(ro(i-1,j)*b(i-1,j)+ro(i,j)*b(i,j))*(dy_Pn(i,j)+dy_sP(i,j))*...
            ut(i-1,j)+...
            -0.5*(ro(i,j)*b(i,j)+ro(i,j+1)*b(i,j+1))*(dx_Pe(i,j)+dx_wP(i,j))*...
            vt(i,j)+...
            0.5*(ro(i,j-1)*b(i,j-1)+ro(i,j)*b(i,j))*(dx_Pe(i,j)+dx_wP(i,j))*...
            vt(i,j-1)+...
            b(i,j)*(dy_Pn(i,j)+dy_sP(i,j))*(dx_Pe(i,j)+dx_wP(i,j))*...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while (ErrorPC>=TolP)&&(niteraPC<=nmaxP) % -> Esto es si se usa metodo iterativo!
    niteraPC=niteraPC+1; pclast=pc;
   % TDMA por Lineas %%%%%%%%%%%%%%%
   % x linea horizontal
    for i=2:Nx+1
        dy=zeros(Ny,1); My=zeros(Ny,Ny);
        j=2; jM=j-1;
        My(jM,jM)=aPpc(i,j); My(jM,jM+1)=-aNpc(i,j);
        dy(jM)=bpc(i,j)+...
            aEpc(i,j)*pc(i+1,j)+aWpc(i,j)*pc(i-1,j)+aNpc(i,j)*pc(i,j-1);

        for j=3:Ny
            jM=j-1;
            dy(jM)=bpc(i,j)+...
                aEpc(i,j)*pc(i+1,j)+aWpc(i,j)*pc(i-1,j); 
            My(jM,jM)=aPpc(i,j); My(jM-1,jM)=-aSpc(i,j); My(jM+1,jM)=-aNpc(i,j);
        end
        j=Ny+1; jM=j-1;
        dy(jM)=bpc(i,j)+...
            aEpc(i,j)*pc(i+1,j)+aWpc(i,j)*pc(i-1,j)+aSpc(i,j)*pc(i,j-1); 
        My(jM,jM)=aPpc(i,j); My(jM,jM-1)=-aSpc(i,j);
        
        % Solucion
        pc(i,2:Ny+1)=TDMA(My,dy);
    end
    % x Linea Vertical
    for j=2:Ny+1
        dx=zeros(Nx,1); Mx=zeros(Nx,Nx);
        
        i=2; iM=i-1;
        Mx(iM,iM)=aPpc(i,j); Mx(iM,iM+1)=-aEpc(i,j);
        dx(iM)=bpc(i,j)+...
            aNpc(i,j)*pc(i,j+1)+aSpc(i,j)*pc(i,j-1)+aWpc(i,j)*pc(i-1,j);

        for i=3:Nx
            iM=i-1;
            dx(iM)=bpc(i,j)+...
                aNpc(i,j)*pc(i,j+1)+aSpc(i,j)*pc(i,j-1); 
            Mx(iM,iM)=aPpc(i,j); Mx(iM-1,iM)=-aWpc(i,j); Mx(iM+1,iM)=-aEpc(i,j);
        end
        i=Nx+1; iM=i-1;
        dx(iM)=bpc(i,j)+...
            aNpc(i,j)*pc(i,j+1)+aSpc(i,j)*pc(i,j-1)+aEpc(i,j)*pc(i+1,j); 
        Mx(iM,iM)=aPpc(i,j); Mx(iM,iM-1)=-aWpc(i,j);
        
        % Solucion
        pc(2:Nx+1,j)=TDMA(Mx,dx);
    end
    
%     ErrorPC=max(max(abs((pc-pclast)./...
%         max(max(abs(pc/100))))))*100;

%     ErrorPC=max(max(abs((pc-pclast)./pc)))*100;
    ErrorPC1=max(max(abs((pc-pclast)./pc)))*100;
    ErrorPC2=max(max(abs((pc-pclast)./mean(mean(abs(pc))))))*100;
%     ErrorPC=max(ErrorPC1,mean([ErrorPC2,ErrorPC1])); 
    ErrorPC=min(ErrorPC1,ErrorPC2*2.5);

end

disp('Sol. Ec. Presion de Correccion')
info=strcat('#Iteraciones realizadas = ',num2str(niteraPC),'_ Error = ',num2str(ErrorPC));
disp(info)


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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%