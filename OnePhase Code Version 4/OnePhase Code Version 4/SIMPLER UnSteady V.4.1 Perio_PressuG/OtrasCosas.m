% =======================================================================
% Otras variables
% =======================================================================

% =======================================================================
% Inicializacion de Otras Variables 
% interespointsu=[126,125,124,123,110,95,80,65,59,37,23,14,10,9,8];
% interespointsv=[125,124,123,122,117,111,104,65,31,30,21,13,11,10,9];
% % interespoints=[124/4;104/4;];
% % u velocity along Vertical Line througth geometry center
% uvels=zeros(length(interespointsu),1); vvels=zeros(length(interespointsv),1);
% yuvels=zeros(length(interespointsu),1); xvvels=zeros(length(interespointsv),1);
% for i=1:length(interespointsu)
%         uvels(i)=0.5*(u(Nx/2+1,interespointsu(i))+u(Nx/2+1,interespointsu(i)+1));
%         vvels(i)=0.5*(v(interespointsv(i),Ny/2+1)+v(interespointsv(i)+1,Ny/2+1));
%         yuvels(i)=0.5*(yu(Nx/2+1,interespointsu(i))+yu(Nx/2+1,interespointsu(i)+1));
%         xvvels(i)=0.5*(xv(interespointsv(i),Ny/2+1)+xv(interespointsv(i)+1,Ny/2+1));
% end
% =======================================================================

% =======================================================================
    % u velocity along Vertical Line througth geometry center
    % v velocity along Horizontal Line througth geometry center
% 
%     for i=1:length(interespoints)
%         uvels(i)=0.5*(u(interespoints(i),Ny/2+1)+u(interespoints(i),Ny/2+2));
%         vvels(i)=0.5*(v(Nx/2+1,interespoints(i))+v(Nx/2+1,interespoints(i)));
%     end
%     
%     maxdiffu=max(abs((uvels(1)-0.84123)/0.84123),max(abs((uvels(:)-uvels0(:))./uvels0(:))));
%     maxdiffv=max(abs((vvels(:)-vvels0(:))./vvels0(:)));
%     maxuv=max(maxdiffu,maxdiffv);
%     if maxuv*100>=0.0001, Nt=Nt+1; Nit=0:Nt;
%         ErrorNit0=ErrorNit; niteraNit0=niteraNit;
%         ErrorNitPE0=ErrorNitPE; niteraNitPE0=niteraNitPE;
%         ErrorNitPC0=ErrorNitPC; niteraNitPC0=ErrorNitPC;
%         uvels10=uvels1; vvels10=vvels1;
%         
%         ErrorNit=zeros(Nt,1); niteraNit=zeros(Nt,1);
%         ErrorNitPE=zeros(Nt,1); niteraNitPE=zeros(Nt,1);
%         ErrorNitPC=zeros(Nt,1); niteraNitPC=zeros(Nt,1);
%         uvels1=zeros(Nt,1); vvels1=zeros(Nt,1);
%         
%         ErrorNit(1:Nt-1)=ErrorNit0; niteraNit(1:Nt-1)=niteraNit0;
%         ErrorNitPE(1:Nt-1)=ErrorNitPE0; niteraNitPE(1:Nt-1)=niteraNitPE0;
%         ErrorNitPC(1:Nt-1)=ErrorNitPC0; niteraNitPC(1:Nt-1)=niteraNitPC0;
%         uvels1(1:Nt-1)=uvels10; vvels1(1:Nt-1)=vvels10;
%     end
%     
%     disp('-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.')
%     info=strcat('maxuv % # = ',num2str(maxuv*100));
%     disp(info)
%     disp('-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.')
%     if it==(Nt-10), dt=dt*50; end

% -----------------------------------------------------------------------
%     for i=1:length(interespoints)
%         uvels(i)=0.5*(u(interespoints(i),Ny/2+1)+u(interespoints(i),Ny/2+2));
%         vvels(i)=0.5*(v(Nx/2+1,interespoints(i))+v(Nx/2+1,interespoints(i)));
%     end
%     for i=1:length(interespointsu)
%         uvels(i)=0.5*(u(Nx/2+1,interespointsu(i))+u(Nx/2+1,interespointsu(i)+1));
%         vvels(i)=0.5*(v(interespointsv(i),Ny/2+1)+v(interespointsv(i)+1,Ny/2+1));
%         yuvels(i)=0.5*(yu(Nx/2+1,interespointsu(i))+yu(Nx/2+1,interespointsu(i)+1));
%         xvvels(i)=0.5*(xv(interespointsv(i),Ny/2+1)+xv(interespointsv(i)+1,Ny/2+1));
%     end
% =======================================================================
