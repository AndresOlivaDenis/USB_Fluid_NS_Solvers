% =======================================================================
% Plot Mesh Results
% =======================================================================

if (PlotM==1)
    figure; 
    subplot(3,3,[3]);
    plot([1:Nx+2]',xP(:,1),'-o',[1:Nx+1]',xu(:,1),'o',[1:Nx+2]',xv(:,1),'o')
    legend('Escalar','u','v'); 
    title('Espaciamiento en x');  xlabel('Nx'); ylabel('x')

    subplot(3,3,[6]);
    plot([1:Ny+2]',yP(1,:),'-o',[1:Ny+2]',yu(1,:),'o',[1:Ny+1]',yv(1,:),'o')
    legend('Escalar','u','v'); 
    title('Espaciamiento en y');  xlabel('Ny'); ylabel('y')

    subplot(3,3,[1 2 4 5 7 8]);
    xx(1:Nx+1,1:Ny+1)=0.5*(xu(1:Nx+1,2:Ny+2)+xu(1:Nx+1,1:Ny+1));
    yy(1:Nx+1,1:Ny+1)=0.5*(yv(2:Nx+2,1:Ny+1)+yv(1:Nx+1,1:Ny+1));
    mesh(xx,yy,zeros(Nx+1, Ny+1))
    view(2)
    
    if ((Plot1 == 1)||(Plot1 == 2))
        figure;
    end
end
