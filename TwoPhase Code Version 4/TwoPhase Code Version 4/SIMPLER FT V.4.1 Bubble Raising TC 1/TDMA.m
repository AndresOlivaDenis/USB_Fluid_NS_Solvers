%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algoritmo de Matrices Tridiagonales, Metodo Thomas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M-> Matriz Tridiagonal
% b-> Vector de Terminos Independientes
function Tsol=TDMA(M,d)
N=length(d); P=zeros(N,1); Q=zeros(N,1);
P(1)=-M(1,2)/M(1,1); Q(1)=d(1)/M(1,1);
for i=2:N-1
    P(i)=-M(i,i+1)/(M(i,i)+M(i,i-1)*P(i-1));
    Q(i)=(-M(i,i-1)*Q(i-1)+d(i))/(M(i,i)+M(i,i-1)*P(i-1));
end
i=N; Q(i)=(-M(i,i-1)*Q(i-1)+d(i))/(M(i,i)+M(i,i-1)*P(i-1));
Tsol(N)=Q(N);
for i=1:N-1
    Tsol(N-i)=P(N-i)*Tsol(N-i+1)+Q(N-i);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
