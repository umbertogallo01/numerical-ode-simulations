%Seconda Prova Metodo di Adams-Moulton

clc
clear
close all

T=60;
N=12000;
t=linspace(0,T,N);
h=T/(N-1);

p=2; %passi del metodo
dim=3; %dimensione del problema

%Definizione dell'approssimazione
U=zeros(dim,N);
U(:,1)=[0.1 0.1 0.1];

%Per trovare l'approssimazione U(:,2) usiamo un metodo che ha ordine di
%conistenza pari a 2 dato che AM2 ha ordine 3 (Ad esempio Heun)
U(:,2)=U(:,1)+(h/2)*(F(t(1),U(:,1))+F(t(2),U(:,1)+h*F(t(1),U(:,1))));

%Definizione dei buffer per le u e per le f
buffer_u=zeros(dim,p);
buffer_f=zeros(dim,p);

for i=1:p
    buffer_u(:,i)=U(:,i);
    buffer_f(:,i)=F(t(i),U(:,i));
end

%Definizione del vettore b del metodo AM2
b=[-1/12 8/12 5/12]; %l'ho invertito


%Definizione dei vettori ausiliari per le iterazioni di punto fisso
aux_k=zeros(dim,1);
aux_kp1=zeros(dim,1);
itmax=70;
tol=1.e-12;

for n=p:N-1
    aux_k=U(:,n);
    it=0;
    err=tol+1;
    while err>tol && it<itmax
        aux_kp1=AM2(buffer_u,buffer_f,t,h,p,b,@F,aux_k,dim,n);
        err=norm(aux_kp1-aux_k,'inf');
        aux_k=aux_kp1;
        it=it+1;
    end

    U(:,n+1)=aux_kp1;

    %Aggiornamento buffer_u
    buffer_u=[buffer_u U(:,n+1)];
    buffer_u(:,1)=[];

    %Aggiornamento buffer_f
    buffer_f=[buffer_f F(t(n+1),U(:,n+1))];
    buffer_f(:,1)=[];
end


%Visualizzazione grafica del risultato
p=plot3(U(1,:),U(2,:),U(3,:));
xlabel('y1(t)');
ylabel('y2(t)');
zlabel('y3(t)');
axis square

%Animazione
waitforbuttonpress
pa=plot3(U(1,1),U(2,1),U(3,1),'k');
title('Tempo= ',num2str(0));
axis([-20 20 -20 40 0 60]);
az=0;

for n=1:N-1
    view(az,30);
    az=az+0.2;
    pa(1).XData=[pa(1).XData U(1,n+1)];
    pa(1).YData=[pa(1).YData U(2,n+1)];
    pa(1).ZData=[pa(1).ZData U(3,n+1)];
    drawnow
    title('Tempo= ',num2str(n));
end



%% Funzioni

%Dinamica
function fx=F(t,u)
    sigma=10;
    pho=28;
    beta=8/3;
    fx=u;
    fx(1)=sigma*(u(2)-u(1));
    fx(2)=pho*u(1)-u(2)-u(1)*u(3);
    fx(3)=u(1)*u(2)-beta*u(3);
end


%Metodo AM2
function unew=AM2(buffer_u,buffer_f,t,h,p,b,F,aux_k,dim,n)
    s=zeros(dim,1);
    for i=1:p
        s=s+b(i)*buffer_f(:,i);
    end
    s=s+b(p+1)*F(t(n+1),aux_k);
    unew=buffer_u(:,p)+h*s;
end


