%Secondo appello Adams-Bashforth (AB2)

clc
clear
close all

T=60;
N=12000;
t=linspace(0,T,N);
h=T/(N-1);
dim=3;
p=2;

%Definizione dell'approssimazione
U=zeros(dim,N);
U(:,1)=[0.1 0.1 0.1];


%Calcolo del valore di innesco con un metodo di ordine 1 (EA) essendo AB2 di
%ordine 2
U(:,2)=U(:,1)+h*F(t(1),U(:,1));



%Definizione dei buffer delle soluzione e delle dinamiche
buffer_u=zeros(dim,p);
buffer_f=zeros(dim,p);

for i=1:p
    buffer_u(:,i)=U(:,i);
    buffer_f(:,i)=F(t(i),U(:,i))
end

%Definizione dei coefficienti del metodo
b=[3/2, -1/2];
b=b(end:-1:1); %inverto il vettore

for n=p:N-1
    U(:,n+1)=AB2(U(:,n),buffer_f,h,p,b,dim);

    %Aggiornamento buffer_u
    buffer_u=[buffer_u U(:,n+1)];
    buffer_u(:,1)=[];

    %Aggiornamento buffer_f
    buffer_f=[buffer_f F(t(n+1),U(:,n+1))];
    buffer_f(:,1)=[];
end

%visualizzazione del risultato
p=plot3(U(1,:),U(2,:),U(3,:),'k');
xlabel('y1(t)')
ylabel('y2(t)')
zlabel('y3(t)')
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

%Metodo AB2
function unew=AB2(Un,buffer_f,h,p,b,dim)
    s=zeros(dim,1);
    for i=1:p
        s(:)=s(:)+b(i)*buffer_f(:,i);
    end
    
    unew=Un+h*s;
end


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



