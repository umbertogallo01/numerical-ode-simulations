%Metodo di Runge-Kutta 4 per il secondo esame pratico

clc
clear
close all

N=12000;
T=60;
t=linspace(0,T,N);
h=T/(N-1);

D=3;
u=zeros(D,N); %approssimazione
u(1,1)=0.1;
u(2,1)=0.1;
u(3,1)=0.1;


for n=1:N-1
    u(:,n+1)=RK4(u(:,n),t(n),h,@Harmonic2D);
end


%visualizzazione grafica
p=plot3(u(1,:),u(2,:),u(3,:),'k');
xlabel('y1(t)');
ylabel('y2(t)');
zlabel('y3(t)');
axis([-20 20 -20 40 0 60]);
axis square


%animazione
waitforbuttonpress
pa=plot3(u(1,1),u(2,1),u(3,1),'k');

xlabel('y1(t)');
ylabel('y2(t)');
zlabel('y3(t)');
axis([-20 20 -20 40 0 60]);
title(['Tempo= ',num2str(0)]);
axis square

az=0;

for n=2:N
    view(az,30);
    az=az+0.2;
    pa(1).XData=[pa.XData u(1,n)];
    pa(1).YData=[pa.YData u(2,n)];
    pa(1).ZData=[pa.ZData u(3,n)];
    title(['Tempo= ',num2str(n)]);
    drawnow;
end





%% Funzioni

%Metodo RK4
function unew=RK4(un,tn,h,f)
    unew=un;
    K1=f(tn,un);
    K2=f(tn+h/2,un+(h/2)*K1);
    K3=f(tn+h/2,un+(h/2)*K2);
    K4=f(tn+h,un+h*K3);
    unew=un+h*((1/6)*K1 + (1/3)*K2 + (1/3)*K3 + (1/6)*K4);
end


%Dinamica equazioni di Navier-Stokes
function fx=Harmonic2D(t,u)
    sigma=10;
    pho=28;
    beta=8/3;
    fx=u;
    fx(1)=sigma*(u(2)-u(1));
    fx(2)=pho*u(1)-u(2)-u(1)*u(3);
    fx(3)=u(1)*u(2)-beta*u(3);
end