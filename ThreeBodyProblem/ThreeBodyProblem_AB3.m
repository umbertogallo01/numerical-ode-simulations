%Terzo Appello MNA (Gallo Umberto 1937779)

clc
clear
close all

T=200;
N=20001;
t=linspace(0,T,N);
h=T/(N-1);
p=3; %numero di passi di AB3
dim=12;
b0=23/12;
b1=-16/12;
b2=5/12;


%Definizione Approssimazione
U=zeros(dim,N);
U(:,1)=[-1 -1.6 -1 1.5 2 0 0.5 0.2 0.4 0.6 0.6 0.3];


%Calcolo di U(:,2) e U(:,3) con un metodo di ordine 2 perché AB3 ha ordine
%3. In particolare usiamo Heun
U(:,2)=U(:,1)+(0.5*h)*(F(U(:,1))+F(U(:,1)+h*F(U(:,1))));
U(:,3)=U(:,2)+(0.5*h)*(F(U(:,2))+F(U(:,2)+h*F(U(:,2))));


%Metodo AB3
for n=p:N-1
    U(:,n+1)=AB3(U,n,h,@F,b0,b1,b2);
end

load MNA_2023_APP_III.mat
norm(U-Ucheck,'inf')


%Visualizzazione dei risultati
p=plot(U(1,:),U(2,:),'r-',U(3,:),U(4,:),'b-',U(5,:),U(6,:),'g-');
axis square
legend('P1','P2','P3');

%Animazione
waitforbuttonpress
pa=plot(U(1,1),U(2,1),'r-',U(3,1),U(4,1),'b-',U(5,1),U(6,1),'g-');
title('Tempo= ',num2str(0));
axis square

%Handle agli assi
a=gca;

for n=1:10:N-1
    pa(1).XData=[pa(1).XData U(1,n+1)];
    pa(1).YData=[pa(1).YData U(2,n+1)];

    pa(2).XData=[pa(2).XData U(3,n+1)];
    pa(2).YData=[pa(2).YData U(4,n+1)];

    pa(3).XData=[pa(3).XData U(5,n+1)];
    pa(3).YData=[pa(3).YData U(6,n+1)];
    
    bx=(U(1,n+1)+U(3,n+1)+U(5,n+1))/3;
    by=(U(2,n+1)+U(4,n+1)+U(6,n+1))/3;

    a.XLim=[bx-5, bx+5];
    a.YLim=[by-5, by+5];

    drawnow;
    title('Tempo= ',num2str(n));
end



%Funzioni

%AB3
function unew=AB3(U,n,h,F,b0,b1,b2)
    unew=U(:,n)+h*(b0*F(U(:,n))+b1*F(U(:,n-1))+b2*F(U(:,n-2)));
end


%Dinamica
function fx=F(u)

    m1=1; m2=2; m3=3; epsilon=2; 
    fx=u;
    r12=[u(1)-u(3),u(2)-u(4)];
    r13=[u(1)-u(5),u(2)-u(6)];
    r23=[u(3)-u(5),u(4)-u(6)];

    R12=(norm(r12)^2+epsilon)^(3/2);
    R13=(norm(r13)^2+epsilon)^(3/2);
    R23=(norm(r23)^2+epsilon)^(3/2);

    fx(1)=u(7);
    fx(2)=u(8);
    fx(3)=u(9);
    fx(4)=u(10);
    fx(5)=u(11);
    fx(6)=u(12);
    fx(7)=-( (m2*(u(1)-u(3)))/R12 + (m3*(u(1)-u(5)))/R13 );
    fx(8)=-( (m2*(u(2)-u(4)))/R12 + (m3*(u(2)-u(6)))/R13 );
    fx(9)=-( (m1*(u(3)-u(1)))/R12 + (m3*(u(3)-u(5)))/R23 );
    fx(10)=-( (m1*(u(4)-u(2)))/R12 + (m3*(u(4)-u(6)))/R23 );
    fx(11)=-( (m1*(u(5)-u(1)))/R13 + (m2*(u(5)-u(3)))/R23 );
    fx(12)=-( (m1*(u(6)-u(2)))/R13 + (m2*(u(6)-u(4)))/R23 );
end

