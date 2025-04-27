%2022.10.27
%chenyixin
clear all
clc
data.Nx=100;
data.dt=0.01;data.c=10;
data.left=-50;data.right=50;
data.x=linspace(data.left,data.right,data.Nx);
data.dx=(data.right-data.left)/data.Nx;
data.Nt=6000;
% [q,Uint,Pint]=initial(data);
% plot(data.x,Pint) 展示初始化inc
% [incleft,incright]=incregion(data);
% U=zeros(1,data.Nx);
% P=zeros(1,data.Nx);
% for tstep=0:1:data.Nt
%     for i=incleft:1:incright
%         [~,Ui,Pi]=buildinc(data,tstep,i);
%         U(i)=Ui;P(i)=Pi;
%     end
%     pause(data.dt)
%     plot(data.x,P)%展示不可压缩场变量
%     axis([data.left data.right -60 60])
% end
[qint,Uint,Pint]=initial(data);
q=qint;
figure(2)
plot(data.x,q(3,:))
pause(2)
tmpq=zeros(3,data.Nx);
timeadvance_RK34(data,q);
function timeadvance_RK34(data,q)
for tstep=1:1:data.Nt
    q0=q;
%     for i=1:1:data.Nx
%         [qfi]=filter(q0,i,data);
%         tmpq(:,i)=qfi;
%     end
%     q0=tmpq;
    [Q0]=makeQ(q0,data,tstep);    
    q1=q+0.169193539*data.dt*Q0;
%     for i=1:1:data.Nx
%         [qfi]=filter(q1,i,data);
%         tmpq(:,i)=qfi;
%     end
%     q1=tmpq;
    [Q1]=makeQ(q1,data,tstep);
    q2=q+0.1874412*data.dt*Q1;
%     for i=1:1:data.Nx
%         [qfi]=filter(q2,i,data);
%         tmpq(:,i)=qfi;
%     end
%     q2=tmpq;
    [Q2]=makeQ(q2,data,tstep);
    q3=q+0.25*data.dt*Q2;
%     for i=1:1:data.Nx
%         [qfi]=filter(q3,i,data);
%         tmpq(:,i)=qfi;
%     end
%     q3=tmpq;
    [Q3]=makeQ(q3,data,tstep);
    q4=q+1.0/3*data.dt*Q3;
%     for i=1:1:data.Nx
%         [qfi]=filter(q4,i,data);
%         tmpq(:,i)=qfi;
%     end
%     q4=tmpq;
    [Q4]=makeQ(q4,data,tstep);
    q5=q+0.5*data.dt*Q4;
%     for i=1:1:data.Nx
%         [qfi]=filter(q5,i,data);
%         tmpq(:,i)=qfi;
%     end
%     q5=tmpq;
    [Q5]=makeQ(q5,data,tstep);
    newq=q+data.dt*Q5;
    for i=1:1:data.Nx
        [qfi]=filter(newq,i,data);
        tmpq(:,i)=qfi;
    end
    newq=tmpq;
    q=newq;
    figure(2)
    plot(data.x,q(3,:))
    axis([data.left data.right -1 1])
    title(['tstep=',num2str(tstep)])
    pause(0.02)
end
end
function [Q]=makeQ(q,data,tstep)
boundi=3;
Q=zeros(3,data.Nx);
[dqdxleft,dqdxright]=drpboundary(q,data);

[dqdx]=compact4order(dqdxleft,dqdxright,q,data);

for i=1:1:data.Nx
if (i>boundi&&i<data.Nx-boundi+1)
    [Qi]=makeright(q,i,data,tstep,dqdx);
    Q(:,i)=Qi;
else
    [Qi]=boundmakeright(q,i,data,tstep);
    Q(:,i)=Qi;
end

end



%%comp  bound

end
function [Q]=makeright(q,i,data,tstep,dqdx)

S=buildS(q,i,data,tstep);
[rho1,U,P]=buildinc(data,tstep,i);
u=U+q(2,i);rho=1+rho1+q(1,i);c=(1.4*P/rho)^0.5;
%A=[u rho;data.c^2/rho u];
A=[u rho 0.0;u^2 2*rho*u 1.0;u*c^2 rho*c^2 0];
%%%
dqdxi=dqdx(:,i);
%Q=S-A*dqdx-[0;artdamp;0];
Q=S-A*dqdxi;
end
function [Q]=boundmakeright(q,i,data,tstep)
[~,U,~]=buildinc(data,tstep,i);
V=U+q(2,i);
%V=V+data.c;
%[dqdxleft,dqdxright]=compboundary(q,data);
[dqdxleft,dqdxright]=drpboundary(q,data);
if i==1 
    dqdx=dqdxleft(:,1);
elseif i==2
    dqdx=dqdxleft(:,2);
elseif i==3
    dqdx=dqdxleft(:,3);
elseif i==data.Nx
    dqdx=dqdxright(:,1);
elseif i==data.Nx-1
    dqdx=dqdxright(:,2);
elseif i==data.Nx-2
    dqdx=dqdxright(:,3);
end
%[~,artdamp]=drp7damp(q,i,data);
%Q=-V*(dqdx+q/(2*data.x(i)))-[0;artdamp];
Q=-V*(dqdx+q(:,i)/(2*abs(data.x(i))));
% S=buildS(q,i,data,tstep);
% [rho1,U,P]=buildinc(data,tstep,i);
% u=U+q(2,i);rho=1+rho1+q(1,i);
%A=[u rho;data.c^2/rho u];
%A=[U 1.0 0.0;0.0 U 1.0;0.0 1.414*P U];
% A=[u rho 0.0;u^2 2*rho*u 1.0;u*data.c^2 rho*data.c^2 0];
% %%%
% Q=S-A*dqdx;
end
function [dqdxleft,dqdxright]=drpboundary(q,data)
a06_0=-2.192280339;a06_1=4.748611401;a06_2=-5.108851915;a06_3=4.461567104;a06_4=-2.833498741;a06_5=1.128328861;a06_6=-0.203876371;
dqdx1=(a06_0*q(:,1)+a06_1*q(:,2)+a06_2*q(:,3)+a06_3*q(:,4)+a06_4*q(:,5)+a06_5*q(:,6)+a06_6*q(:,7))/data.dx;
a60_0=2.192280339;a60_1=-4.748611401;a60_2=5.108851915;a60_3=-4.461567104;a60_4=2.833498741;a60_5=-1.128328861;a60_6=0.203876371;
dqdxNx=(a60_0*q(:,data.Nx)+a60_1*q(:,data.Nx-1)+a60_2*q(:,data.Nx-2)+a60_3*q(:,data.Nx-3)+a60_4*q(:,data.Nx-4)+a60_5*q(:,data.Nx-5)+a60_6*q(:,data.Nx-6))/data.dx;
a15f1=-0.209337622;a15_0=-1.084875676;a15_1=2.147776050;a15_2=-1.388928322;a15_3=0.768949766;a15_4=-0.281814650;a15_5=0.048230454;
dqdx2=(a15f1*q(:,1)+a15_0*q(:,2)+a15_1*q(:,3)+a15_2*q(:,4)+a15_3*q(:,5)+a15_4*q(:,6)+a15_5*q(:,7))/data.dx;
a51_1=0.209337622;a51_0=1.084875676;a51f1=-2.147776050;a51f2=1.388928322;a51f3=-0.768949766;a51f4=0.281814650;a51f5=-0.048230454;
dqdxNx1=(a51_1*q(:,data.Nx)+a51_0*q(:,data.Nx-1)+a51f1*q(:,data.Nx-2)+a51f2*q(:,data.Nx-3)+a51f3*q(:,data.Nx-4)+a51f4*q(:,data.Nx-5)+a51f5*q(:,data.Nx-6))/data.dx;
a24f2=0.049041958;a24f1=-0.468840357;a24_0=-0.474760914;a24_1=1.273274737;a24_2=-0.518484526;a24_3=0.166138533;a24_4=-0.026369431;
dqdx3=(a24f2*q(:,1)+a24f1*q(:,2)+a24_0*q(:,3)+a24_1*q(:,4)+a24_2*q(:,5)+a24_3*q(:,6)+a24_4*q(:,7))/data.dx;
a42_2=-0.049041958;a42_1=0.468840357;a42_0=0.474760914;a42f1=-1.273274737;a42f2=0.518484526;a42f3=-0.166138533;a42f4=0.026369431;
dqdxNx2=(a42_2*q(:,data.Nx)+a42_1*q(:,data.Nx-1)+a42_0*q(:,data.Nx-2)+a42f1*q(:,data.Nx-3)+a42f2*q(:,data.Nx-4)+a42f3*q(:,data.Nx-5)+a42f4*q(:,data.Nx-6))/data.dx;


a0=0;a1=0.77088238051822552;a2=-0.166705904414580469;a3=0.02084314277031176;
dqdx4=(a0*q(:,4)+a1*q(:,5)+a2*q(:,6)+a3*q(:,7)-a1*q(:,3)-a2*q(:,2)-a3*q(:,1))/data.dx;
dqdxNx3=(a0*q(:,data.Nx-3)+a1*q(:,data.Nx-2)+a2*q(:,data.Nx-1)+a3*q(:,data.Nx)-a1*q(:,data.Nx-4)-a2*q(:,data.Nx-5)-a3*q(:,data.Nx-6))/data.dx;

dqdxleft(:,1)=dqdx1;dqdxleft(:,2)=dqdx2;dqdxleft(:,3)=dqdx3;dqdxleft(:,4)=dqdx4;
dqdxright(:,1)=dqdxNx;dqdxright(:,2)=dqdxNx1;dqdxright(:,3)=dqdxNx2;dqdxright(:,4)=dqdxNx3;
end

function [dqdx]=compact4order(dqdxleft,dqdxright,q,data)
lambda=0.3821038098462933;
a=0.7940346032820977;
b=0.0440346032820977;
    dqdx(:,1)=dqdxleft(:,1);
    dqdx(:,2)=dqdxleft(:,2);
    dqdx(:,3)=dqdxleft(:,3);
    dqdx(:,4)=dqdxleft(:,4);
    dqdx(:,data.Nx)=dqdxright(:,1);

    dqdx(:,data.Nx-1)=dqdxright(:,2);

    dqdx(:,data.Nx-2)=dqdxright(:,3);
    dqdx(:,data.Nx-3)=dqdxright(:,4);
A_x=zeros(data.Nx-8,data.Nx-8);
d_x=zeros(3,data.Nx-8);    
for j=1:1:data.Nx-8
    if j==1
        A_x(1,1)=1;
        A_x(1,2)=lambda;
    elseif j==data.Nx-8
        A_x(data.Nx-8,data.Nx-8)=1;
        A_x(data.Nx-8,data.Nx-9)=lambda;    
    else
        A_x(j,j)=1;
        A_x(j,j+1)=lambda;
        A_x(j,j-1)=lambda;
    end
    
end

for j=1:1:data.Nx-8
          if j==1
                d_x(:,1)=a*(q(:,6)-q(:,4))/data.dx+b*(q(:,7)-q(:,3))/data.dx-lambda*dqdx(:,4);
          elseif j==data.Nx-8
                d_x(:,data.Nx-8)=a*(q(:,data.Nx-3)-q(:,data.Nx-5))/data.dx+b*(q(:,data.Nx-2)-q(:,data.Nx-6))/data.dx-lambda*dqdx(:,data.Nx-3);          
          else
                d_x(:,j)=a*(q(:,j+5)-q(:,j+3))/data.dx+b*(q(:,j+6)-q(:,j+2))/data.dx;
          end
end   
     for k=1:1:3
     b=d_x(k,:)';
     x = tridiagonalSolver(A_x, b);
     dqdx(k,5:data.Nx-4)=x;
     end
    
    
    

end

function x = tridiagonalSolver(A, b)
    % 获取矩阵 A 的维度
    n = size(A, 1);
    
    % 初始化变量
    alpha = zeros(n, 1);
    beta = zeros(n, 1);
    x = zeros(n, 1);
    
    % 追赶法的前向消元
    alpha(1) = -A(1, 2) / A(1, 1);
    beta(1) = b(1) / A(1, 1);
    for i = 2:n-1
        alpha(i) = -A(i, i+1) / (A(i, i) + A(i, i-1) * alpha(i-1));
        beta(i) = (b(i) - A(i, i-1) * beta(i-1)) / (A(i, i) + A(i, i-1) * alpha(i-1));
    end
    
    % 追赶法的回代求解
    x(n) = (b(n) - A(n, n-1) * beta(n-1)) / (A(n, n) + A(n, n-1) * alpha(n-1));
    for i = n-1:-1:1
        x(i) = alpha(i) * x(i+1) + beta(i);
    end
end




function [dqdx]=drp7damp(q,i,data)
a0=0;a1=0.77088238051822552;a2=-0.166705904414580469;a3=0.02084314277031176;
miua=0.5*(data.dx^2);
d0=0.287392842460;d1=-0.226146951809;d2=0.106303578770;d3=-0.023853048191;
dqdx=(a0*q(:,i)+a1*q(:,i+1)+a2*q(:,i+2)+a3*q(:,i+3)-a1*q(:,i-1)-a2*q(:,i-2)-a3*q(:,i-3))/data.dx;
%artdamp=miua*(d0*q(2,i)+d1*q(2,i+1)+d2*q(2,i+2)+d3*q(2,i+3)+d1*q(2,i-1)+d2*q(2,i-2)+d3*q(2,i-3))/(data.dx^2);
end
function S=buildS(q,i,data,tstep)
% w=0.1;
% incleft=floor(data.Nx/2)-4;
% incright=floor(data.Nx/2)+4-1;
% t=tstep*data.dt;
% [rho1,U,~]=buildinc(data,tstep,i);
% rho0=1;rho=rho1+rho0+q(1,i);
% if i<incright&&i>incleft
% dPdt=0.0;
% dPdx=0;
% else
%     dPdt=0.0;
% end
% u=U+q(2,i);
% % coe_rho=(rho1+q(1,i))/rho;
% % drho1dt=w^2/(data.c^2)*data.x(i)*sin(w*t);
% % drho1dx=-w*cos(w*t)/data.c^2;
% % dUdt=w*cos(w*t);
% % dUdx=0;
% %S1=drho1dt+u*drho1dx;
% %S2=coe_rho*dUdt+(coe_rho*U+q(2,i))*dUdx;
% S3=-dPdt;
S=[0;0;0];
end
function [qfi]=filter(q,i,data)
xigema=0.5;
qghostNx=0.5*(q(:,data.Nx)+q(:,data.Nx-1));
qghost1=0.5*(q(:,1)+q(:,2));
if i>5&&i<data.Nx-4
    Df=xigema*(63/256*q(:,i)-105/513*q(:,i+1)-105/513*q(:,i-1)+15/128*q(:,i+2)+15/128*q(:,i-2)-45/1024*q(:,i+3)-45/1024*q(:,i-3)+5/512*q(:,i+4)+5/512*q(:,i-4)-1/1024*q(:,i+5)-1/1024*q(:,i-5));
    qfi=q(:,i)-Df;
elseif i==5||i==data.Nx-4
    Df=xigema*(35/128*q(:,i)-7/32*q(:,i+1)-7/32*q(:,i-1)+7/64*q(:,i+2)+7/64*q(:,i-2)-1.0/32*q(:,i+3)-1.0/32*q(:,i-3)+1/256*q(:,i+4)+1/256*q(:,i-4));
    qfi=q(:,i)-Df;
elseif i==4||i==data.Nx-3
    Df=xigema*(5/16*q(:,i)-15/64*q(:,i+1)-15/64*q(:,i-1)+3/32*q(:,i+2)+3/32*q(:,i-2)-1.0/64*q(:,i+3)-1.0/64*q(:,i-3));
    qfi=q(:,i)-Df;
elseif i==3||i==data.Nx-2
    Df=xigema*(3/8*q(:,i)-1/4*q(:,i+1)-1/4*q(:,i-1)+1/6*q(:,i+2)+1/6*q(:,i-2));
    qfi=q(:,i)-Df;
elseif i==2||i==data.Nx-1
    Df=xigema*(1/2*q(:,i)-1/4*q(:,i+1)-1/4*q(:,i-1));
    qfi=q(:,i)-Df;
elseif i==data.Nx
    Df=xigema*(1/2*q(:,i)-1/4*qghostNx-1/4*q(:,i-1));
    qfi=q(:,i)-Df;
elseif i==1
    Df=xigema*(1/2*q(:,i)-1/4*qghost1-1/4*q(:,i+1));
    qfi=q(:,i)-Df;
end
end
function [incleft,incright]=incregion(data)
incleft=floor(data.Nx/2)-4;
incright=floor(data.Nx/2)+4-1;
end
function [q,Uint,Pint]=initial(data) %初始化U,P,rho1,q=(rho' u') %
q=zeros(3,data.Nx);
Uint=zeros(1,data.Nx);
Pint=zeros(1,data.Nx);
[incleft,incright]=incregion(data);
for i=1:1:data.Nx
    [rho1,Uinti,Pinti]=buildinc(data,0,i);
    Uint(i)=Uinti;Pint(i)=Pinti;% 应该给赋值-rho1吗？存在间断
%     q(1,i)=exp(-0.693147*data.x(i)^2/9)/data.c^2+exp(-0.693147*(data.x(i)-50)^2/9)/data.c^2;
%     q(3,i)=exp(-0.693147*(data.x(i)-50)^2/9)+exp(-0.693147*data.x(i)^2/9);
     q(1,i)=exp(-0.693147*data.x(i)^2/9)/data.c^2;
     q(3,i)=exp(-0.693147*data.x(i)^2/9);
end

end
function [rho1,U,P]=buildinc(data,tstep,i) %输出一个点的某一时刻的U P rho1%
[incleft,incright]=incregion(data);
P0=data.c^2*1.0/1.4;
w=0.1;t=tstep*data.dt;
if i<incright&&i>incleft
    P=P0;
    U=0.0;
    rho1=0;
else
    P=P0;
    U=0;
    rho1=0;
end


%这是一个满足一维不可压缩NS方程的U和P解




end