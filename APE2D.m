%from 2022.11.15 chenyixin
%======================mesh set===============================%
clear all
clc
data.Nx=50;data.Ny=50;%100列，100行
data.dt=1.0;data.c=1.5;data.r0=0.5;
a=400*data.r0;
data.gamma=0.1*pi;
data.top=a;data.bottom=-a;data.left=-a;data.right=a;
data.x=linspace(data.left,data.right,data.Nx);
data.y=linspace(data.bottom,data.top,data.Ny);
data.var=3;
data.xx=zeros(data.Ny,data.Nx);data.yy=zeros(data.Ny,data.Nx);

for k=1:1:data.Nx
data.yy(:,k)=data.y;
end
for k=1:1:data.Ny
data.xx(k,:)=data.x;
end
data.dx=(data.right-data.left)/data.Nx;
data.dy=(data.top-data.bottom)/data.Ny;
data.Nt=20000;CFL=data.c*data.dt/data.dx;
bei=(data.gamma/(0.2*pi))^4/(data.r0^4*data.c^2);
%===========================================================%
% [q]=initial(data);
% figure(2)
% surf(data.xx,data.yy,q(:,:,4))
% shading interp
% [qf]=filter(q,data);
% surf(data.xx,data.yy,qf(:,:,4))
% shading interp
[q]=initial(data);
surf(data.xx,data.yy,q(:,:,1))
shading interp

figure(1)
pause(1)
%a=-0.05;b=0.05;
 a=-3e-5*bei;b=3e-5*bei ;
for tstep=1:1:data.Nt
   %[newq]=LDDRK(q,data,tstep);
     [newq]=RK(q,data,tstep);
    q=newq;
    figure(1)
    surf(data.xx,data.yy,q(:,:,1))
    shading interp
    axis([data.left data.right data.bottom data.top 100*a 100*b])
%     colorbar('YTickLabel',{'Freezing','Cold','Cool',...
%     'Neutral','Warm','Hot','Burning','Nuclear'})
    set(gca,'CLim',[a b])
    title(['tstep=',num2str(tstep)])
    view(0,90)
    xlabel('X')
    ylabel('Y')
    pause(0.0001)
    
end
% function [newq]=LDDRK_LPCE(q,data,tstep)
%     
%     q0=q;
%     [Q0]=makeQ_LPCE(q0,data,tstep);    
%     q1=q+0.169193539*data.dt*Q0;
%     [Q1]=makeQ_LPCE(q1,data,tstep);
%     q2=q+0.1874412*data.dt*Q1;
%     [Q2]=makeQ_LPCE(q2,data,tstep);
%     q3=q+0.25*data.dt*Q2;
%     [Q3]=makeQ_LPCE(q3,data,tstep);
%     q4=q+1.0/3*data.dt*Q3;
%     [Q4]=makeQ_LPCE(q4,data,tstep);
%     q5=q+0.5*data.dt*Q4;
%     [Q5]=makeQ_LPCE(q5,data,tstep);
%     newq=q+data.dt*Q5;
%     [tmpq1]=filter_1(newq,data);
%     %[tmpq2]=filter_x(tmpq1,data);
%     newq=tmpq1;
% end
function [newq]=RK(q,data,tstep)
    
    q0=q;
    [Q0]=makeQ(q0,data,tstep);    
    q1=q+0.5*data.dt*Q0;
    [Q1]=makeQ(q1,data,tstep);
    q2=q+0.5*data.dt*Q1;
    [Q2]=makeQ(q2,data,tstep);
    q3=q+1.0*data.dt*Q2;
    [Q3]=makeQ(q3,data,tstep);
    newq=q+1.0/6*data.dt*(Q0+2*Q1+2*Q2+Q3);
    [tmpq1]=filter(newq,data);
    %[tmpq2]=filter_x(tmpq1,data);
    newq=tmpq1;
end
function [newq]=LDDRK(q,data,tstep)
    
    q0=q;
    [Q0]=makeQ(q0,data,tstep);    
    q1=q+0.169193539*data.dt*Q0;
    [Q1]=makeQ(q1,data,tstep);
    q2=q+0.1874412*data.dt*Q1;
    [Q2]=makeQ(q2,data,tstep);
    q3=q+0.25*data.dt*Q2;
    [Q3]=makeQ(q3,data,tstep);
    q4=q+1.0/3*data.dt*Q3;
    [Q4]=makeQ(q4,data,tstep);
    q5=q+0.5*data.dt*Q4;
    [Q5]=makeQ(q5,data,tstep);
    newq=q+data.dt*Q5;
    [tmpq1]=filter(newq,data);
    %[tmpq2]=filter_x(tmpq1,data);
    newq=tmpq1;
end
function [Q]=makeQ(q,data,tstep)
boundi=5;
Q=zeros(data.Ny,data.Nx,data.var);
[dqdx,dqdy]=drp7boundary(q,data);%
%[dqdx,dqdy]=compact4order(dqdx,dqdy,q,data);
%[dqdx,dqdy]=compact6order(dqdx,dqdy,q,data);
[dqdx,dqdy]=drp7all(q,data,dqdx,dqdy);
%[dqdx,dqdy]=drp7_comp6bound(q,data);
for i=1:1:data.Ny
    for j=1:1:data.Nx
        [Uij,Vij,Pij]=buildinc(data,i,j,tstep);
        rho0=1.0;
        u=Uij+q(i,j,2);v=Vij+q(i,j,3);p=Pij+q(i,j,1);
        ratio=1;%c=(1.4*(abs(Pij+q(i,j,4)))/rho)^0.5;
        if (i>boundi&&i<data.Ny-boundi+1)&&(j>boundi&&j<data.Nx-boundi+1)
            [S]=buildS(q,i,j,data,tstep);
%             A41=ratio*Uij*p/rho^2;
%             B41=ratio*Vij*p/rho^2;
%             A44=ratio*q(i,j,2)+ratio*Uij*q(i,j,1)/rho;
%             B44=ratio*q(i,j,3)+ratio*Vij*q(i,j,1)/rho;
%              A=[u rho 0 0;0 u 0 1.0/rho;0 0 u 0;A41 ratio*p 0 A44;];
%              B=[v 0 rho 0;0 v 0 0;0 0 v 1.0/rho;B41 0 ratio*p B44;];
%             A=[u rho 0 0;0 u 0 1.0/rho;0 0 u 0;0 ratio*p 0 u;];
%             B=[v 0 rho 0;0 v 0 0;0 0 v 1.0/rho;0 0 ratio*p v;];
%            A=[Uij rho0*data.c^2 0;1/rho0 Uij 0;0  0 Uij;];
%            B=[Vij 0 rho0*data.c^2;0 Vij 0;1.0/rho0 0 Vij;];
           A=[0 rho0*data.c^2 0;1/rho0 0 0;0 0 0;];
           B=[0 0 rho0*data.c^2;0 0 0;1.0/rho0 0 0;];
%             A=[Uij rho0 0 0;0 Uij 0 1.0/rho0;0 0 Uij 0;0 ratio*Pij 0 Uij;];
%             B=[Vij 0 rho0 0;0 Vij 0 0;0 0 Vij 1.0/rho0;0 0 ratio*Pij Vij;];
            pqpx=[dqdx(i,j,1);dqdx(i,j,2);dqdx(i,j,3);];
            pqpy=[dqdy(i,j,1);dqdy(i,j,2);dqdy(i,j,3);];
            Q(i,j,:)=S-A*pqpx-B*pqpy;
        else
            r=(data.x(j)^2+data.y(i)^2)^0.5;
            cos=data.x(j)/r;sin=data.y(i)/r;
            vel=(u)*cos+(v)*sin+data.c;
           % vel=Uij*cos+Vij*sin;
            Q(i,j,:)=vel*(-cos*dqdx(i,j,:)-sin*dqdy(i,j,:)-q(i,j,:)/(2*r));
        end
    end
end
end
% function [Q]=makeQ_LPCE(q,data,tstep)
% boundi=4;ratio=1.67;rho0=1.0;
% Q=zeros(data.Ny,data.Nx,data.var);
% [dqdx,dqdy]=drp7all(q,data);
% %[dqdx,dqdy]=drp7_comp6bound(q,data);
% for i=1:1:data.Ny
%     for j=1:1:data.Nx
%         [Uij,Vij,Pij]=buildinc(data,i,j,tstep);
%         u=Uij+q(i,j,2);v=Vij+q(i,j,3);rho=1+q(i,j,1);%c=(1.0*(abs(Pij+q(i,j,4)))/rho)^0.5;
%         if (i>boundi&&i<data.Ny-boundi+1)&&(j>boundi&&j<data.Nx-boundi+1)
%             [S]=buildS_LPCE(q,i,j,data,tstep);
%             A=[Uij rho0 0 0;0 Uij 0 1.0/rho0;0 0 Uij 0;0 ratio*Pij 0 Uij;];
%             B=[Vij 0 rho0 0;0 Vij 0 0;0 0 Vij 1.0/rho0;0 0 ratio*Pij Vij;];
%             pqpx=[dqdx(i,j,1);dqdx(i,j,2);dqdx(i,j,3);dqdx(i,j,data.var)];
%             pqpy=[dqdy(i,j,1);dqdy(i,j,2);dqdy(i,j,3);dqdy(i,j,data.var)];
%             Q(i,j,:)=S-A*pqpx-B*pqpy;
%         else
%             r=(data.x(j)^2+data.y(i)^2)^0.5;
%             cos=data.x(j)/r;sin=data.y(i)/r;
%             %vel=u*cos+v*sin;
%             vel=data.c*sign(u*cos+v*sin)+u*cos+v*sin;
%             Q(i,j,:)=vel*(-cos*dqdx(i,j,:)-sin*dqdy(i,j,:)-q(i,j,:)/(2*r));
%         end
%     end
% end
% end
% function [dqdxbound,dqdybound]=drpboundary(q,data)
% 
% end
% function [dqdx,dqdy]=drp7(q,i,j,data)    
% a0=0;a1=0.77088238051822552;a2=-0.166705904414580469;a3=0.02084314277031176;
% dqdx=(a0*q(i,j,:)+a1*q(i,j+1,:)+a2*q(i,j+2,:)+a3*q(i,j+3,:)-a1*q(i,j-1,:)-a2*q(i,j-2,:)-a3*q(i,j-3,:))/data.dx;
% dqdy=(a0*q(i,j,:)+a1*q(i+1,j,:)+a2*q(i+2,j,:)+a3*q(i+3,j,:)-a1*q(i-1,j,:)-a2*q(i-2,j,:)-a3*q(i-3,j,:))/data.dx;
% end
function [dqdx,dqdy]=drp7all(q,data,dqdx,dqdy)

%%%dqdx%%%
for i=1:1:data.Ny
    for j=1:1:data.Nx
        if j>3&&j<data.Nx-2
            a0=0.0;a1=0.77088238051822552;a2=-0.166705904414580469;a3=0.02084314277031176;
           % a0=0;a1=0.7992664269741557;a2=-0.1894131415793246;a3=0.0265199520614978;
            %
            dqdx(i,j,:)=(a0*q(i,j,:)+a1*q(i,j+1,:)+a2*q(i,j+2,:)+a3*q(i,j+3,:)-a1*q(i,j-1,:)-a2*q(i,j-2,:)-a3*q(i,j-3,:))/data.dx;
        
        end
    end 
end
%%%dqdy%%%
for j=1:1:data.Nx
    for i=1:1:data.Ny
        if i>3&&i<data.Nx-2
            a0=0.0;a1=0.77088238051822552;a2=-0.166705904414580469;a3=0.02084314277031176;
            %a0=0;a1=0.7992664269741557;a2=-0.1894131415793246;a3=0.0265199520614978;
            dqdy(i,j,:)=(a0*q(i,j,:)+a1*q(i+1,j,:)+a2*q(i+2,j,:)+a3*q(i+3,j,:)-a1*q(i-1,j,:)-a2*q(i-2,j,:)-a3*q(i-3,j,:))/data.dy;
        
        end
    end
end

end

function [dqdx,dqdy]=drp7boundary(q,data)
dqdx=zeros(data.Ny,data.Nx,data.var);
dqdy=zeros(data.Ny,data.Nx,data.var);
%%%dqdx%%%
for i=1:1:data.Ny
        
           a0=0.0;a1=0.77088238051822552;a2=-0.166705904414580469;a3=0.02084314277031176;%
           %  a1=0.7992664269741557;a2=-0.1894131415793246;a3=0.0265199520614978;
            
            dqdx(i,4,:)=(a1*q(i,5,:)+a2*q(i,6,:)+a3*q(i,7,:)-a1*q(i,3,:)-a2*q(i,2,:)-a3*q(i,1,:))/data.dx;
            dqdx(i,data.Nx-3,:)=(a1*q(i,data.Nx-2,:)+a2*q(i,data.Nx-1,:)+a3*q(i,data.Nx,:)-a1*q(i,data.Nx-4,:)-a2*q(i,data.Nx-5,:)-a3*q(i,data.Nx-6,:))/data.dx;
            a06_0=-2.192280339;a06_1=4.748611401;a06_2=-5.108851915;a06_3=4.461567104;a06_4=-2.833498741;a06_5=1.128328861;a06_6=-0.203876371;
            dqdx(i,1,:)=(a06_0*q(i,1,:)+a06_1*q(i,2,:)+a06_2*q(i,3,:)+a06_3*q(i,4,:)+a06_4*q(i,5,:)+a06_5*q(i,6,:)+a06_6*q(i,7,:))/data.dx;
            a60_0=2.192280339;a60_1=-4.748611401;a60_2=5.108851915;a60_3=-4.461567104;a60_4=2.833498741;a60_5=-1.128328861;a60_6=0.203876371;
            dqdx(i,data.Nx,:)=(a60_0*q(i,data.Nx,:)+a60_1*q(i,data.Nx-1,:)+a60_2*q(i,data.Nx-2,:)+a60_3*q(i,data.Nx-3,:)+a60_4*q(i,data.Nx-4,:)+a60_5*q(i,data.Nx-5,:)+a60_6*q(i,data.Nx-6,:))/data.dx;
            a15f1=-0.209337622;a15_0=-1.084875676;a15_1=2.147776050;a15_2=-1.388928322;a15_3=0.768949766;a15_4=-0.281814650;a15_5=0.048230454;
            dqdx(i,2,:)=(a15f1*q(i,1,:)+a15_0*q(i,2,:)+a15_1*q(i,3,:)+a15_2*q(i,4,:)+a15_3*q(i,5,:)+a15_4*q(i,6,:)+a15_5*q(i,7,:))/data.dx;
            a51_1=0.209337622;a51_0=1.084875676;a51f1=-2.147776050;a51f2=1.388928322;a51f3=-0.768949766;a51f4=0.281814650;a51f5=-0.048230454;
            dqdx(i,data.Nx-1,:)=(a51_1*q(i,data.Nx,:)+a51_0*q(i,data.Nx-1,:)+a51f1*q(i,data.Nx-2,:)+a51f2*q(i,data.Nx-3,:)+a51f3*q(i,data.Nx-4,:)+a51f4*q(i,data.Nx-5,:)+a51f5*q(i,data.Nx-6,:))/data.dx;
            a24f2=0.049041958;a24f1=-0.468840357;a24_0=-0.474760914;a24_1=1.273274737;a24_2=-0.518484526;a24_3=0.166138533;a24_4=-0.026369431;
            dqdx(i,3,:)=(a24f2*q(i,1,:)+a24f1*q(i,2,:)+a24_0*q(i,3,:)+a24_1*q(i,4,:)+a24_2*q(i,5,:)+a24_3*q(i,6,:)+a24_4*q(i,7,:))/data.dx;
            a42_2=-0.049041958;a42_1=0.468840357;a42_0=0.474760914;a42f1=-1.273274737;a42f2=0.518484526;a42f3=-0.166138533;a42f4=0.026369431;
            dqdx(i,data.Nx-2,:)=(a42_2*q(i,data.Nx,:)+a42_1*q(i,data.Nx-1,:)+a42_0*q(i,data.Nx-2,:)+a42f1*q(i,data.Nx-3,:)+a42f2*q(i,data.Nx-4,:)+a42f3*q(i,data.Nx-5,:)+a42f4*q(i,data.Nx-6,:))/data.dx;
end
%%%dqdy%%%
for j=1:1:data.Nx
        
            a1=0.77088238051822552;a2=-0.166705904414580469;a3=0.02084314277031176;
            %a1=0.7992664269741557;a2=-0.1894131415793246;a3=0.0265199520614978;
            dqdy(4,j,:)=(a1*q(5,j,:)+a2*q(6,j,:)+a3*q(7,j,:)-a1*q(3,j,:)-a2*q(2,j,:)-a3*q(1,j,:))/data.dy;
            dqdy(data.Ny-3,j,:)=(a1*q(data.Ny-2,j,:)+a2*q(data.Ny-1,j,:)+a3*q(data.Ny,j,:)-a1*q(data.Ny-4,j,:)-a2*q(data.Ny-5,j,:)-a3*q(data.Ny-6,j,:))/data.dy;

            a06_0=-2.192280339;a06_1=4.748611401;a06_2=-5.108851915;a06_3=4.461567104;a06_4=-2.833498741;a06_5=1.128328861;a06_6=-0.203876371;
            dqdy(1,j,:)=(a06_0*q(1,j,:)+a06_1*q(2,j,:)+a06_2*q(3,j,:)+a06_3*q(4,j,:)+a06_4*q(5,j,:)+a06_5*q(6,j,:)+a06_6*q(7,j,:))/data.dy;
            a60_0=2.192280339;a60_1=-4.748611401;a60_2=5.108851915;a60_3=-4.461567104;a60_4=2.833498741;a60_5=-1.128328861;a60_6=0.203876371;
            dqdy(data.Ny,j,:)=(a60_0*q(data.Ny,j,:)+a60_1*q(data.Ny-1,j,:)+a60_2*q(data.Ny-2,j,:)+a60_3*q(data.Ny-3,j,:)+a60_4*q(data.Ny-4,j,:)+a60_5*q(data.Ny-5,j,:)+a60_6*q(data.Ny-6,j,:))/data.dy;
            a15f1=-0.209337622;a15_0=-1.084875676;a15_1=2.147776050;a15_2=-1.388928322;a15_3=0.768949766;a15_4=-0.281814650;a15_5=0.048230454;
            dqdy(2,j,:)=(a15f1*q(1,j,:)+a15_0*q(2,j,:)+a15_1*q(3,j,:)+a15_2*q(4,j,:)+a15_3*q(5,j,:)+a15_4*q(6,j,:)+a15_5*q(7,j,:))/data.dy;
            a51_1=0.209337622;a51_0=1.084875676;a51f1=-2.147776050;a51f2=1.388928322;a51f3=-0.768949766;a51f4=0.281814650;a51f5=-0.048230454;
            dqdy(data.Ny-1,j,:)=(a51_1*q(data.Ny,j,:)+a51_0*q(data.Ny-1,j,:)+a51f1*q(data.Ny-2,j,:)+a51f2*q(data.Ny-3,j,:)+a51f3*q(data.Ny-4,j,:)+a51f4*q(data.Ny-5,j,:)+a51f5*q(data.Ny-6,j,:))/data.dy;
            a24f2=0.049041958;a24f1=-0.468840357;a24_0=-0.474760914;a24_1=1.273274737;a24_2=-0.518484526;a24_3=0.166138533;a24_4=-0.026369431;
            dqdy(3,j,:)=(a24f2*q(1,j,:)+a24f1*q(2,j,:)+a24_0*q(3,j,:)+a24_1*q(4,j,:)+a24_2*q(5,j,:)+a24_3*q(6,j,:)+a24_4*q(7,j,:))/data.dy;
            a42_2=-0.049041958;a42_1=0.468840357;a42_0=0.474760914;a42f1=-1.273274737;a42f2=0.518484526;a42f3=-0.166138533;a42f4=0.026369431;
            dqdy(data.Ny-2,j,:)=(a42_2*q(data.Ny,j,:)+a42_1*q(data.Ny-1,j,:)+a42_0*q(data.Ny-2,j,:)+a42f1*q(data.Ny-3,j,:)+a42f2*q(data.Ny-4,j,:)+a42f3*q(data.Ny-5,j,:)+a42f4*q(data.Ny-6,j,:))/data.dy;
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

function [dqdx,dqdy]=compact4order(dqdx,dqdy,q,data)
lambda=0.3821038098462933;
a=0.7940346032820977;
b=0.0440346032820977;
A_x=zeros(data.Nx-8,data.Nx-8);
d_x=zeros(data.Nx-8,data.var);
A_y=zeros(data.Ny-8,data.Ny-8);
d_y=zeros(data.Ny-8,data.var);
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
for j=1:1:data.Ny-8
    if j==1
        A_y(1,1)=1;
        A_y(1,2)=lambda;
    elseif j==data.Ny-8
        A_y(data.Ny-8,data.Ny-8)=1;
        A_y(data.Ny-8,data.Ny-9)=lambda;    
    else
        A_y(j,j)=1;
        A_y(j,j+1)=lambda;
        A_y(j,j-1)=lambda;
    end
end

%%%dqdx%%%
for i=1:1:data.Ny
     for j=1:1:data.Nx-8
          if j==1
                d_x(1,:)=a*(q(i,6,:)-q(i,4,:))/data.dx+b*(q(i,7,:)-q(i,3,:))/data.dx-lambda*dqdx(i,4,:);
          elseif j==data.Nx-8
                d_x(data.Nx-8,:)=a*(q(i,data.Nx-3,:)-q(i,data.Nx-5,:))/data.dx+b*(q(i,data.Nx-2,:)-q(i,data.Nx-6,:))/data.dx-lambda*dqdx(i,data.Nx-3,:);          
          else
                d_x(j,:)=a*(q(i,j+5,:)-q(i,j+3,:))/data.dx+b*(q(i,j+6,:)-q(i,j+2,:))/data.dx;
          end
     end   
     for k=1:1:data.var
     x = tridiagonalSolver(A_x, d_x(:,k));
     dqdx(i,5:data.Nx-4,k)=x;
     end
end
%%%dqdy%%%
for j=1:1:data.Nx
    for i=1:1:data.Ny-8
          if i==1
             d_y(1,:)=a*(q(6,j,:)-q(4,j,:))/data.dy+b*(q(7,j,:)-q(3,j,:))/data.dy-lambda*dqdy(4,j,:);
          elseif i==data.Ny-8
              d_y(data.Ny-8,:)=a*(q(i+5,j,:)-q(i+3,j,:))/data.dy+b*(q(i+6,j,:)-q(i+2,j,:))/data.dy-lambda*dqdy(data.Ny-3,j,:);             
          else              
              d_y(i,:)=a*(q(i+5,j,:)-q(i+3,j,:))/data.dy+b*(q(i+6,j,:)-q(i+2,j,:))/data.dy;        
          end
      
    end
     for k=1:1:data.var
     y = tridiagonalSolver(A_y, d_y(:,k));
     dqdy(5:data.Ny-4,j,k)=y;
     end 
end

end

function [dqdx,dqdy]=compact6order(dqdx,dqdy,q,data)
lambda=0.4111403764203249;
a=0.7842616980320271;
b=0.0692748674241733;
c=-0.0038903521543495;
A_x=zeros(data.Nx-8,data.Nx-8);
d_x=zeros(data.Nx-8,data.var);
A_y=zeros(data.Ny-8,data.Ny-8);
d_y=zeros(data.Ny-8,data.var);
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
for j=1:1:data.Ny-8
    if j==1
        A_y(1,1)=1;
        A_y(1,2)=lambda;
    elseif j==data.Ny-8
        A_y(data.Ny-8,data.Ny-8)=1;
        A_y(data.Ny-8,data.Ny-9)=lambda;    
    else
        A_y(j,j)=1;
        A_y(j,j+1)=lambda;
        A_y(j,j-1)=lambda;
    end
end

%%%dqdx%%%
for i=1:1:data.Ny
     for j=1:1:data.Nx-8
          if j==1
                d_x(1,:)=a*(q(i,6,:)-q(i,4,:))/data.dx+b*(q(i,7,:)-q(i,3,:))/data.dx+c*(q(i,8,:)-q(i,2,:))/data.dx-lambda*dqdx(i,4,:);
          elseif j==data.Nx-8
                d_x(data.Nx-8,:)=a*(q(i,data.Nx-3,:)-q(i,data.Nx-5,:))/data.dx+b*(q(i,data.Nx-2,:)-q(i,data.Nx-6,:))/data.dx+c*(q(i,data.Nx-1,:)-q(i,data.Nx-7,:))/data.dx-lambda*dqdx(i,data.Nx-3,:);          
          else
                d_x(j,:)=a*(q(i,j+5,:)-q(i,j+3,:))/data.dx+b*(q(i,j+6,:)-q(i,j+2,:))/data.dx+c*(q(i,j+7,:)-q(i,j+1,:))/data.dx;
          end
     end   
     for k=1:1:data.var
     x = tridiagonalSolver(A_x, d_x(:,k));
     dqdx(i,5:data.Nx-4,k)=x;
     end
end
%%%dqdy%%%
for j=1:1:data.Nx
    for i=1:1:data.Ny-8
          if i==1
             d_y(1,:)=a*(q(6,j,:)-q(4,j,:))/data.dy+b*(q(7,j,:)-q(3,j,:))/data.dy+c*(q(8,j,:)-q(2,j,:))/data.dy-lambda*dqdy(4,j,:);
          elseif i==data.Ny-8
              d_y(data.Ny-8,:)=a*(q(i+5,j,:)-q(i+3,j,:))/data.dy+b*(q(i+6,j,:)-q(i+2,j,:))/data.dy+c*(q(i+7,j,:)-q(i+1,j,:))/data.dy-lambda*dqdy(data.Ny-3,j,:);             
          else              
              d_y(i,:)=a*(q(i+5,j,:)-q(i+3,j,:))/data.dy+b*(q(i+6,j,:)-q(i+2,j,:))/data.dy+c*(q(i+7,j,:)-q(i+1,j,:))/data.dy;        
          end
      
    end
     for k=1:1:data.var
     y = tridiagonalSolver(A_y, d_y(:,k));
     dqdy(5:data.Ny-4,j,k)=y;
     end 
end

end

function [S]=buildS(q,i,j,data,tstep)
ratio=1;
    if tstep==1
        dUdt=0.0;
        dVdt=0.0;
        dPdt=0.0;
    elseif tstep>1&&tstep<data.Nt
        [Uij1,Vij1,Pij1]=buildinc(data,i,j,tstep-1);
        [Uij2,Vij2,Pij2]=buildinc(data,i,j,tstep+1);
        
        dPdt=(Pij2-Pij1)/(2*data.dt);
        dUdt=(Uij2-Uij1)/(2*data.dt);
        dVdt=(Vij2-Vij1)/(2*data.dt);
    elseif tstep==data.Nt
        [Uij1,Vij1,Pij1]=buildinc(data,i,j,tstep-1);
        [Uij2,Vij2,Pij2]=buildinc(data,i,j,tstep);
        dPdt=(Pij2-Pij1)/(data.dt);
        dUdt=(Uij2-Uij1)/(data.dt);
        dVdt=(Vij2-Vij1)/(data.dt);
    end
     [Uij,Vij,~]=buildinc(data,i,j,tstep);
    u=Uij+q(i,j,2);
    v=Vij+q(i,j,3);
    [Uidown1j,Vidown1j,Pidown1j]=buildinc(data,i-1,j,tstep);
    [Uiup1j,Viup1j,Piup1j]=buildinc(data,i+1,j,tstep);
    
    dUdy=(Uiup1j-Uidown1j)/(2*data.dy);
    dVdy=(Viup1j-Vidown1j)/(2*data.dy);
    dPdy=(Piup1j-Pidown1j)/(2*data.dy);
    
    [Uijdown1,Vijdown1,Pijdown1]=buildinc(data,i,j-1,tstep);
    [Uijup1,Vijup1,Pijup1]=buildinc(data,i,j+1,tstep);
    
    dUdx=(Uijup1-Uijdown1)/(2*data.dx);
    dVdx=(Vijup1-Vijdown1)/(2*data.dx);
    dPdx=(Pijup1-Pijdown1)/(2*data.dx);
    rho=q(i,j,1)+1.0;
    %S4=-dPdt-ratio*(q(i,j,2)+Uij*q(i,j,1)/rho)*dPdx-ratio*(q(i,j,3)+Vij*q(i,j,1)/rho)*dPdy;
%     S4=-dPdt;
%     S2=-q(i,j,1)/rho*dUdt-(u-1/rho*Uij)*dUdx-(v-1/rho*Vij)*dUdy;
%     S3=-q(i,j,1)/rho*dVdt-(u-1/rho*Uij)*dVdx-(v-1/rho*Vij)*dVdy;
%     S4=-dPdt-ratio*(q(i,j,2)+q(i,j,1)/rho*Uij)*dPdx-ratio*(q(i,j,3)+q(i,j,1)/rho*Vij)*dPdy;
%     S2=-q(i,j,2)*dUdx-q(i,j,3)*dUdy;
%     S3=-q(i,j,2)*dVdx-q(i,j,3)*dVdy;  
%     S4=-dPdt-u*dPdx-v*dPdy;
    S1=dPdx;S2=dPdy;
     S=[-dPdt;0;0;];
 

end
% function [S]=buildS_LPCE(q,i,j,data,tstep)
%     [Uij,Vij,~]=buildinc(data,i,j,tstep);
%     u=Uij+q(i,j,2);
%     v=Vij+q(i,j,3);
%     [Uidown1j,Vidown1j,Pidown1j]=buildinc(data,i-1,j,tstep);
%     [Uiup1j,Viup1j,Piup1j]=buildinc(data,i+1,j,tstep);
%     dUdy=(Uiup1j-Uidown1j)/(2*data.dy);
%     dVdy=(Viup1j-Vidown1j)/(2*data.dy);
%     dPdy=(Piup1j-Pidown1j)/(2*data.dy);
%     [Uijdown1,Vijdown1,Pijdown1]=buildinc(data,i,j-1,tstep);
%     [Uijup1,Vijup1,Pijup1]=buildinc(data,i,j+1,tstep);
%     dUdx=(Uijup1-Uijdown1)/(2*data.dx);
%     dVdx=(Vijup1-Vijdown1)/(2*data.dx);
%     dPdx=(Pijup1-Pijdown1)/(2*data.dx);
%    
%     if tstep==1
%         dPdt=0.0;
%     elseif tstep>1&&tstep<data.Nt
%         [~,~,Pij1]=buildinc(data,i,j,tstep-1);
%         [~,~,Pij2]=buildinc(data,i,j,tstep+1);
%         dPdt=(Pij2-Pij1)/(2*data.dt);
%     elseif tstep==data.Nt
%         [~,~,Pij1]=buildinc(data,i,j,tstep-1);
%         [~,~,Pij2]=buildinc(data,i,j,tstep+1);
%         dPdt=(Pij2-Pij1)/(data.dt);
%     end
%     S2=-q(i,j,2)*dUdx-q(i,j,3)*dUdy;
%     S3=-q(i,j,2)*dVdx-q(i,j,3)*dVdy;  
%     S4=-dPdt-u*dPdx-v*dPdy;
%     S=[0;S2;S3;S4;];
% 
% end
function [qf]=filter(q,data)
qtmp=zeros(data.Ny,data.Nx,data.var);
xigema=0.75;
qghostright=0.5*(q(:,data.Nx,:)+q(:,data.Nx-1,:));
qghostleft=0.5*(q(:,1,:)+q(:,2,:));
qghosttop=0.5*(q(data.Nx,:,:)+q(data.Nx-1,:,:));
qghostbottom=0.5*(q(1,:,:)+q(2,:,:));
for i=1:1:data.Ny
    for j=1:1:data.Nx
if i>5&&i<data.Ny-4&&j>5&&j<data.Nx-4
    Dfx=xigema*(63/256*q(i,j,:)-105/512*q(i,j+1,:)-105/512*q(i,j-1,:)+15/128*q(i,j+2,:)+15/128*q(i,j-2,:)-45/1024*q(i,j+3,:)-45/1024*q(i,j-3,:)+5/512*q(i,j+4,:)+5/512*q(i,j-4,:)-1/1024*q(i,j+5,:)-1/1024*q(i,j-5,:));
    Dfy=xigema*(63/256*q(i,j,:)-105/512*q(i+1,j,:)-105/512*q(i-1,j,:)+15/128*q(i+2,j,:)+15/128*q(i-2,j,:)-45/1024*q(i+3,j,:)-45/1024*q(i-3,j,:)+5/512*q(i+4,j,:)+5/512*q(i-4,j,:)-1/1024*q(i+5,j,:)-1/1024*q(i-5,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (i==5&&j>4&&j<data.Nx-3)||(i==data.Ny-4&&j>4&&j<data.Nx-3)||(j==5&&i>4&&i<data.Nx-3)||(j==data.Ny-4&&i>4&&i<data.Nx-3)
    Dfx=xigema*(35/128*q(i,j,:)-7/32*q(i,j+1,:)-7/32*q(i,j-1,:)+7/64*q(i,j+2,:)+7/64*q(i,j-2,:)-1/32*q(i,j+3,:)-1/32*q(i,j-3,:)+1/256*q(i,j+4,:)+1/256*q(i,j-4,:));
    Dfy=xigema*(35/128*q(i,j,:)-7/32*q(i+1,j,:)-7/32*q(i-1,j,:)+7/64*q(i+2,j,:)+7/64*q(i-2,j,:)-1/32*q(i+3,j,:)-1/32*q(i-3,j,:)+1/256*q(i+4,j,:)+1/256*q(i-4,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (i==4&&j>3&&j<data.Nx-2)||(i==data.Ny-3&&j>3&&j<data.Nx-2)||(j==4&&i>3&&i<data.Nx-2)||(j==data.Ny-3&&i>3&&i<data.Nx-2)
    Dfx=xigema*(5/16*q(i,j,:)-15/64*q(i,j+1,:)-15/64*q(i,j-1,:)+3/32*q(i,j+2,:)+3/32*q(i,j-2,:)-1.0/64*q(i,j+3,:)-1.0/64*q(i,j-3,:));
    Dfy=xigema*(5/16*q(i,j,:)-15/64*q(i+1,j,:)-15/64*q(i-1,j,:)+3/32*q(i+2,j,:)+3/32*q(i-2,j,:)-1.0/64*q(i+3,j,:)-1.0/64*q(i-3,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (i==3&&j>2&&j<data.Nx-1)||(i==data.Ny-2&&j>2&&j<data.Nx-1)||(j==3&&i>2&&i<data.Nx-1)||(j==data.Ny-2&&i>2&&i<data.Nx-1)
    Dfx=xigema*(3/8*q(i,j,:)-1/4*q(i,j+1,:)-1/4*q(i,j-1,:)+1/6*q(i,j+2,:)+1/6*q(i,j-2,:));
    Dfy=xigema*(3/8*q(i,j,:)-1/4*q(i+1,j,:)-1/4*q(i-1,j,:)+1/6*q(i+2,j,:)+1/6*q(i-2,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (i==2&&j>1&&j<data.Nx)||(i==data.Ny-1&&j>1&&j<data.Nx)||(j==2&&i>1&&i<data.Nx)||(j==data.Ny-1&&i>1&&i<data.Nx)
    Dfx=xigema*(1/2*q(i,j,:)-1/4*q(i,j+1,:)-1/4*q(i,j-1,:));
    Dfy=xigema*(1/2*q(i,j,:)-1/4*q(i+1,j,:)-1/4*q(i-1,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (i==data.Ny&&j>1&&j<data.Nx) %top
    Dfx=xigema*(1/2*q(i,j,:)-1/4*q(i,j+1,:)-1/4*q(i,j-1,:));
    Dfy=xigema*(1/2*q(i,j,:)-1/4*qghosttop(1,j,:)-1/4*q(i-1,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (i==1&&j>1&&j<data.Nx) %bottom
    Dfx=xigema*(1/2*q(i,j,:)-1/4*q(i,j+1,:)-1/4*q(i,j-1,:));
    Dfy=xigema*(1/2*q(i,j,:)-1/4*qghostbottom(1,j,:)-1/4*q(i+1,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (j==1&&i>1&&i<data.Ny) %left
    Dfx=xigema*(1/2*q(i,j,:)-1/4*qghostleft(i,1,:)-1/4*q(i,j+1,:));
    Dfy=xigema*(1/2*q(i,j,:)-1/4*q(i-1,j,:)-1/4*q(i+1,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (j==data.Nx&&i>1&&i<data.Ny) %right
    Dfx=xigema*(1/2*q(i,j,:)-1/4*qghostright(i,1,:)-1/4*q(i,j-1,:));
    Dfy=xigema*(1/2*q(i,j,:)-1/4*q(i-1,j,:)-1/4*q(i+1,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (i==1&&j==1)
    Dfx=xigema*(1/2*q(i,j,:)-1/4*qghostleft(i,1,:)-1/4*q(i,j+1,:));
    Dfy=xigema*(1/2*q(i,j,:)-1/4*qghostbottom(1,j,:)-1/4*q(i+1,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (i==1&&j==data.Nx) 
    Dfx=xigema*(1/2*q(i,j,:)-1/4*qghostright(i,1,:)-1/4*q(i,j-1,:));
    Dfy=xigema*(1/2*q(i,j,:)-1/4*qghostbottom(1,j,:)-1/4*q(i+1,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (i==data.Ny&&j==1) 
    Dfx=xigema*(1/2*q(i,j,:)-1/4*qghostleft(i,1,:)-1/4*q(i,j+1,:));
    Dfy=xigema*(1/2*q(i,j,:)-1/4*qghosttop(1,j,:)-1/4*q(i-1,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (i==data.Ny&&j==data.Nx) 
    Dfx=xigema*(1/2*q(i,j,:)-1/4*qghostright(i,1,:)-1/4*q(i,j-1,:));
    Dfy=xigema*(1/2*q(i,j,:)-1/4*qghosttop(1,j,:)-1/4*q(i-1,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
end   
    end
end
qf=qtmp;
end

function [qf]=filter_opt(q,data)
qtmp=zeros(data.Ny,data.Nx,data.var);
xigema=0.25;
qghostright=0.5*(q(:,data.Nx,:)+q(:,data.Nx-1,:));
qghostleft=0.5*(q(:,1,:)+q(:,2,:));
qghosttop=0.5*(q(data.Nx,:,:)+q(data.Nx-1,:,:));
qghostbottom=0.5*(q(1,:,:)+q(2,:,:));
o11p0=0.215044884112;   o9p0=0.243527493120;
o11p1=-0.187772883589;  o9p1=-0.204788880640;
o11p2=0.123755948787;   o9p2=0.120007591680;
o11p3=-0.059227575576;  o9p3=-0.045211119360;
o11p4=0.018721609157;   o9p4=0.008228661760;
o11p5=-0.002999540835;
for i=1:1:data.Ny
    for j=1:1:data.Nx
if i>5&&i<data.Ny-4&&j>5&&j<data.Nx-4
    Dfx=xigema*(o11p0*q(i,j,:)+o11p1*q(i,j+1,:)+o11p1*q(i,j-1,:)+o11p2*q(i,j+2,:)+o11p2*q(i,j-2,:)+o11p3*q(i,j+3,:)+o11p3*q(i,j-3,:)+o11p4*q(i,j+4,:)+o11p4*q(i,j-4,:)+o11p5*q(i,j+5,:)+o11p5*q(i,j-5,:));
    Dfy=xigema*(o11p0*q(i,j,:)+o11p1*q(i+1,j,:)+o11p1*q(i-1,j,:)+o11p2*q(i+2,j,:)+o11p2*q(i-2,j,:)+o11p3*q(i+3,j,:)+o11p3*q(i-3,j,:)+o11p4*q(i+4,j,:)+o11p4*q(i-4,j,:)+o11p5*q(i+5,j,:)+o11p5*q(i-5,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (i==5&&j>4&&j<data.Nx-3)||(i==data.Ny-4&&j>4&&j<data.Nx-3)||(j==5&&i>4&&i<data.Nx-3)||(j==data.Ny-4&&i>4&&i<data.Nx-3)
    Dfx=xigema*(o9p0*q(i,j,:)+o9p1*q(i,j+1,:)+o9p1*q(i,j-1,:)+o9p2*q(i,j+2,:)+o9p2*q(i,j-2,:)+o9p3*q(i,j+3,:)+o9p3*q(i,j-3,:)+o9p4*q(i,j+4,:)+o9p4*q(i,j-4,:));
    Dfy=xigema*(o9p0*q(i,j,:)+o9p1*q(i+1,j,:)+o9p1*q(i-1,j,:)+o9p2*q(i+2,j,:)+o9p2*q(i-2,j,:)+o9p3*q(i+3,j,:)+o9p3*q(i-3,j,:)+o9p4*q(i+4,j,:)+o9p4*q(i-4,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (i==4&&j>3&&j<data.Nx-2)||(i==data.Ny-3&&j>3&&j<data.Nx-2)||(j==4&&i>3&&i<data.Nx-2)||(j==data.Ny-3&&i>3&&i<data.Nx-2)
    Dfx=xigema*(5/16*q(i,j,:)-15/64*q(i,j+1,:)-15/64*q(i,j-1,:)+3/32*q(i,j+2,:)+3/32*q(i,j-2,:)-1.0/64*q(i,j+3,:)-1.0/64*q(i,j-3,:));
    Dfy=xigema*(5/16*q(i,j,:)-15/64*q(i+1,j,:)-15/64*q(i-1,j,:)+3/32*q(i+2,j,:)+3/32*q(i-2,j,:)-1.0/64*q(i+3,j,:)-1.0/64*q(i-3,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (i==3&&j>2&&j<data.Nx-1)||(i==data.Ny-2&&j>2&&j<data.Nx-1)||(j==3&&i>2&&i<data.Nx-1)||(j==data.Ny-2&&i>2&&i<data.Nx-1)
    Dfx=xigema*(3/8*q(i,j,:)-1/4*q(i,j+1,:)-1/4*q(i,j-1,:)+1/6*q(i,j+2,:)+1/6*q(i,j-2,:));
    Dfy=xigema*(3/8*q(i,j,:)-1/4*q(i+1,j,:)-1/4*q(i-1,j,:)+1/6*q(i+2,j,:)+1/6*q(i-2,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (i==2&&j>1&&j<data.Nx)||(i==data.Ny-1&&j>1&&j<data.Nx)||(j==2&&i>1&&i<data.Nx)||(j==data.Ny-1&&i>1&&i<data.Nx)
    Dfx=xigema*(1/2*q(i,j,:)-1/4*q(i,j+1,:)-1/4*q(i,j-1,:));
    Dfy=xigema*(1/2*q(i,j,:)-1/4*q(i+1,j,:)-1/4*q(i-1,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (i==data.Ny&&j>1&&j<data.Nx) %top
    Dfx=xigema*(1/2*q(i,j,:)-1/4*q(i,j+1,:)-1/4*q(i,j-1,:));
    Dfy=xigema*(1/2*q(i,j,:)-1/4*qghosttop(1,j,:)-1/4*q(i-1,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (i==1&&j>1&&j<data.Nx) %bottom
    Dfx=xigema*(1/2*q(i,j,:)-1/4*q(i,j+1,:)-1/4*q(i,j-1,:));
    Dfy=xigema*(1/2*q(i,j,:)-1/4*qghostbottom(1,j,:)-1/4*q(i+1,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (j==1&&i>1&&i<data.Ny) %left
    Dfx=xigema*(1/2*q(i,j,:)-1/4*qghostleft(i,1,:)-1/4*q(i,j+1,:));
    Dfy=xigema*(1/2*q(i,j,:)-1/4*q(i-1,j,:)-1/4*q(i+1,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (j==data.Nx&&i>1&&i<data.Ny) %right
    Dfx=xigema*(1/2*q(i,j,:)-1/4*qghostright(i,1,:)-1/4*q(i,j-1,:));
    Dfy=xigema*(1/2*q(i,j,:)-1/4*q(i-1,j,:)-1/4*q(i+1,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (i==1&&j==1)
    Dfx=xigema*(1/2*q(i,j,:)-1/4*qghostleft(i,1,:)-1/4*q(i,j+1,:));
    Dfy=xigema*(1/2*q(i,j,:)-1/4*qghostbottom(1,j,:)-1/4*q(i+1,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (i==1&&j==data.Nx) 
    Dfx=xigema*(1/2*q(i,j,:)-1/4*qghostright(i,1,:)-1/4*q(i,j-1,:));
    Dfy=xigema*(1/2*q(i,j,:)-1/4*qghostbottom(1,j,:)-1/4*q(i+1,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (i==data.Ny&&j==1) 
    Dfx=xigema*(1/2*q(i,j,:)-1/4*qghostleft(i,1,:)-1/4*q(i,j+1,:));
    Dfy=xigema*(1/2*q(i,j,:)-1/4*qghosttop(1,j,:)-1/4*q(i-1,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
elseif (i==data.Ny&&j==data.Nx) 
    Dfx=xigema*(1/2*q(i,j,:)-1/4*qghostright(i,1,:)-1/4*q(i,j-1,:));
    Dfy=xigema*(1/2*q(i,j,:)-1/4*qghosttop(1,j,:)-1/4*q(i-1,j,:));
    qtmp(i,j,:)=q(i,j,:)-0.5*Dfx-0.5*Dfy;
end   
    end
end
qf=qtmp;
end

function [qf]=filter_1(q,data)
qtmp=zeros(data.Ny,data.Nx,data.var);
qtmp2=zeros(data.Ny,data.Nx,data.var);
xigema=0.75;
qghostright=0.5*(q(:,data.Nx,:)+q(:,data.Nx-1,:));
qghostleft=0.5*(q(:,1,:)+q(:,2,:));
qghosttop=0.5*(q(data.Nx,:,:)+q(data.Nx-1,:,:));
qghostbottom=0.5*(q(1,:,:)+q(2,:,:));
for i=1:1:data.Ny
    for j=1:1:data.Nx
if j>5&&j<data.Nx-4
    Dfx=xigema*(0.215044884112*q(i,j,:)-0.187772883589*q(i,j+1,:)-0.187772883589*q(i,j-1,:)+0.123755948787*q(i,j+2,:)+0.123755948787*q(i,j-2,:)-0.059227575576*q(i,j+3,:)-0.059227575576*q(i,j-3,:)+0.018721609157*q(i,j+4,:)+0.018721609157*q(i,j-4,:)-0.002999540835*q(i,j+5,:)-0.002999540835*q(i,j-5,:));
%    Dfy=xigema*(o11p0*q(i,j,:)+o11p1*q(i+1,j,:)+o11p1*q(i-1,j,:)+o11p2*q(i+2,j,:)+o11p2*q(i-2,j,:)+o11p3*q(i+3,j,:)+o11p3*q(i-3,j,:)+o11p4*q(i+4,j,:)+o11p4*q(i-4,j,:)+o11p5*q(i+5,j,:)+o11p5*q(i-5,j,:));
    qtmp(i,j,:)=q(i,j,:)-Dfx;
elseif (j==5)||(j==data.Ny-4)
    Dfx=xigema*(0.243527493120*q(i,j,:)-0.204788880640*q(i,j+1,:)-0.204788880640*q(i,j-1,:)+0.120007591680*q(i,j+2,:)+0.120007591680*q(i,j-2,:)-0.045211119360*q(i,j+3,:)-0.045211119360*q(i,j-3,:)+0.00822866176*q(i,j+4,:)+0.00822866176*q(i,j-4,:));
%    Dfy=xigema*(o9p0*q(i,j,:)+o9p1*q(i+1,j,:)+o9p1*q(i-1,j,:)+o9p2*q(i+2,j,:)+o9p2*q(i-2,j,:)+o9p3*q(i+3,j,:)+o9p3*q(i-3,j,:)+o9p4*q(i+4,j,:)+o9p4*q(i-4,j,:));
    qtmp(i,j,:)=q(i,j,:)-Dfx;
elseif (j==4)||(j==data.Ny-3)
    Dfx=xigema*(5/16*q(i,j,:)-15/64*q(i,j+1,:)-15/64*q(i,j-1,:)+3/32*q(i,j+2,:)+3/32*q(i,j-2,:)-1.0/64*q(i,j+3,:)-1.0/64*q(i,j-3,:));
%    Dfy=xigema*(5/16*q(i,j,:)-15/64*q(i+1,j,:)-15/64*q(i-1,j,:)+3/32*q(i+2,j,:)+3/32*q(i-2,j,:)-1.0/64*q(i+3,j,:)-1.0/64*q(i-3,j,:));
    qtmp(i,j,:)=q(i,j,:)-Dfx;
elseif (j==3)||(j==data.Ny-2)
    Dfx=xigema*(3/8*q(i,j,:)-1/4*q(i,j+1,:)-1/4*q(i,j-1,:)+1/6*q(i,j+2,:)+1/6*q(i,j-2,:));
%    Dfy=xigema*(3/8*q(i,j,:)-1/4*q(i+1,j,:)-1/4*q(i-1,j,:)+1/6*q(i+2,j,:)+1/6*q(i-2,j,:));
    qtmp(i,j,:)=q(i,j,:)-Dfx;
elseif (j==2)||(j==data.Ny-1)
    Dfx=xigema*(1/2*q(i,j,:)-1/4*q(i,j+1,:)-1/4*q(i,j-1,:));
%    Dfy=xigema*(1/2*q(i,j,:)-1/4*q(i+1,j,:)-1/4*q(i-1,j,:));
    qtmp(i,j,:)=q(i,j,:)-Dfx;
elseif (j==1) %left
    Dfx=xigema*(1/2*q(i,j,:)-1/4*qghostleft(i,1,:)-1/4*q(i,j+1,:));
%    Dfy=xigema*(1/2*q(i,j,:)-1/4*q(i-1,j,:)-1/4*q(i+1,j,:));
    qtmp(i,j,:)=q(i,j,:)-Dfx;
elseif (j==data.Nx) %right
    Dfx=xigema*(1/2*q(i,j,:)-1/4*qghostright(i,1,:)-1/4*q(i,j-1,:));
%    Dfy=xigema*(1/2*q(i,j,:)-1/4*q(i-1,j,:)-1/4*q(i+1,j,:));
    qtmp(i,j,:)=q(i,j,:)-Dfx;
end   
    end
end

for i=1:1:data.Ny
    for j=1:1:data.Nx
if i>5&&i<data.Ny-4
   % Dfx=xigema*(o11p0*q(i,j,:)+o11p1*q(i,j+1,:)+o11p1*q(i,j-1,:)+o11p2*q(i,j+2,:)+o11p2*q(i,j-2,:)+o11p3*q(i,j+3,:)+o11p3*q(i,j-3,:)+o11p4*q(i,j+4,:)+o11p4*q(i,j-4,:)+o11p5*q(i,j+5,:)+o11p5*q(i,j-5,:));
    Dfy=xigema*(0.215044884112*qtmp(i,j,:)-0.187772883589*qtmp(i+1,j,:)-0.187772883589*qtmp(i-1,j,:)+0.123755948787*qtmp(i+2,j,:)+0.123755948787*qtmp(i-2,j,:)-0.059227575576*qtmp(i+3,j,:)-0.059227575576*qtmp(i-3,j,:)+0.018721609157*qtmp(i+4,j,:)+0.018721609157*qtmp(i-4,j,:)-0.002999540835*qtmp(i+5,j,:)-0.002999540835*qtmp(i-5,j,:));
    qtmp2(i,j,:)=qtmp(i,j,:)-Dfy;
elseif (i==5)||(i==data.Ny-4)
  %  Dfx=xigema*(o9p0*q(i,j,:)+o9p1*q(i,j+1,:)+o9p1*q(i,j-1,:)+o9p2*q(i,j+2,:)+o9p2*q(i,j-2,:)+o9p3*q(i,j+3,:)+o9p3*q(i,j-3,:)+o9p4*q(i,j+4,:)+o9p4*q(i,j-4,:));
    Dfy=xigema*(0.243527493120*qtmp(i,j,:)-0.204788880640*qtmp(i+1,j,:)-0.204788880640*qtmp(i-1,j,:)+0.120007591680*qtmp(i+2,j,:)+0.120007591680*qtmp(i-2,j,:)-0.045211119360*qtmp(i+3,j,:)-0.045211119360*qtmp(i-3,j,:)+0.00822866176*qtmp(i+4,j,:)+0.00822866176*qtmp(i-4,j,:));
    qtmp2(i,j,:)=qtmp(i,j,:)-Dfy;
elseif (i==4)||(i==data.Ny-3)
   % Dfx=xigema*(5/16*q(i,j,:)-15/64*q(i,j+1,:)-15/64*q(i,j-1,:)+3/32*q(i,j+2,:)+3/32*q(i,j-2,:)-1.0/64*q(i,j+3,:)-1.0/64*q(i,j-3,:));
    Dfy=xigema*(5/16*qtmp(i,j,:)-15/64*qtmp(i+1,j,:)-15/64*qtmp(i-1,j,:)+3/32*qtmp(i+2,j,:)+3/32*qtmp(i-2,j,:)-1.0/64*qtmp(i+3,j,:)-1.0/64*qtmp(i-3,j,:));
    qtmp2(i,j,:)=qtmp(i,j,:)-Dfy;
elseif (i==3)||(i==data.Ny-2)
%    Dfx=xigema*(3/8*q(i,j,:)-1/4*q(i,j+1,:)-1/4*q(i,j-1,:)+1/6*q(i,j+2,:)+1/6*q(i,j-2,:));
    Dfy=xigema*(3/8*qtmp(i,j,:)-1/4*qtmp(i+1,j,:)-1/4*qtmp(i-1,j,:)+1/6*qtmp(i+2,j,:)+1/6*qtmp(i-2,j,:));
    qtmp2(i,j,:)=qtmp(i,j,:)-Dfy;
elseif (i==2)||(i==data.Ny-1)
%    Dfx=xigema*(1/2*q(i,j,:)-1/4*q(i,j+1,:)-1/4*q(i,j-1,:));
    Dfy=xigema*(1/2*qtmp(i,j,:)-1/4*qtmp(i+1,j,:)-1/4*qtmp(i-1,j,:));
    qtmp2(i,j,:)=qtmp(i,j,:)-Dfy;
elseif (i==data.Ny) %top
%    Dfx=xigema*(1/2*q(i,j,:)-1/4*q(i,j+1,:)-1/4*q(i,j-1,:));
    Dfy=xigema*(1/2*qtmp(i,j,:)-1/4*qghosttop(1,j,:)-1/4*qtmp(i-1,j,:));
    qtmp2(i,j,:)=qtmp(i,j,:)-Dfy;
elseif (i==1) %bottom
 %   Dfx=xigema*(1/2*q(i,j,:)-1/4*q(i,j+1,:)-1/4*q(i,j-1,:));
    Dfy=xigema*(1/2*qtmp(i,j,:)-1/4*qghostbottom(1,j,:)-1/4*qtmp(i+1,j,:));
    qtmp2(i,j,:)=qtmp(i,j,:)-Dfy;
end   
    end
end
qf=qtmp2;
end

function [qf]=filter_x(q,data)
qtmp=zeros(data.Ny,data.Nx,data.var);
xigema=0.1;
qghostright=0.5*(q(:,data.Nx,:)+q(:,data.Nx-1,:));
qghostleft=0.5*(q(:,1,:)+q(:,2,:));
qghosttop=0.5*(q(data.Nx,:,:)+q(data.Nx-1,:,:));
qghostbottom=0.5*(q(1,:,:)+q(2,:,:));
for i=1:1:data.Ny
    for j=1:1:data.Nx
if i>=6&&i<=data.Ny-5&&j>=6&&j<=data.Nx-5
    Df1=xigema*(63/256*q(i,j,:)-105/512*q(i+1,j+1,:)-105/512*q(i-1,j-1,:)+15/128*q(i+2,j+2,:)+15/128*q(i-2,j-2,:)-45/1024*q(i+3,j+3,:)-45/1024*q(i-3,j-3,:)+5/512*q(i+4,j+4,:)+5/512*q(i-4,j-4,:)-1/1024*q(i+5,j+5,:)-1/1024*q(i-5,j-5,:));
    Df2=xigema*(63/256*q(i,j,:)-105/512*q(i-1,j+1,:)-105/512*q(i+1,j-1,:)+15/128*q(i+2,j-2,:)+15/128*q(i-2,j+2,:)-45/1024*q(i+3,j-3,:)-45/1024*q(i-3,j+3,:)+5/512*q(i+4,j-4,:)+5/512*q(i-4,j+4,:)-1/1024*q(i+5,j-5,:)-1/1024*q(i-5,j+5,:));
    qtmp(i,j,:)=q(i,j,:)-Df1*0.5-Df2*0.5;
elseif (i==5&&j>4&&j<data.Nx-3)||(i==data.Ny-4&&j>4&&j<data.Nx-3)||(j==5&&i>4&&i<data.Nx-3)||(j==data.Ny-4&&i>4&&i<data.Nx-3)
    Df1=xigema*(35/128*q(i,j,:)-7/32*q(i+1,j+1,:)-7/32*q(i-1,j-1,:)+7/64*q(i+2,j+2,:)+7/64*q(i-2,j-2,:)-1/32*q(i+3,j+3,:)-1/32*q(i-3,j-3,:)+1/256*q(i+4,j+4,:)+1/256*q(i-4,j-4,:));
    Df2=xigema*(35/128*q(i,j,:)-7/32*q(i+1,j-1,:)-7/32*q(i-1,j+1,:)+7/64*q(i+2,j-2,:)+7/64*q(i-2,j+2,:)-1/32*q(i+3,j-3,:)-1/32*q(i-3,j+3,:)+1/256*q(i+4,j-4,:)+1/256*q(i-4,j+4,:));
    qtmp(i,j,:)=q(i,j,:)-Df1*0.5-Df2*0.5;
elseif (i==4&&j>3&&j<data.Nx-2)||(i==data.Ny-3&&j>3&&j<data.Nx-2)||(j==4&&i>3&&i<data.Nx-2)||(j==data.Ny-3&&i>3&&i<data.Nx-2)
    Df1=xigema*(5/16*q(i,j,:)-15/64*q(i+1,j+1,:)-15/64*q(i-1,j-1,:)+3/32*q(i+2,j+2,:)+3/32*q(i-2,j-2,:)-1.0/64*q(i+3,j+3,:)-1.0/64*q(i-3,j-3,:));
    Df2=xigema*(5/16*q(i,j,:)-15/64*q(i+1,j-1,:)-15/64*q(i-1,j+1,:)+3/32*q(i+2,j-2,:)+3/32*q(i-2,j+2,:)-1.0/64*q(i+3,j-3,:)-1.0/64*q(i-3,j+3,:));
    qtmp(i,j,:)=q(i,j,:)-Df1*0.5-Df2*0.5;
elseif (i==3&&j>2&&j<data.Nx-1)||(i==data.Ny-2&&j>2&&j<data.Nx-1)||(j==3&&i>2&&i<data.Nx-1)||(j==data.Ny-2&&i>2&&i<data.Nx-1)
    Df1=xigema*(3/8*q(i,j,:)-1/4*q(i+1,j+1,:)-1/4*q(i-1,j-1,:)+1/6*q(i+2,j+2,:)+1/6*q(i-2,j-2,:));
    Df2=xigema*(3/8*q(i,j,:)-1/4*q(i+1,j-1,:)-1/4*q(i-1,j+1,:)+1/6*q(i+2,j-2,:)+1/6*q(i-2,j+2,:));
    qtmp(i,j,:)=q(i,j,:)-Df1*0.5-Df2*0.5;
elseif (i==2&&j>1&&j<data.Nx)||(i==data.Ny-1&&j>1&&j<data.Nx)||(j==2&&i>1&&i<data.Nx)||(j==data.Ny-1&&i>1&&i<data.Nx)
    Df1=xigema*(1/2*q(i,j,:)-1/4*q(i+1,j+1,:)-1/4*q(i-1,j-1,:));
    Df2=xigema*(1/2*q(i,j,:)-1/4*q(i+1,j-1,:)-1/4*q(i-1,j+1,:));
    qtmp(i,j,:)=q(i,j,:)-Df1*0.5-Df2*0.5;
elseif (i==data.Ny&&j>1&&j<data.Nx) %top
    Dfx=xigema*(1/2*q(i,j,:)-1/4*q(i,j+1,:)-1/4*q(i,j-1,:));
    %Dfy=xigema*(1/2*q(i,j,:)-1/4*qghosttop(1,j,:)-1/4*q(i-1,j,:));
    qtmp(i,j,:)=q(i,j,:)-Dfx;
elseif (i==1&&j>1&&j<data.Nx) %bottom
    Dfx=xigema*(1/2*q(i,j,:)-1/4*q(i,j+1,:)-1/4*q(i,j-1,:));
    %Dfy=xigema*(1/2*q(i,j,:)-1/4*qghostbottom(1,j,:)-1/4*q(i+1,j,:));
    qtmp(i,j,:)=q(i,j,:)-Dfx;
elseif (j==1&&i>1&&i<data.Ny) %left
    %Dfx=xigema*(1/2*q(i,j,:)-1/4*qghostleft(i,1,:)-1/4*q(i,j+1,:));
    Dfy=xigema*(1/2*q(i,j,:)-1/4*q(i-1,j,:)-1/4*q(i+1,j,:));
    qtmp(i,j,:)=q(i,j,:)-Dfy;
elseif (j==data.Nx&&i>1&&i<data.Ny) %right
    %Dfx=xigema*(1/2*q(i,j,:)-1/4*qghostright(i,1,:)-1/4*q(i,j-1,:));
    Dfy=xigema*(1/2*q(i,j,:)-1/4*q(i-1,j,:)-1/4*q(i+1,j,:));
    qtmp(i,j,:)=q(i,j,:)-Dfy;
elseif (i==1&&j==1)
    %Dfx=xigema*(1/2*q(i,j,:)-1/4*qghostleft(i,1,:)-1/4*q(i,j+1,:));
    %Dfy=xigema*(1/2*q(i,j,:)-1/4*qghostbottom(1,j,:)-1/4*q(i+1,j,:));
    qtmp(i,j,:)=q(i,j,:);
elseif (i==1&&j==data.Nx) 
%     Dfx=xigema*(1/2*q(i,j,:)-1/4*qghostright(i,1,:)-1/4*q(i,j-1,:));
%     Dfy=xigema*(1/2*q(i,j,:)-1/4*qghostbottom(1,j,:)-1/4*q(i+1,j,:));
    qtmp(i,j,:)=q(i,j,:);
elseif (i==data.Ny&&j==1) 
%     Dfx=xigema*(1/2*q(i,j,:)-1/4*qghostleft(i,1,:)-1/4*q(i,j+1,:));
%     Dfy=xigema*(1/2*q(i,j,:)-1/4*qghosttop(1,j,:)-1/4*q(i-1,j,:));
    qtmp(i,j,:)=q(i,j,:);
elseif (i==data.Ny&&j==data.Nx) 
%     Dfx=xigema*(1/2*q(i,j,:)-1/4*qghostright(i,1,:)-1/4*q(i,j-1,:));
%     Dfy=xigema*(1/2*q(i,j,:)-1/4*qghosttop(1,j,:)-1/4*q(i-1,j,:));
    qtmp(i,j,:)=q(i,j,:);
end   
    end
end
qf=qtmp;
end

function [q]=initial(data)
q=zeros(data.Ny,data.Nx,data.var);
R=287.053;
T=273+20;
% for i=1:1:data.Ny
%     for j=1:1:data.Nx
%         %         %q(i,j,data.var)=exp(-log(2)/9*((data.x(j))^2+data.y(i)^2))+exp(-log(2)/9*((data.x(j)-40)^2+data.y(i)^2))+exp(-log(2)/9*((data.x(j)+40)^2+data.y(i)^2));
% 
%         q(i,j,1)=exp(-log(2)/9*((data.x(j))^2+data.y(i)^2));
%         
%         
%     end
% end
end
function [Uij,Vij,Pij]=buildinc(data,i,j,tstep)
[incleft,incright,inctop,incbottom]=incregion(data);
gamma=data.gamma;
r0=data.r0;
w=gamma/(4*pi*r0^2);rho0=1.0;P0=rho0*data.c^2/1;
if (i>incbottom&&i<inctop)&&(j>incleft&&j<incright)
% % %if (data.y(i)^2+data.x(j)^2<=625)
% %     
    %%%vortex pair%%%
    t=tstep*data.dt;%theta=atan(data.y(i)/data.x(j));
%     rc=0.1;
%     r=(data.x(j)^2+data.y(i)^2)^0.5;
%     veltheta=gamma*r/(2*pi*(rc^2+r^2));
%     Uij=-veltheta*data.x(j)/r;
%     Vij=veltheta*data.y(i)/r;
    z=data.x(j)+1i*data.y(i);b=r0*exp(1i*w*t);
%     RR=real(b^2/(z^2-b^2));
%    Pij=P0-0.5*rho0*(Uij^2+Vij^2)+rho0*gamma*w/pi*RR;%;
%     Pij=P0+rho0*gamma*w/pi*real(RR-gamma*r/(2*pi*w*sqrt(z^2-b^2)));%;

    Uij=real(gamma/(1i*pi)*z/(z^2-b^2));
    Vij=-imag(gamma/(1i*pi)*z/(z^2-b^2));
    Pij=P0+rho0*gamma*w/pi*real(b^2/(z^2-b^2))-0.5*rho0*(Uij^2+Vij^2);
%============================================%
%    Uij=0.0;
%    %  Uij=0.5*data.c;
%     Vij=0.0;
%     Pij=P0;
% %============================================%
 %   Uij=0.5*data.c*(data.y(i)/abs(data.y(i)));
%    Vij=0.0;
%    Pij=0.0;
else
     %        Uij=0.5*data.c;
       Uij=0.0;
%   %   Uij=0.5*data.c*(data.y(i)/abs(data.y(i)));
        Vij=0.0;
        Pij=P0;
end
end
function [incleft,incright,inctop,incbottom]=incregion(data)
dis=floor(data.Nx/2)-3; %region=2*(dis-1)^2
incleft=floor(data.Nx/2)-dis;
incright=floor(data.Nx/2)+dis;
inctop=floor(data.Ny/2)+dis;
incbottom=floor(data.Ny/2)-dis;
end