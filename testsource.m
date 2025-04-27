clear all
clc
data.Nx=100;data.Ny=100;%100列，100行
data.dt=0.5;data.c=1.0;data.r0=1.0;data.gamma=0.2*pi;
a=20*data.r0;

data.top=a;data.bottom=-a;data.left=-a;data.right=a;
data.x=linspace(data.left,data.right,data.Nx);
data.y=linspace(data.bottom,data.top,data.Ny);
data.var=4;
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

%a=-0.5;b=0.5;
 a=-3e-5*bei;b=3e-5*bei ;
 Source=zeros(data.Ny,data.Nx,4);
for tstep=1:1:data.Nt
  for i=2:1:data.Ny-1
      for j=2:1:data.Nx-1
          Source(i,j,:)=buildS(i,j,data,tstep);
      end
  end
    
    
    
    figure(1)
    surf(data.xx,data.yy,Source(:,:,4))
    
    
    
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

function [tao_1,tao_2]=viscoustensor(i,j,dqdx,dqdy,data)
miu=1*10e-3;

       [dudxdx,dudxdy]=center4order(i,j,dqdx(:,:,2),data);
       [dvdxdx,dvdxdy]=center4order(i,j,dqdx(:,:,3),data);
       [dudydx,dudydy]=center4order(i,j,dqdy(:,:,2),data);
       [dvdydx,dvdydy]=center4order(i,j,dqdy(:,:,3),data);
       tao_1=miu*(4/3*dudxdx+dudydy+1/3*(0.5*dvdxdy+0.5*dvdydx));
       tao_2=miu*(4/3*dvdydy+dvdxdx+1/3*(0.5*dudxdy+0.5*dudydx));
       

end
function [differfX,differfY]=center4order(i,j,f,data)
differfX=(8*(f(i,j+1,1)-f(i,j-1,1))-(f(i,j+2,1)-f(i,j-2,1)))/(12*data.dx);
differfY=(8*(f(i+1,j,1)-f(i-1,j,1))-(f(i+2,j,1)-f(i-2,j,1)))/(12*data.dy);
end
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
function [dqdx,dqdy]=center4orderall(q,data,dqdx,dqdy)

%%%dqdx%%%
for i=1:1:data.Ny
    for j=1:1:data.Nx
        if j>3&&j<data.Nx-2
           
            dqdx(i,j,:)=(8*(q(i,j+1,:)-q(i,j-1,:))-(q(i,j+2,:)-q(i,j-2,:)))/(12*data.dx);
        
        end
    end 
end
%%%dqdy%%%
for j=1:1:data.Nx
    for i=1:1:data.Ny
        if i>3&&i<data.Nx-2
            
            dqdy(i,j,:)=(8*(q(i+1,j,:)-q(i-1,j,:))-(q(i+2,j,:)-q(i-2,j,:)))/(12*data.dy);
        
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
function [dqdx,dqdy]=compact10order(dqdx,dqdy,q,data)
%%%6order%%%
lambda_6=0.4111403764203249;
a_6=0.7842616980320271;
b_6=0.0692748674241733;
c_6=-0.0038903521543495;
%%%10order%%%
lambda_10=0.4388871532438393;
a_10=0.7748150462341548;
b_10=0.0962949739000680;
c_10=-0.0110116258189503;
d_10=0.0012257054792086;
e_10=-0.0000771570500869;


A_x=zeros(data.Nx-8,data.Nx-8);
d_x=zeros(data.Nx-8,data.var);
A_y=zeros(data.Ny-8,data.Ny-8);
d_y=zeros(data.Ny-8,data.var);
for j=1:1:data.Nx-8
    if j==1
        A_x(1,1)=1;
        A_x(1,2)=lambda_6;
    elseif j==data.Nx-8
        A_x(data.Nx-8,data.Nx-8)=1;
        A_x(data.Nx-8,data.Nx-9)=lambda_6;    
    elseif j==2||j==3||j==data.Nx-7||j==data.Nx-6
        A_x(j,j)=1;
        A_x(j,j+1)=lambda_6;
        A_x(j,j-1)=lambda_6;
    else
        A_x(j,j)=1;
        A_x(j,j+1)=lambda_10;
        A_x(j,j-1)=lambda_10;        
    end
end
for j=1:1:data.Ny-8
    if j==1
        A_y(1,1)=1;
        A_y(1,2)=lambda_6;
    elseif j==data.Ny-8
        A_y(data.Ny-8,data.Ny-8)=1;
        A_y(data.Ny-8,data.Ny-9)=lambda_6;    
    elseif j==2||j==3||j==data.Ny-7||j==data.Ny-6
        A_y(j,j)=1;
        A_y(j,j+1)=lambda_6;
        A_y(j,j-1)=lambda_6;
    else
        A_y(j,j)=1;
        A_y(j,j+1)=lambda_10;
        A_y(j,j-1)=lambda_10;  
    end
end

%%%dqdx%%%
for i=1:1:data.Ny
     for j=1:1:data.Nx-8
          if j==1
                d_x(1,:)=a_6*(q(i,6,:)-q(i,4,:))/data.dx+b_6*(q(i,7,:)-q(i,3,:))/data.dx+c_6*(q(i,8,:)-q(i,2,:))/data.dx-lambda_6*dqdx(i,4,:);
          elseif j==data.Nx-8
                d_x(data.Nx-8,:)=a_6*(q(i,data.Nx-3,:)-q(i,data.Nx-5,:))/data.dx+b_6*(q(i,data.Nx-2,:)-q(i,data.Nx-6,:))/data.dx+c_6*(q(i,data.Nx-1,:)-q(i,data.Nx-7,:))/data.dx-lambda_6*dqdx(i,data.Nx-3,:);          
          elseif j==2||j==3||j==data.Nx-7||j==data.Nx-6
                d_x(j,:)=a_6*(q(i,j+5,:)-q(i,j+3,:))/data.dx+b_6*(q(i,j+6,:)-q(i,j+2,:))/data.dx+c_6*(q(i,j+7,:)-q(i,j+1,:))/data.dx;
          else
              d_x(j,:)=a_10*(q(i,j+5,:)-q(i,j+3,:))/data.dx+b_10*(q(i,j+6,:)-q(i,j+2,:))/data.dx+c_10*(q(i,j+7,:)-q(i,j+1,:))/data.dx+d_10*(q(i,j+8,:)-q(i,j,:))/data.dx+e_10*(q(i,j+9,:)-q(i,j-1,:))/data.dx;
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
             d_y(1,:)=a_6*(q(6,j,:)-q(4,j,:))/data.dy+b_6*(q(7,j,:)-q(3,j,:))/data.dy+c_6*(q(8,j,:)-q(2,j,:))/data.dy-lambda_6*dqdy(4,j,:);
          elseif i==data.Ny-8
              d_y(data.Ny-8,:)=a_6*(q(i+5,j,:)-q(i+3,j,:))/data.dy+b_6*(q(i+6,j,:)-q(i+2,j,:))/data.dy+c_6*(q(i+7,j,:)-q(i+1,j,:))/data.dy-lambda_6*dqdy(data.Ny-3,j,:);             
          elseif i==2||i==3||i==data.Nx-7||i==data.Nx-6              
              d_y(i,:)=a_6*(q(i+5,j,:)-q(i+3,j,:))/data.dy+b_6*(q(i+6,j,:)-q(i+2,j,:))/data.dy+c_6*(q(i+7,j,:)-q(i+1,j,:))/data.dy;
          else
              d_y(i,:)=a_10*(q(i+5,j,:)-q(i+3,j,:))/data.dy+b_10*(q(i+6,j,:)-q(i+2,j,:))/data.dy+c_10*(q(i+7,j,:)-q(i+1,j,:))/data.dy+d_10*(q(i+8,j,:)-q(i,j,:))/data.dy+e_10*(q(i+9,j,:)-q(i-1,j,:))/data.dy;
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

function [S]=buildS(i,j,data,tstep)
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
    %S4=-dPdt-ratio*(q(i,j,2)+Uij*q(i,j,1)/rho)*dPdx-ratio*(q(i,j,3)+Vij*q(i,j,1)/rho)*dPdy;
    S4=-dPdt-Uij*dPdx-Vij*dPdy;
%     S2=-q(i,j,1)/rho*(dUdt+Uij*dUdx+Vij*dUdy)-q(i,j,2)*dUdx-q(i,j,3)*dUdy;
%     S3=-q(i,j,1)/rho*(dVdt+Uij*dVdx+Vij*dVdy)-q(i,j,2)*dVdx-q(i,j,3)*dVdy;
%     S4=-dPdt-ratio*(q(i,j,2)+q(i,j,1)/rho*Uij)*dPdx-ratio*(q(i,j,3)+q(i,j,1)/rho*Vij)*dPdy;
    
%     S4=-dPdt-u*dPdx-v*dPdy;
    % S=[0;0;0;S4;];
    %S4=0.0;
   S=[0;0;0;S4;];
    %S=[0;0;0;0];

end


function [Uij,Vij,Pij]=buildinc(data,i,j,tstep)
gamma=data.gamma;
r0=data.r0;
w=gamma/(4*pi*r0^2);rho0=1.0;P0=rho0*data.c^2/1;
%if (i>incbottom&&i<inctop)&&(j>incleft&&j<incright)
if (data.y(i)^2+data.x(j)^2<=(300*data.r0)^2)
%     
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
%   %   Uij=0.5*data.c;
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
