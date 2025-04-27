% syms u v rho ratio p M;  % 声明符号变量
% 
% A = [v 0 rho 0;0 v 0 0;0 0 v 1/(ratio*M^2*rho);0 0 ratio*p v];  % 符号矩阵
% 
% eigenvalues = eig(A);  % 计算特征值
% 
% disp(eigenvalues);  % 显示特征值
clear all
 data.Nx=60;
    data.Ny=60;
    data.incregion_x=2.4;
    data.incregion_y=2;
    data.dt=0.005;
    data.Nt=100;
    data.x=linspace(0,data.incregion_x,data.Nx);
    data.y=linspace(0,data.incregion_y,data.Ny);
    for k=1:1:data.Nx
    data.yy(:,k)=data.y;
    end
    for k=1:1:data.Ny
    data.xx(k,:)=data.x;
    end
for time=0.005:0.005:0.5
    path=['E:\jet\jet_290x260\' num2str(time) '\coarse_Pij.dat'];
    P=importdata(path);
    surf(data.xx,data.yy,P)
shading interp
view(0,90)
xlabel('X')
ylabel('Y')
    set(gca,'CLim',[-0.2 0.2])
    title(['time=',num2str(time)])
pause(0.25)



end