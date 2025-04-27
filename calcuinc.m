clear all
clc
delete(gcp('nocreate'));
parpool('local', 10);
dt=0.2;
endt=0.2;
path='E:\jet\jet_290x260\';


fine_positions=readmesh(path) ; 
data=coarse(); 
parfor tstep=1:1:100
 %time=0;dt;endt
time=tstep*0.005;

velData=readvel(path,time) ;
num_fine = size(velData, 1);
pressureData=readpress(path,time,num_fine) ;
radius=0.03;
[coarse_data_p, coarse_data_U, coarse_data_V]= radial_basis_interpolation(fine_positions, pressureData,velData, radius,data);

surf(data.xx,data.yy,coarse_data_p)
shading interp
view(0,90)
xlabel('X')
ylabel('Y')

filename = [path num2str(time) '\coarse_Pij.dat'];
output(coarse_data_p, filename)

filename = [path num2str(time) '\coarse_Uij.dat'];
output(coarse_data_U, filename)

filename = [path num2str(time) '\coarse_Vij.dat'];
output(coarse_data_V, filename)
disp([num2str(time) 'is completed']);
end
delete(gcp);

function data=coarse()
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

end
function fine_position=readmesh(path) 
    fileID = fopen([path 'C.txt'], 'r');
    headerLines = 0; % 跳过文件开头的行数
    for i = 1:headerLines
    fgetl(fileID);
    end
    
    data = textscan(fileID, '(%f %f %f)', 'CollectOutput', true);
fclose(fileID);

fine_position = data{1};
    
end
function pressureData=readpress(path,time,num_fine)
% 打开文件
fileID = fopen([path num2str(time) '\p'], 'r');
% 跳过文件头部信息
while ~feof(fileID)
    line = fgetl(fileID);
    if contains(line, 'internalField')
        break;
    end
end
% 读取压力数据
pressureData = zeros(num_fine,1);num=1;
while ~feof(fileID)
    line = fgetl(fileID);
    if strcmp(line, '(')
        % 读取压力数据直到遇到')'为止
        while ~feof(fileID)
            line = fgetl(fileID);
            if strcmp(line, ')')
                break;
            else
                pressureData(num,1) = str2double(line);
                num=num+1;
            end
        end
        break;
    end
end
% 关闭文件
fclose(fileID);
% 压力数据现在存储在pressureData变量中，可以进行后续处理或分析
end
function velocityData=readvel(path,time)
fileID = fopen([path num2str(time) '\U'], 'r');
headerLines = 23; % 跳过文件开头的行数
for i = 1:headerLines
    fgetl(fileID);
end

% 读取速度数据
data = textscan(fileID, '(%f %f %f)', 'CollectOutput', true);
fclose(fileID);

% 获取速度数据
velocityData = data{1};
% 压力数据现在存储在pressureData变量中，可以进行后续处理或分析
end
function [coarse_data_p, coarse_data_U, coarse_data_V]= radial_basis_interpolation(fine_positions, fine_p,fine_U, radius,data)

    num_fine = size(fine_positions, 1);
    coarse_data_p = zeros(data.Ny, data.Nx);
    coarse_data_U = zeros(data.Ny, data.Nx);
    coarse_data_V = zeros(data.Ny, data.Nx);

    % 遍历每个粗网格节点
    for i = 1:data.Ny
        for  j=1:data.Nx
        % 当前粗网格节点的位置
        coarse_pos = [data.x(j) data.y(i)];

        % 遍历每个细网格节点
        weighted_sum_p = 0;
        total_weight_p = 0;
        weighted_sum_U = 0;
        total_weight_U = 0;
        weighted_sum_V = 0;
        total_weight_V = 0;
        for nodenum = 1:num_fine
            % 当前细网格节点的位置和物理量 x & z
            fine_pos = [fine_positions(nodenum, 1) fine_positions(nodenum, 3)];
            fine_value_p = fine_p(nodenum);
            fine_value_U = fine_U(nodenum,1);
            fine_value_V = fine_U(nodenum,3);
            % 计算细网格节点到粗网格节点的距离
            distance = norm(fine_pos - coarse_pos);

            % 高斯径向基函数
            weight = exp(-(distance / radius)^2);

            % 更新加权和和总权重
            weighted_sum_p = weighted_sum_p + weight * fine_value_p;
            total_weight_p = total_weight_p + weight;
            weighted_sum_U = weighted_sum_U + weight * fine_value_U;
            total_weight_U = total_weight_U + weight;
            weighted_sum_V = weighted_sum_V + weight * fine_value_V;
            total_weight_V = total_weight_V + weight;
        end
        

        % 对粗网格节点上的物理量进行插值
            coarse_data_p(i,j) = weighted_sum_p / total_weight_p;
            coarse_data_U(i,j) = weighted_sum_U / total_weight_U;
            coarse_data_V(i,j) = weighted_sum_V / total_weight_V;
        
        end
    end
end
function output(matrix, filename)
     % 打开 ".dat" 文件以写入数据
    fid = fopen(filename, 'w');

    % 获取矩阵的大小
    [rows, cols] = size(matrix);

    % 将矩阵数据按照指定格式写入文件
    for row = 1:rows
        % 将每一行的数据转换为带有指定精度和间隔的字符串
        rowStr = sprintf('%.15e ', matrix(row, 1:cols));
        rowStr = [rowStr, sprintf('\n')]; % 每行数据后添加换行符

        % 将带有指定格式的字符串写入文件
        fprintf(fid, '%s', rowStr);
    end

    % 关闭文件
    fclose(fid);
end
