function [LCM_vector] = Seed_proc(dataf,Seedpos)
%SEED_PROC 循环处理所有种子点
%模板
template = [pi/4,  pi/2,  3*pi/4;
            0   ,   0  ,    pi  ;
           -pi/4, -pi/2,-3*pi/4];
% t_sim = [1,1,1;
%          1,0,1;
%          1,1,1];

list = 1 : 3 : 16;
% [ylist,xlist] = meshgrid(list);

vector_matrix = zeros(25,9);
vector_drt = zeros(25,9);
Ecorr = zeros(25,1);
Elist = [1,2,3,4,5,6,10,11,15,16,20,21,22,23,24,25];

Bcorr = zeros(25,1);
Blist = [7,8,9,12,14,17,18,19];

t_num = size(Seedpos,1);
LCM_vector = zeros(t_num,1);
LCM8_vector = zeros(t_num,1);
LCM9_vector = zeros(t_num,1);
% 填充dataf 避免目标出现在边界
sz = 7;
% datap = padarray(dataf,[sz sz],'both','replicate');
datap = padarray(dataf,[sz sz],1);
%


xsec = -sz:sz;
ysec = -sz:sz;
for t = 1 : t_num
    posx = Seedpos(t,2);
    posy = Seedpos(t,1);
    % 在图像中取出以种子点为中心的15*15区域
        Larea = datap(posx+7+xsec,posy+ysec+7); 
    % 计算梯度图像
        [Gx, Gy] = gradient(Larea);
    %     Gmag = sqrt(Gx.^2 + Gy.^2); %梯度图
        Gdir = atan2(Gy, Gx); % 方向图   
        % 向量化
        for x = 1 : 5
            for y = 1 : 5
                patch = Larea(list(x):list(x+1)-1,list(y):list(y+1)-1);
                vector_matrix((x-1)*5+y,:) = patch(:); %图像向量化
                patch = Gdir(list(x):list(x+1)-1,list(y):list(y+1)-1);
                vector_drt((x-1)*5+y,:) = patch(:); %方向向量化
            end
        end

    % k_means聚类    
    % 重塑矩阵为一个列向量
    data = reshape(Larea, [], 1);
    tic
    % 运行K-means聚类，分为3类
    k = 3;  % 聚类的目标数
    [idx, C] = kmeans(data, k);
    toc
    % 重塑聚类结果到原始矩阵的大小
    clusteredMatrix = reshape(idx, 15, 15);
    
    close all 
    figure
    % 显示聚类结果
    imagesc(clusteredMatrix);  % 使用颜色编码显示聚类结果
    colormap('jet');  % 设置颜色图
    colorbar;  % 显示颜色条
    title('Clustered Image');

    figure
    imshow(Larea, [])
%     % 基于目标梯度方向各向异性估计目标能量
%         factorlist = cos( vector_drt- template(:)');
%         factor = sum(factorlist,2)/9; %方向加权因子
% 
%         patch_mean = mean(vector_matrix,2); %块均值
% 
%         patch_factor = patch_mean.*factor; %块均值加权
%             
%         Ecorr(Elist) = patch_factor(Elist);
%         Emax = find(Ecorr == max(Ecorr));
%         E_v = Ecorr(Emax);
% 
% 
% 
%         Bcorr(Blist) = patch_factor(Blist);
%         Bmax = find(Bcorr == max(Bcorr));
%         B_v = Bcorr(Bmax);
%         
%         T_v = patch_factor(13);
%         % 计算能量
%         LCM9_vector(t) = T_v.^2/(E_v*B_v);
    % 基于相似度度量估计边缘区域和背景区域能量z

%         vector_t = vector_matrix(13,:)'; %目标向量
%         % 计算皮尔森相关系数
%         % 13 为自相似性
%         % [1,2,3,4,5,6,10,11,15,16,20,21,22,23,24,25]
%         % [7,8,9,12,14,17,18,19]
%         pcorr = Cossim(vector_t,vector_matrix);
% 
%         Ecorr(Elist) = pcorr(Elist);
%         Emax = find(abs(Ecorr) == max(abs(Ecorr)));
% 
%         Bcorr(Blist) = pcorr(Blist);
%         Bmax = find(abs(Bcorr) == max(abs(Bcorr)));
%         
%         E_v = mean(vector_matrix(Emax,:));
%         B_v = mean(vector_matrix(Bmax,:));
%         
% 
%         LCMp_vector(t) = T_v.^2/(E_v*B_v);

%         %%
%         pcorr = mean(vector_matrix,2);
% 
%         Ecorr(Elist) = pcorr(Elist);
%         Emax = find(abs(Ecorr) == max(abs(Ecorr)));
% 
%         Bcorr(Blist) = pcorr(Blist);
%         Bmax = find(abs(Bcorr) == max(abs(Bcorr)));
%         
%         E_v = mean(vector_matrix(Emax(1),:));
%         B_v = mean(vector_matrix(Bmax(1),:));
%     
%         LCM8_vector(t) = T_v.^2/(E_v*B_v);
%         %%
%         pcorr = mean(vector_matrix,2);
% 
%         Ecorr(Elist) = pcorr(Elist);
%         Emax = find(abs(Ecorr) == max(abs(Ecorr)));
% 
%         Bcorr(Blist) = pcorr(Blist);
%         Bmax = find(abs(Bcorr) == max(abs(Bcorr)));
%         
%         E_v = mean(vector_matrix(Emax(1),:));
%         B_v = mean(vector_matrix(Bmax(1),:));
%         T_v = mean(vector_matrix(13,:));
%         LCM9_vector(t) = T_v.^2/(E_v*B_v);
end
end

