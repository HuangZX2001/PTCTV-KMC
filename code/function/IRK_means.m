function [clusteredMatrix] = IRK_means(Larea,k)
%IRK_MEANS 此处显示有关此函数的摘要

    % 重塑矩阵为一个列向量
    data = reshape(Larea, [], 1);
  % 运行K-means聚类，分为3类
%     k = k;  % 聚类的目标数
    [idx, C] = kmeans(data, k,'Start','plus');
  % 重塑聚类结果到原始矩阵的大小
    clusteredMatrix = reshape(idx, 15, 15);
    
end

