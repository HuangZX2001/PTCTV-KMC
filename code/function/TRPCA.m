function [tarImg] = TRPCA(dataf,rho)
% 执行张量分离
opts.lambda = rho;
 % 导入超参数  
patch = 50;
step = 40;

% 分解图像
tenImg = double(gen_patch_ten(dataf, patch, step));
[~, E, ~, ~] = PTCTV_KMC_Model(tenImg,opts);

tarImg = res_patch_ten_mean(E, dataf, patch, step);
tarImg(tarImg<0) = 0;   

end

