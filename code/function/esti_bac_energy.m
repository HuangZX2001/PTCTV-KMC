function [B_energy] = esti_bac_energy(Larea,mask,num)
%ESTI_BAC_ENERGY 估计背景区域能量  大于9个像素时，取前9个最大的


B = Larea.*mask;

tnum = ceil(num/2);

if num <= 9
    % 直接计算均值
    B_energy = sum(B,'all')/num;
elseif num > 9
    B = B(:);
    B = sort(B,'descend');
    B_energy = sum(B(1:9))/9;
end


% if num <= tnum
%     % 直接计算均值
%     B_energy = sum(B,'all')/num;
% elseif num > tnum
%     B = B(:);
%     B = sort(B,'descend');
%     B_energy = sum(B(1:tnum))/tnum;
% end