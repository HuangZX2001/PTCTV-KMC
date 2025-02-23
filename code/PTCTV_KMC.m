function [res] = PTCTV_KMC(img)
%PTCTV-KMC

%% 预处理
I = double(img);

if size(I,3) == 3
    dataf = rgb2gray(I / 255);
else
    dataf = I / 255;
end   

%% 检测
%—————张量分离—————
[tarImg] = TRPCA(dataf,0.85);
%——————密度峰聚类———————
[Seedpos,~] = SearchPeak(dataf,tarImg);

    if size(Seedpos,1) < 1
%          disp(['0'])
         %—————张量主成分分析—————
          [tarImg] = TRPCA(dataf,0.85*0.6);
         %——————密度峰聚类———————
         [Seedpos,~] = SearchPeak(dataf,tarImg);
    end

%———————聚类分析———————
[~,Weight] = seed_kmeans(dataf,Seedpos);
Weight(Weight<=1) = 0;
%—————当前帧检测结果———
resf = tarImg.*Weight;
resf(resf<0) = 0;
res = resf; 
end

