function [LCM,Weight] = seed_kmeans(dataf,Seedpos)
%SEED_KMEANS 
[dim1,dim2] = size(dataf);
% Weight = ones(size(dataf));

t_num = size(Seedpos,1); %目标数量


sz = 7;
datap = padarray(dataf,[sz sz],1); % 填充dataf 避免目标出现在边界
sw = 1;

xsec = -sz:sz;
ysec = -sz:sz;
wsec = -sw:sw;

mask = zeros(15,15);
mask(7:9,7:9) = 1;
maske = zeros(15,15);
maske(4:12,4:12) = 1; 

W5 = zeros(15,15);
W5(6:10,6:10) = 1;


mask_set = find(mask == 1);
LCM = zeros(t_num,1);
Wt = ones(dim1,dim2,t_num);
Weightp = ones(dim1,dim2);
Weightp = padarray(Weightp,[sw sw],1);
Wt = padarray(Wt,[sw sw],1); % 填充dataf 避免目标出现在边界
%循环处理每一个种子点
for t = 1 : t_num
    posx = Seedpos(t,2);
    posy = Seedpos(t,1);
    % 在图像中取出以种子点为中心的15*15区域
        Larea = datap(posx+7+xsec,posy+ysec+7); 
    % k_means聚类
        [cluster] = IRK_means(Larea,3);
    % 显示聚类结果
%         close all 
%         figure
%         imagesc(cluster);  % 使用颜色编码显示聚类结果
%         colormap('jet');  % 设置颜色图
%         colorbar;  % 显示颜色条
%         title('Clustered Image');
%         figure
%         imshow(Larea, [])
   % 读取种子点标签
        
    % 分配mode
       % 场景0 种子点位置的连通域只包含其自身像素，考虑其为单像素噪声
       % 场景1 种子点位置的连通域处于3*3窗口内
       % 场景2 种子点位置的连通域超出3*3窗口
        seed_c = cluster(8,8);
        cluster_seed = zeros(15,15);
        cluster_seed(cluster == seed_c) = 1; %提取所有与种子点同类别的像素

        L = bwlabel(cluster_seed,4);%连通域分析 4邻域 
        seed_L = L(8,8);
        L_seed = find(L == seed_L); %提取目标区域

        seed_L_num = numel(L_seed);
        L_mask = zeros(15,15);
        L_mask(L == seed_L) = 1;
         
        if seed_L_num > 81
            mode = -1;
        end
       
        if seed_L_num == 1     

            test = Larea(7:9,7:9);
            test(5) = [];
            test(find(test == min(test))) = [];
            if std(test) <= 0.025 && sum(cluster_seed.*W5,'all') == 1
                mode = -1;
            else
                mode = 1; 
            end
        elseif seed_L_num > 1
            % 如果连通域完全在窗口中心，subsetCheck = 0，反之为1
            subsetCheck = ~all(ismember(L_seed, mask_set));
            % 将subsetCheck = 0 和 mode 对应
            mode = subsetCheck + 1;
        end
        
   % 计算LCM

   mask_center = zeros(15,15);
   mask_B1 = zeros(15,15);
   mask_B2 = zeros(15,15);
   mask_B3 = zeros(15,15);
        switch mode
            case -1
                LCM(t) = 0;
            case 0
                % 计算目标能量
                mask_center(L == seed_L) = 1; 
                T = mask_center.*Larea;
                T_num = sum(T>0,'all'); 
                T_energy = esti_maxnum_energy(Larea,mask);

                % 统计三种标签的背景能量 取最大值作为最终背景
                mask_B1(cluster == 1) = 1; mask_B1 = mask_B1 .*(~mask); num_B1 = sum(mask_B1,'all');
                mask_B2(cluster == 2) = 1; mask_B2 = mask_B2 .*(~mask); num_B2 = sum(mask_B2,'all');
                mask_B3(cluster == 3) = 1; mask_B3 = mask_B3 .*(~mask); num_B3 = sum(mask_B3,'all');

                eval(['mask_B',num2str(seed_c),'=','mask_B',num2str(seed_c),'.*','~mask_center','.*','~mask',';'])
                eval(['num_B',num2str(seed_c),'=','sum(','mask_B',num2str(seed_c),',','''all''',')',';'])
               
                B1_energy = esti_bac_energy(Larea,mask_B1,num_B1);
                B2_energy = esti_bac_energy(Larea,mask_B2,num_B2);
                B3_energy = esti_bac_energy(Larea,mask_B3,num_B3);
                B_energy_list = [B1_energy,B2_energy,B3_energy];
                B_energy = max(B_energy_list); 

                % 计算对比度
                LCM(t) = T_energy/B_energy;
            case 1
                % 计算目标能量
                mask_center(L == seed_L) = 1; 
                T = mask_center.*Larea;
                T_num = sum(T>0,'all'); 
                T_energy = max(T,[],'all'); %% 修改
                % 统计三种标签的背景能量 取最大值作为最终背景
                mask_B1(cluster == 1) = 1;  num_B1 = sum(mask_B1,'all');
                mask_B2(cluster == 2) = 1;  num_B2 = sum(mask_B2,'all');
                mask_B3(cluster == 3) = 1;  num_B3 = sum(mask_B3,'all');

                eval(['mask_B',num2str(seed_c),'=','mask_B',num2str(seed_c),'.*','~mask_center',';'])
                eval(['num_B',num2str(seed_c),'=','sum(','mask_B',num2str(seed_c),',','''all''',')',';'])
               
                B1_energy = esti_bac_energy(Larea,mask_B1,num_B1);
                B2_energy = esti_bac_energy(Larea,mask_B2,num_B2);
                B3_energy = esti_bac_energy(Larea,mask_B3,num_B3);

            %                 remove center
                eval(['CNUM','=','num_B',num2str(seed_c),';' ])
                if CNUM <= 3
                    eval(['B',num2str(seed_c),'_energy','=','0',';'])
 
                end

                B_energy_list = [B1_energy,B2_energy,B3_energy];
                B_energy = max(B_energy_list); 
    
                % 计算对比度
                LCM(t) = T_energy/B_energy;
            case 2
                % 计算目标能量
                mask_center(L == seed_L) = 1; 
                maskts = mask_center.*mask;
                T = maskts.*Larea;
                T_num = sum(T>0,'all');  
                %
                T_energy = max(T,[],'all');
                %
%                 T_energy = sum(T,'all')/T_num;
                
                % 统计三种标签的背景能量 取最大值作为最终背景
                mask_B1(cluster == 1) = 1; mask_B1 = mask_B1 .*(~mask); num_B1 = sum(mask_B1,'all');
                mask_B2(cluster == 2) = 1; mask_B2 = mask_B2 .*(~mask); num_B2 = sum(mask_B2,'all');
                mask_B3(cluster == 3) = 1; mask_B3 = mask_B3 .*(~mask); num_B3 = sum(mask_B3,'all');

                % 统计边缘能量      maske
                eval(['mask_E','=','mask_B',num2str(seed_c),'.*','L_mask','.*','maske',';']);num_E = sum(mask_E,'all');
                eval(['mask_B',num2str(seed_c),'=','mask_B',num2str(seed_c),'.*','~mask_E','.*','~maskts',';'])
                eval(['num_B',num2str(seed_c),'=','sum(','mask_B',num2str(seed_c),',','''all''',')',';'])
                

                E_energy = esti_bac_energy(Larea,mask_E,num_E);
                
                if num_B1 == 0 
                    B1_energy = 0;
                else
                    B1_energy = esti_bac_energy(Larea,mask_B1,num_B1);
                end

                if num_B2 == 0 
                    B2_energy = 0;
                else
                    B2_energy = esti_bac_energy(Larea,mask_B2,num_B2);
                end


                if num_B3 == 0 
                    B3_energy = 0;
                else
                    B3_energy = esti_bac_energy(Larea,mask_B3,num_B3);
                end

                eval(['CNUM','=','num_B',num2str(seed_c),';' ])
                if CNUM <= 3
                    eval(['B',num2str(seed_c),'_energy','=','0',';'])
                end
                

                B_energy_list = [B1_energy,B2_energy,B3_energy];
                B_energy = max(B_energy_list);
    
                
                %
                t1 = T_energy/E_energy;
                t1(t1>1) = 1;
                LCM(t) = (T_energy/B_energy).*t1;
        end

    % 构建权重矩阵
    wsec = 0;
    Wt(posx+sw+wsec,posy+wsec+sw,t) = LCM(t);%当前目标的权重矩阵
    
end
    for t = 1 : t_num
        Weightp = Weightp.*Wt(:,:,t);
    end
    % 裁剪权重矩阵
    [dim1,dim2] = size(Weightp);
    Weight = Weightp(2:dim1-1,2:dim2-1);
end


