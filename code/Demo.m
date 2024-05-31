clear
clc
dbstop if error
%——————Demo———————

img = imread("0.bmp");

res = PTCTV_KMC(img);

imshow(res,[])



