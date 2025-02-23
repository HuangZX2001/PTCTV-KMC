clear
clc
dbstop if error
addpath('function\')
addpath('lib\')
%——————Demo———————

img = imread("0.bmp");

res = PTCTV_KMC(img);

imshow(res,[])



