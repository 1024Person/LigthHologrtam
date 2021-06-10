clc;clear all;close all;
disp('开始计算正衍射')
Uz = ClassicalDiffraction();
figure
imshow(abs(Uz),[]),colormap(gray);
title('衍射图样')
temp = abs(Uz)/max(max(abs(Uz))) * 255;
imwrite(temp,gray,'./result/temp.bmp');
% disp('开始计算逆衍射')
% U0 = InverseClassicalDiffraction(Uz);
Uzz = ClassicalDiffraction();


figure
imshow(abs(Uzz),[]),colormap(gray);
title('复现图样')