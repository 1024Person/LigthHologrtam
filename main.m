clc;clear all;close all;
disp('开始计算正衍射')
Uz = ClassicalDiffraction();
figure
imshow(abs(Uz),[]),colormap(gray);
title('衍射图样')

disp('开始计算逆衍射')
U0 = InverseClassicalDiffraction(Uz);
figure
imshow(abs(U0),[]),colormap(gray);
title('复现图样')