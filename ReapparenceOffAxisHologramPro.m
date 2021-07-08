% 数字全息重建
%% 准备工作
clear all; clc;close all;
[fileName,filePath] = uigetfile('*.*','输入全息图','./result/');
I1 = imread([filePath,fileName]);
figure(1),imshow(I1),colormap(gray);title('全息图');

f = double(I1);
[N1,N2] = size(f);
N = max([N1,N2]);
Ih = zeros(N,N);
Ih(N/2 - N1/2+1:N/2+N1/2,N/2 - N2/2 + 1:N/2+N2/2) = f(1:N1,1:N2);

%% 输入数据

h = input('波长=?(mm)');
L = input('CCD尺寸=?(mm)');
z0= input('重建距离=?(mm)');
k = 2*pi / h;
L0 = h*k*z0/L;  

%% 1-FFT重建开始
% 离散化坐标(物平面)
n = 1:N;
x = -L/2 + L/N * (n-1);
y = x;
% 网格化数据
[xx,yy] = meshgrid(x,y);

fresenl = exp(1j*k/2/z0*(xx.^2+yy.^2));
f2  = Ih.*fresenl;
% 菲涅耳衍射的积分部分
Uf = fft2(f2,N,N);
Uf = fftshift(Uf);

% 再现平面的真实离散坐标
x = -L0/2+L0/N*(n-1);
y = x;
[xx,yy] = meshgrid(x,y);


% 菲涅耳衍射1-FFT前面的相位部分
phase = exp(1j*k*z0) / (1j*h*z0) * exp(1j*k/2/z0*(xx.^2+yy.^2));
% 再现图的复振幅表达式
Uz = Uf.*phase;
figure(2),imshow(abs(Uz) / max(max(abs(Uz))));
%% 全息图再现
If = Uz .* conj(Uz);
Imax = max(max(abs(If)));
figure(3),imshow(If/Imax),colormap(gray);
title('再现图');

