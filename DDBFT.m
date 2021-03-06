%% DDBFT算法重建
close all;clc;clear all;
%% 准备工作

[fileName,filePath] = uigetfile(['*.*'],'输入全息图',100,100);

I1 = imread([filePath,fileName]);

figure(1),imshow(I1),title('数字全息图');
%% 初始化数据工作

f0 = double(I1);

[N1,N2] = size(f0);
N = min(N1,N2);
h = input('波长：？(mm)');
L = input('CCD尺寸：?(mm)');
z0 = input('衍射距离：？(mm)');
k = 2*pi/h;     % 波矢
L0 = h*z0*N/L;

Ih = zeros(N,N);
Ih(1:N,1:N) = f0(1:N,1:N);
%% 1-FFT重建开始
n = 1:N;
x = -L/2 + L/N * (n-1);
y = x;

[xx,yy] = meshgrid(x,y);
Fresnel = exp(1j*k/2/z0 *(xx.^2+yy.^2));
f2 = Ih.*Fresnel;
Uf = fft2(f2,N,N);
Uf = fftshift(Uf);
x = -L0 / 2 + L0/N *(n-1);
y = x;
[xx,yy] = meshgrid(x,y);
phase = exp(1j*k*z0) / (1j*h*z0) * exp(1j*k/2/z0*(xx.^2+yy.^2));
Uf = Uf.*phase;
% 1-FFT重建结束
%% 显示
figstr0 = strcat('物重建平面宽度=',num2str(L0),'mm');
figure(2),
Gmax = max(max(abs(Uf)));
Gmin = min(min(abs(Uf)));
imshow(abs(Uf)+eps ,[Gmin Gmax / 10]),xlabel(figstr0);
title('1-FFT物平面重建图像');
%% 开始“接力”？？？

[xc,yc] = ginput(1);        % 鼠标左肩选择重建局部图像中心
Lx = L0 / N*(N-xc);
Ly = L0 /N*(N-yc);
r = sqrt(Lx*Lx+z0*z0);
Qx = Lx/r;
Qy = Ly / r;

Ufi = zeros(N,N);
dx = input('X向半宽度（像素）=？');
dy = input('Y向半宽度（像素）=？');
dxy = max(dx,dy);
figstr = strcat('重建物平面宽度=',num2str(L0/N*dxy*2),'mm');
if xc - dx < 1
    xb0 = 1;
else 
    xb0 = xc - dx;
end
if yc - dy<1
    yb0 = 1;
else
    yb0 = yc-dy;
end



