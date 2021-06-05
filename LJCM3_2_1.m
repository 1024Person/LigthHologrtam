% --------------------------------------------------------------
%   
%   功能：      菲涅尔衍射的S-FFT计算
%   执行：调入一个图像，计算平面波照射下给定波长以及距离的衍射场振幅图像
%   主要变量：
%   h ----- 波长
%   z0 ---- 衍射距离
%   U0 ---- 初始光波场的复振幅
%   L0 ---- 初始场宽度
%   Uf ---- 衍射光波场的复振幅
%   L  ---- 衍射场宽度
% --------------------------------------------------------------

clear all;clc;close all;

mm = 1;
um = 1e-3;
nm = 1e-6;

rootPath = '~/';
[nom,chemin] = uigetfile([rootPath,'*.*'],['输入初始图像'],100,100);
[XRGB,MAP] = imread([chemin,nom]);
X = rgb2gray(XRGB); 
z0 = 1000*mm;
z0 = input('衍射距离z0=？(mm)');
h = 532*nm;
k = 2*pi/h;
[M,N] = size(X);
N = min([M,N]);
U0 = double(X);
L0 = sqrt(h*z0*N);
L0 = input('初始场宽度=？(mm)');
Uh = [0:N-1] - N/2;Vh = [0:N-1] - N/2;
[mh, nh] = meshgrid(Uh, Vh);
figstr = strcat('初始图像宽度=',num2str(L0),'mm');
figure(1),imshow(X,[]),colormap(gray);set(get(gca, 'YLabel'), 'String', '衍射计算');
xlabel(figstr);title('S-FFT方法计算衍射');
% -------------------------------- 菲涅耳衍射的S-FFT计算起始
n = 1:N;
x = -L0/2+L0/N*(n-1);
y = x;
[yy, xx] = meshgrid(y,x);
Fresnel = exp(1j*k/2/z0*(xx.^2+yy.^2));
f2 = U0.*Fresnel;
Uf = fft2(f2,N,N);
Uf = fftshift(Uf);
L = h*z0*N/L0;
x = -L/2 + L/N*(n-1);
y = x;
[yy, xx] = meshgrid(y,x);
phase = exp(1j*k*z0)/(1j*h*z0)*exp(1j*k/2/z0*(xx.^2+yy.^2));    % 计算菲涅尔积分前方的相位因子
Uf = Uf.*phase;
T = L0/N;
Uf = Uf*T*T;
% ------------------------------ 菲涅尔衍射的S-FFT计算结束
If = Uf.*conj(Uf);
figstr = strcat('衍射场图像宽度=',num2str(L),'mm');
figure(2),imshow(abs(Uf),[]),colormap(gray);set(get(gca, 'YLabel'), 'String', '衍射计算');
xlabel(figstr);title('S-FFT计算衍射');
