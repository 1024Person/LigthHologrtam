%--------------- 消除零级衍射的离轴数字全息图
%           1-FFT重建
%-------------------------------------------
%% 准备工作
clc;close all;
clear all;

[fileName,filePath] = uigetfile('*.*','制作全息图','./res/');
[X0,cmap] = imread([filePath,fileName]);
if length(size(X0)) > 2
    X0 = rgb2gray(X0);
end
X0 = X0 / max(max(abs(X0)));

[r,c] = size(X0);
N1 = max([r,c]);
N = 1024;   % 抽样点数
X1 = imresize(X0,N/2/N1); % 将最大的边长限制在了N/2,菲涅尔衍射前后光场的尺寸不会发生变化。

[M1,N1] = size(X1);
X = zeros(N,N);
% 这里加1是因为，matlab中数组下标是从1开始的
X(N/2 - M1 /2 + 1:N/2 + M1/2,N/2 - N1/2 + 1:N/2+N1/2) = X1(1:M1,1:N1);
figure(1);
imshow(X),colormap(gray);
title('扩充之后的原图');
%% 输入数据
h = input('波长：？(mm)');
k = 2*pi / h;
L = input('CCD面阵宽度：?(mm)');
z0 = input('衍射距离：?(mm)');
L0 = h*k*z0/L;
%% 开始菲涅尔衍射 1-FFT
% 随机相位
phase = rand(N).*2*pi;
f = double(X).*exp(1j.*phase);

n = 1:N;
% 真实的物理坐标
% 物平面真实的坐标
x = -L0/2+L0/N*(n-1);
y = x;
% 网格化数据
[xx,yy] = meshgrid(x,y);
frensel = exp(1j*k/2/z0*(xx.^2+yy.^2));
f = fftshift(f);
f2 = f.*frensel;
% 菲涅耳衍射的积分部分
Uf = fft2(f2,N,N);

% CCD平面真实的坐标
x = L/2 - L/N * (n-1);
y = x;
[xx,yy] = meshgrid(x,y);
% 菲涅耳衍射前面的相位部分
phase = exp(1j*k*z0) / (1j*h*z0) * exp(1j*k/2/z0*(xx.^2+yy.^2));
% 在CCD平面上成的相
Uz = Uf.*phase;
figure(2),imshow(abs(Uz)/max(max(abs(Uz)))),colormap(gray);
title('物光衍射到CCD平面上的振幅')
%% 开始计算参考光的衍射
rou = z0 * h * N /2/L0/L;
Qx  = rou*L0/2/z0;
Qy = Qx;

x = -L/2:L/N:L/2 - L/N;
y = x;
% 网格化数据
[X,Y] = meshgrid(x,y);  

% 参考光的振幅
Ar = max(max(abs(Uz)));
Ur = Ar.*exp(1j*k*(X.*Qx + Y.*Qy));

% Wh = abs(Ur.*conj(Ur)) + abs(Uz.*conj(Uz));

Uh = Ur+Uz;
Wh = Uh.*conj(Uh);  % 干涉强度条纹
Wh = Wh - Ur.*conj(Ur) - Uz.*conj(Uz);
Imax = max(max(Wh));
Ih = uint8(Wh/Imax*255);
figure(3),imshow(Ih),colormap(gray);
imwrite(Ih,'./result/proRe.bmp');
title('生成的全息图');
