% -----------------------------
% 功能：离轴数字全息
%    
%
%   
%
%
% -------------------------------
clear all;close all;clc;
rootPath = './res/'; % 设置根路径
[fileName,filePath] = uigetfile([rootPath,'*.*'],'输入图像',100,100);
[XRGB,MAP] = imread([filePath,fileName]);
if length(size(XRGB)) > 2
    X0 = rgb2gray(XRGB);
else
    X0 = XRGB;
end

figure(1),imshow(X0,[]),colormap(gray);
xlabel('原始图片');title('原始图片');
% --------------------------------------------------------- 将输入的图像调整成NxN分辨率
[M0,N0] = size(X0);
N1 = min([M0,N0]);
N = 500;   % 模拟形成的全息图的取样数
X1 = imresize(X0, N/4/N1);  % 这里为什么要除以4？搞不明白
[M1,N1] = size(X1);         % 获取缩放之后的大小
X = zeros(N,N);             % 定义一个大的图片
% 将缩放之后的图片放大打的图片中间去
X(N/2-M1/2+1:N/2+M1/2,N/2-N1/2+1:N/2+N1/2) = X1(1:M1,1:N1);
% -------------------------------------------------------------------- 定义初始化变量
h = 0.532e-3;   % 波长
k = 2*pi/h;     % 波数
pix = 0.0046;   % CCD像素宽度(mm)
L = N*pix;      % CCD宽度(mm)
z0 = 1000;      % 衍射距离
L0 = h*N*z0/L;  % 物平面宽度
% ------------------------------------------------------------------------ 开始处理
% Y = double(X);
% b = rand(N,N)*2*pi;
% 添加随机相位是为了让图像的频谱变得更加平滑。这部分应该是图像处理部分的知识，实在是不会啊！
% f = Y.*exp(1j*b);   % 添加随机相位噪声，形成振幅正比于图像的初始复振幅
% figstr = strcat('初始物平面宽度=',num2str(L0),'mm');
% figure(2),imshow(Y,[]),colormap(gray);
% xlabel(figstr);title('物平面图像');
% --------------------------------------------------------- 菲涅尔衍射的S-FFT计算开始
n = 1:N;
x0= h*(-L0/2+L0/N*(n-1));
y0 = x0;
[fx,fy] = meshgrid(x0, y0);
trans = exp(1j*k*z0*(1-(fx.^2+fy.^2)/2));
Uf = fft2(X,N,N);
UF = Uf.*trans;
Ar = ones(size(UF));

x = -L0/2+L0/N*(n-1);
y = x;
[xx,yy] = meshgrid(x, y);
theta = pi/20;
Rf = Ar.*exp(-1j*k*sin(theta)*yy);

I = (Rf+UF).*conj(Rf+UF);   % 计算出光强分布
imshow(I,[]),colormap(gray);
xlabel(strcat('衍射角度\pi/20，衍射距离：',num2str(z0)));
title('离轴全息图');



% [yy,xx] = meshgrid(y, x);
% Fresnel = exp(1j*k/2/z0*(xx.^2+yy.^2));
% f2 = X.*Fresnel;
% Uf = fft2(f2,N,N);
% Uf = fftshift(Uf);
% x = -L/2+L/N*(n-1); % CCD宽度取样(mm)
% y = x;
% [yy,xx] = meshgrid(y, x);
% 菲涅尔衍射积分前方的相位因子
% phase = exp(1j*k*z0)/(1j*h*z0)*exp(1j*k/2/z0*(xx.^2+yy.^2));
% Uf = Uf.*phase;
% --------------------------------------------------------- 菲涅尔衍射的S-FFT计算结束
% figstr = strcat('模拟CCD宽度=',num2str(L),'mm');
% figure(3),imshow(abs(Uf),[]),colormap(gray);
% xlabel(figstr);title('到达CCD平面的物光振幅分布');
% ------------------------------------------------------- 形成0-255灰度级的数字全息图

% fex = N/L;
% Qx = (4-2.5)*L0/8/z0;       % 按照优化设计定义参考光方向余弦
% Qy = Qx;
% x = -L/2:L/N:L/2-L/N;       % 看不懂
% y = x;
% [X,Y] = meshgrid(x, y);
% Ar = max(max(abs(Uf)));
% Ur = Ar*exp(1j*k*(X.*Qx+y.*Qy));
% Uh = Ur+Uf;
% Wh = Uh.*conj(Uh);
% Imax = max(max(Wh));
% Ih = uint8(Wh./Imax*255);
% imwrite(Ih, './result/Ih.tif');
% figstr = strcat('全息图宽度=',num2str(L),'mm');
% figure(4),imshow(Ih,[]),colormap(gray);
% xlabel(figstr);title('模拟形成全息图');
