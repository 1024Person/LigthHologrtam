% ------------------MyLJCM_3_2_4---------------------------------------
% 功能：经典衍射公式的D-FFT算法 (就是各种标量衍射公式的离散傅里叶变换计算方法) 
% 执行：输入一张图片输出一张衍射图像
% 主要变量：
%   h   ------ 波长(mm)
%   z0  ------ 衍射距离(mm)
%   U0  ------ 初始光波场的复振幅
%   L0  ------ 初始场及衍射场宽度(mm)
%   Uf  ------ 衍射光波场的复振幅
% ----------------------------------------------------------------------
clear all;clc;close all;

[fileName,filePath] = uigetfile(['./','*.*'],['输入图片'],100,100);
[XRGB,MAP] = imread([filePath,fileName]);   % 读取图片

if(length(size(XRGB)) >= 3)
    X0 = rgb2gray(XRGB);
else
    X0 = XRGB;
end
z0 = input('衍射距离z0=？(mm)');

h = 0.532e-3;   % 波长mm
k = 2*pi/h;     % 波数
[M0,N0] = size(X0);

N1 = max(M0,N0);
% N难道说的是抽样点数？
N = 400;        % ?????????这一步是在干什么？，，，图片的实际宽度吗？
% imresize的第二个参数是缩放！！！扩大或缩小到原来的N/N1倍
X1 = imresize(X0, N/N1);    %  这一步又是在干什么？原文：“修改为最大宽度为N点的图像”，你tm语文没学好吧！！

% 显然这一步中max([M1,N1]) = N,这么做也就保证了下面的赋值操作不会出现下表是负数的情况了！！
% 但是N代表的究竟是什么呢？真的是抽样点数吗？
[M1,N1] = size(X1); % --- M1是四倍的M0
X = zeros(N,N);
X(N/2 - M1/2+1:N/2+M1/2,N/2-N1/2+1:N/2+N1/2) = X1(1:M1,1:N1);  % 修改成N*N的点的图像

U0 = double(X); % 初始场复振幅
% 不是很明白，不应该是抽样频率还有抽样间隔满足初始物光场的宽度吗？为什么这里要通过抽样点数N来确定初始物光场的宽度？！！
L0 = sqrt(h*z0*N);  % 计算同时满足振幅和相位抽样的物光场宽度
L0 = input('初始场宽度L0=？(mm)');

cal = input('角谱1，菲涅尔解析2，菲涅尔FFT 3，基尔霍夫4，瑞利-索末菲5？');
% Uh = 0:N-1 - N/2;Vh = 0:N-1 - N/2;
% [mh,nh] = meshgrid([Uh,Vh]);

figstr = strcat('初始像宽度=',num2str(L0),'mm');
figure(1),imshow(X,[]),colormap(gray);xlabel(figstr);title('初始振幅场分布');

% D-FFT计算
Uf = fft2(U0,N,N);  % 傅里叶变换中的N,N,如果U0不满足NxN的大小的话，那么就在做傅里叶逆变换之前进行补0
Uf = fftshift(Uf);
II = Uf.*conj(Uf);
Isum = sum(sum(II))/N/N*L0/N*L0/N;
Gmax = max(max(II));
Gmin = min(min(II));

figstr = strcat('频谱宽度=',num2str(N/L0),'/mm');

figure(2),imshow(II,[Gmin,Gmax/1000]),colormap(gray);xlabel(figstr);
title('FFT计算方法');

% 下面的switch就是用来计算每一种衍射的传递函数的
switch cal
case 1  % ----------------- 角谱衍射传递函数
    method = '角谱衍射传递函数计算衍射';
    n = 1:N;
    % 看看下面的这个公式中的   “h”  ！！！这就说明了下面的trans中为什么没有波长了
    % 下面生成的频谱取值范围：yy和xx就是h*fy和h*fx！！！但是为什么确定了频谱的取值范围是怎么样的？
    x = h*(-N/L0/2+1/L0*(n-1));
    y = x;
    [yy,xx] = meshgrid(y,x);
    trans = exp(1j*k*z0*sqrt(1-xx.^2-yy.^2));

    f2 = Uf.*trans;

    U = ifft2(f2,N,N); % 对f2做傅里叶逆变换得到衍射距离为z0的角谱衍射
case 2  % -------------- 菲涅尔解析传递函数
    method = '菲涅尔解析传递函数的计算';
    n = 1:N;
    x = h*(-N/L0/2+1+1/L0*(n-1));   
    y = x;
    [yy,xx] = meshgrid(y,x);    % 生成频率序列

    trans = exp(1j*k*z0*(1-(xx.^2+yy.^2)/2));
    f2 = Uf.*trans;
    U = ifft2(f2,N,N);

case 3  % ------------- 菲涅尔传递函数的FFT计算
    method = '菲涅尔传递函数的FFT计算';
    n = 1:N;    % 这里的n是从1开始的！！！

    x = -L0/2+L0/N*(n-1);
    y = x;

    [yy,xx] = meshgrid(y,x);
    f2 = exp(1j*k*z0)/(1j*h*z0);
    f2 = f2.*exp(1j*k/2/z0*(xx.^2+yy.^2));
    f2 = fft2(f2,N,N);
    trans = fftshift(f2);

    f2 = Uf.*trans;
    U = ifft2(f2,N,N);     %  这两步和课本上的不太一样
    U = fftshift(U);      %  课本上的是下面注释括起来的
    % ------------------------------------------------------------
    % Uf = fft2(f2,N,N);    % 对N*N点的离散函数f2座FFT计算
    % U = fftshift(Uf);    % 将FFT计算结果进行整序
    % U = imrotate(U,180);% 旋转180度
    % -----------------------------------------------------------
    % 上面的代码有两个不明白的，首先就是为什么要做傅里叶正变换？
    % 虽然说做傅里叶正变换照样能从频域变换到空域上，但是严格来说就应该是逆变换！
    % 第二个不明白的就是，下面的旋转180度，这又是在干什么？

case 4 % ------------------------------ 基尔霍夫传递函数的D-FFT计算
    method = '基尔霍夫传递函数的计算';

    n = 1:N;
    x = -L0/2+L0/N*(n-1);
    y = x;
    [yy,xx] = meshgrid(y,x);

    f2 = exp(1j*k*sqrt(z0^2+xx.^2+yy.^2));
    f2 = f2./(1j*2*h*(z0^2 + xx.^2+yy.^2));
    f2 = f2.*(sqrt(z0^2+xx.^2+yy.^2));
    f2 = fft2(f2);
    trans = fftshift(f2);

    f2 = Uf.*trans;
    U = ifft2(f2,N,N);
    U = fftshift(U);


    % -------------------------------------------
    % U = fft2(f2,N,N);
    % U = fftshift(U);
    % U = imrotate(U,180);

case 5  % ------------------- 瑞利-索末菲传递函数的D-FFT计算
    method = '瑞利-索末菲传递函数计算衍射';

    n = 1:N;
    x = -L0/2+L0/N*(n-1);
    y = x;
    [yy,xx] = meshgrid(y,x);
    f2 = exp(1j*k*sqrt(z0^2+xx.^2+yy.^2));
    f2 = f2./(1j*h*(z0^2+xx.^2+yy.^2));
    f2 = fft2(f2,N,N);
    trans = fftshift(f2);

    Uf = Uf.*trans;

    U = ifft2(Uf);
    U = ifftshift(U);
    % U = imrotate(U,180);

end
figstr = strcat('衍射场图像宽度',num2str(L0),'mm');
figure(3),imshow(abs(U),[]),colormap(gray);xlabel(figstr);
title(method);imwrite(abs(U), gray,['./result/',method,'_',num2str(z0),'_',num2str(L0),'.bmp']);



trans_ = ifft(trans);
trans_ = fftshift(trans_);
figure(4),imshow(abs(trans_)/max(max(abs(trans_))),[]),colormap(gray);



