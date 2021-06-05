% ------------------ 计算经典衍射积分的逆运算
% 通过之前生成的衍射图样，将之输入到程序中来，
% 就能计算出最开始的原始图样了
% 主要变量：
%   h  ------ 波长
%   L  ------ 衍射场宽度(mm)
%   L0 ------ 初始场的宽度(mm)
%   z0 ------ 衍射距离(mm)
% -------------------------------------------
clear all;clc;close all;

[fileName,filePath] = uigetfile(['./res/','*.*'],'输入衍射图片',200,200);

[XRGB,MAP] = imread([filePath,fileName]);

if length(size(XRGB))>=3

    X0 = rgb2gray(XRGB);
else
    X0 = XRGB;
end
% X0便是衍射之后的图片
% 将X0转换成double类型的
X0 = double(X0);



h = 532e-6; % 波长
k = 2*pi/h; % 波数

[M0,N0] = size(X0);
N1 = max([M0,N0]);
% N是抽样点数
N = 400;
N = input('抽样点数=？');
X1 = imresize(X0, N/N1);

[M1,N1] = size(X1);

X = zeros(N,N);
% 将衍射图样调整到NxN的画布上
X(N/2 - M1/2 + 1:N/2+M1/2,N/2-N1/2+1:N/2+N1/2) = X1(1:M1,1:N1);
z0 = 1000;
z0 = input('衍射距离=？(mm)');
L = sqrt(h*z0*N);   % 光场宽度
L = input('衍射场宽度L0=？(mm)');





cal = input('衍射计算公式：\n1、角谱衍射D-FFT，2、菲涅尔衍射S-FFT 3、菲涅尔衍射D-FFT');

X = double(X);
Uf = fft2(X,N,N);
Uf = fftshift(Uf);

switch cal
case 1  % ----------- 角谱衍射的D-FFT逆运算

    method = '角谱衍射D-FFT逆运算';

    n = 1:N;
    % 在这里就先将波长乘到频率序列中，下面trans就不至于括号那么多了
    x = h*(-L/2+L/N*(n-1));
    y = x;
    [xx,yy] = meshgrid(x,y);
    trans = exp(-1j*k*z0*sqrt(1-xx.^2-yy.^2));
    FF = Uf.*trans; % 初始平面的傅里叶空间表达式
case 2  % ---------- 菲涅尔衍射的S-FFT逆运算

    method = '菲涅尔衍射的S-FFT逆运算';
    n = 1:N;
    x = (-L/2+L/N*(n-1));
    y = x;
    [xx,yy] = meshgrid(x,y);
    C1 = exp(-1j*k*z0)/(-1j*h*z0);
    C2 = exp(-1j*k/(2*z0)*(xx.^2+yy.^2));
    C3 = exp(-1j*k/2/z0*(xx.^2+yy.^2));
    trans = ifft2(Uf.*C3);
    U = C1*C2.*trans;

case 3  % ---------- 菲涅尔衍射的D-FFT逆运算
    method = '菲涅尔衍射的D-FFT逆运算';
    n = 1:N;
    x = h*(-L/2+L/N*(n-1));
    y = x;
    [yy,xx] = meshgrid(y,x);
    trans = exp(-1j*k*z0*(1-(xx.^2+yy.^2)/2));
    FF = Uf.*trans; % 初始平面的傅里叶空间数值解
end


figstr = strcat('衍射场宽度：',num2str(L),'mm  衍射距离：',num2str(z0),'mm');

figure(1),imshow(X,[]),colormap(gray);
xlabel(figstr),title('衍射图样');


if cal~=2
    U = ifft2(FF,N,N);
end
I = U.*conj(U);
figure(2),imshow(I,[]),colormap(gray);
xlabel(method),title('初始图样');
