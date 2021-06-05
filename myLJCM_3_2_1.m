% ----------------- 通过自己的理解写出来S-FFT算法
%   根据菲涅尔衍射的傅里叶计算公式，我们需要的参数有：
%       lambda------- 波长
%       d   --------- 衍射距离
%       U0  --------- 初始光波的复振幅
%       Uf  --------- 衍射光波长的复振幅
%       L   --------- 衍射场宽度
% -------------------------------------------

clear all;
clc;
close all;

% 定义单位
mm = 1;
um = 1e-3;
nm = 1e-6;

rootPath = './res/';
[nom,chemin] = uigetfile([rootPath,'*.*'],['输入初始图像'],100,100);
[XRGB,MAP] = imread([chemin,nom]);
X = rgb2gray(XRGB);
z0 = 1000*mm;
z0 = input('请输入衍射距离=？(单位：mm)');
lambda = 532*nm;    % 波长
k = 2*pi/lambda;     % 波数
[M,N] = size(X);
N = max([M,N]);
X = imresize(X, [N,N]);
U0 = double(X);
L0 = sqrt(lambda*z0*N);  % 计算出同时满足振幅和相位抽样定理的物光场的宽度
% L0 = input('物广场') 这里感觉不应该改变啊！！否则的话就会不能同时满足抽样定理了
% 频率空间上的抽样，就是教材上的p，q
Uh = [0:N-1] - N/2;Vh = [0:N-1] - N/2;
[mh,nh] = meshgrid(Uh,Vh);
figstr = strcat('初始像宽度=',num2str(L0),'mm');
subplot(121),imshow(X,[]),colormap(gray);
ylabel('衍射计算');xlabel(figstr);title('S-FFT方法计算衍射');

% ----------------------------------- 开始菲涅尔衍射积分运算
n = 1:N;
% 这一步的作用应该是被做傅里叶变换的函数的抽样吧
% 抽样间距是L0/N
x = -L0/2+L0/N*(n-1);   
y = x;
[yy,xx] = meshgrid(y,x);    % 生成抽样网格
% phase = exp(1j*k*z0);   % 
Fresnel = exp(1j*k/2/z0*(xx.^2+yy.^2));  % 计算傅里叶变换部分的相位部分
f2 = U0.*Fresnel;       % 整个傅里叶变换部分
Uf = fft2(f2);
Uf = fftshift(Uf);
L = lambda*z0*N/L0;          % 计算衍射场的宽度
[yy,xx] = meshgrid(y,x);% 衍射场空域的抽样网格
% 计算傅里叶变换前面的部分
phase = exp(1j*k*z0) / (1j*lambda*z0)*exp(1j*k/2/z0);

Uf = Uf.*phase; % 整个的菲涅尔衍射积分
% 感觉下面的T加不加影响不大
T = L0/N;       % 二位离散变换量值补偿？！什么意思
Uf = Uf*T*T;
% ---------------------------------- 菲涅尔衍射计算的S-FFT计算结束
If = Uf.*conj(Uf);  % 形成衍射场强度分布
figstr = strcat('衍射场图像宽度',num2str(L),'mm');
subplot(122);
imshow(If,[]);

ylabel('衍射计算');xlabel(figstr);title('S-FFT计算衍射');


