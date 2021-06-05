% ---------------------------------------------------------------
%   离轴全息图的复现工作 1-FFT重建
%       
%       输入数字全息图，
%       输入重建距离，可以通过大于1的参数p（rou）的选择进行限幅放大显示
%       重建图像
%
%   主要变量：
%       h ----- 波长                Th ----- 数字全息图
%       L ----- 数字全息图的尺寸      L0 ----- 重建平面宽度
%       z0 ---- 记录全息图的距离      U0 ----- 重建物平面光波长的复振幅
% ---------------------------------------------------------------

clear all;clc;close all;
% ---------------------------------------------- 读入图片
fileFillter = './result/*.*';
[fileName,filePath] = uigetfile(fileFillter,...
'选择数字全息图',100,100);
I1 = imread([filePath,fileName]);

figure(1);imshow(I1);title('数字全息图');
% ---------------------------------------------- 初始化变量
f0 = double(I1);
[N1,N2] = size(f0);
N = min(N1,N2);

h = input('波长：？(mm)');   % 波长(mm)
pix = input('CCD像素宽度：?(mm)');  % CCD像素宽度
% z0 = 500;      % 重建距离
z0 = input('重建距离z0=?(mm)');
L = pix*N;          % CCD的物理宽度
Ih= zeros(N,N);
% 这里就相当于是经过了一次窗口函数
Ih(1:N,1:N) = f0(1:N,1:N); % 感觉这里有错误！！！
% --------------------------------------------- 1-FFT重建开始
n = 1:N;
x = -L/2+L/N*(n-1);
y = x;
[fx,fy] = meshgrid(x,y);
k = 2*pi/h;
Fresnel = exp(1j*k/2/z0*(fx.^2+fy.^2));
f2 = Ih.*Fresnel;
Uf = fft2(f2,N,N);
Uf = fftshift(Uf);
L0 = h*z0*N/L;
x = -L0/2+L0/N*(n-1);
y = x;
[xx,yy] = meshgrid(x,y);
phase = exp(1j*k*z0)/(1j*h*z0)*exp(1j*k/2/z0*(xx.^2+yy.^2));
U0 = Uf.*phase;
% ----------------------------------------------- 1-FFT重建结束

If = U0*conj(U0);
Gmax = max(max(abs(U0)));
Gmin = min(min(abs(U0)));
figstr = strcat('重建物平面宽度',num2str(L0),'mm');
figure(2),imshow(abs(U0),[Gmin,Gmax/1]),colormap(gray);
xlabel(figstr);title('1-FFT物平面重建');
p = 10;
while p
    figure(3);
    imshow(abs(Uf),[Gmin,Gmax/p]),colormap(gray);
    xlabel(figstr);
    title('1-FFT物平面重建图像');
    p = input('Gmax/p,p=10?');
end;
