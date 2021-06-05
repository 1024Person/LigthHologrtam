% ---------------------圆孔夫琅禾费衍射----------------
% 已经手算出来了最后的公式了，
% 最后将那个公式从代码中敲出来。
% 编码注意事项：分母一定不能是0
% ------------------------------------------------------------------
clear all;clc;close all;

% ----------------------- 初始化 ----------------------------------
mm = 1;                        % 定义单位
um = 1e-3;
nm = 1e-6;

w = 2*mm;                      % 圆孔半径% 半径越小，衍射的光斑越大
lambda  = 532*nm;              % 光波波长
d = 1000*mm;                   % 透镜焦距
k = 2*pi/lambda;               % 波数
L = 1*mm;                      % 观察屏的尺寸
N = 500;           
n = 1:N;                       % 抽样序列

x = -L/2 + L/N*(n-1);
y = x;
[yy, xx] = meshgrid(y, x);     % 生成二维网格
r = sqrt(yy.^2+xx.^2);           % 其实就是r

arg = k*w*r/d;                 % 一阶贝塞尔函数的参数

J1 = besselj(1, arg);          % 计算一阶贝塞尔函数
arg(N/2+1,N/2+1) = 1;          % 从上面可以看出来，生成的网格中间是0，
% 那么arg就不可能的避免有一个0
CJ = (2*J1./arg);           % 计算贝塞尔方程的部分
CJ(N/2+1,N/2+1) = 1;           % 另中心值为1，恢复准确值
C = (pi*w^2/lambda/d)^2;       % 前面的常数部分
I = (C*CJ.^2);
Imax = max(max(I));
figstr = strcat('圆孔半径：',num2str(w),'mm,重建物平面宽度=',num2str(L),'mm');
p=1;
while p
    strp = num2str(p);
    s=figure(1);
    imshow(I,[0 Imax/p]);
    title(strcat('Imax/p',strp,'限幅显示的衍射强度图像'));
    xlabel(figstr);
    p = input('Imax/p,限幅显示，p=？（按Enter键结束）');
end
D = 1.22*lambda*d/w;
strD = num2str(D);
Ix = I(N/2,:);
figure,plot(x,Ix);
title(strcat('第一环的直径',strD,'mm'));
xlabel(figstr);