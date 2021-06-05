% -----------------------------
% 功能：离轴数字全息
%    
%   
%
%
% -------------------------------------------------------------------------- 从电脑中读入图像文件
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
% -------------------------------------------------------------------------- 将输入的图像调整成NxN分辨率
[M0,N0] = size(X0);
N1 = min([M0,N0]);
N = 1000;   % 模拟形成的全息图的取样数
X1 = imresize(X0, N/4/N1);  % 这里为什么要除以4？搞不明白
[M1,N1] = size(X1);         % 获取缩放之后的大小
X = zeros(N,N);             % 定义一个大的图片
% 将缩放之后的图片放大打的图片中间去
X(N/2-M1/2+1:N/2+M1/2,N/2-N1/2+1:N/2+N1/2) = X1(1:M1,1:N1);
% --------------------------------------------------------------------------- 定义初始化变量
h =input('波长');   % 波长
k = 2*pi/h;     % 波数
pix = input('CCD像素宽度');   % CCD像素宽度(mm)
L = N*pix;      % CCD宽度(mm)
z0 = input('衍射距离：?(mm)');      % 衍射距离
L0 = h*N*z0/L;  % 物平面宽度
% ---------------------------------------------------------------------------- 开始处理
Y = double(X);
b = rand(N,N)*2*pi;
% 添加随机相位是为了让图像的频谱变得更加平滑。这部分应该是图像处理部分的知识，实在是不会啊！
f = Y.*exp(1j*b);   % 添加随机相位噪声，形成振幅正比于图像的初始复振幅
figstr = strcat('初始物平面宽度=',num2str(L0),'mm');
figure(2),imshow(f,[]),colormap(gray);
xlabel(figstr);title('物平面图像');
% ---------------------------------------------------------------------------- 菲涅尔衍射的S-FFT计算开始
n = 1:N;
x= -L0/2+L0/N*(n-1);        % 物光波取样，坐标原点再中心
y = x;
[yy,xx] = meshgrid(y, x);   % 网格化数据
Fresnel = exp(1j*k/2/z0*(xx.^2+yy.^2));
f2 = f.*Fresnel;
Uf = fft2(f2,N,N);
Uf = fftshift(Uf);
% 观察平面
x = -L/2+L/N*(n-1); % CCD宽度取样(mm)
y = x;
[yy,xx] = meshgrid(y, x);   % 网格化数据
% 菲涅尔衍射积分前方的相位因子
phase = exp(1j*k*z0)/(1j*h*z0)*exp(1j*k/2/z0*(xx.^2+yy.^2));
Uz = Uf.*phase;

% ------------------------------------------------------------------------- 菲涅尔衍射的S-FFT计算结束
figstr = strcat('模拟CCD宽度=',num2str(L),'mm');
figure(3),imshow(abs(Uz),[]),colormap(gray);
xlabel(figstr);title('到达CCD平面的物光振幅分布');
% ------------------------------------------------------------------------- 形成0-255灰度级的数字全息图

% 根据下面的公式来看的话，这个Qx,Qy应该就是x，y方向的方向/余弦了 
% Qx = input('方向与县');       % 按照优化设计定义参考光方向余弦
Qx = input('x方向余弦：？')     % 通过测试发现这里的方向余弦选的越小，未来成像的三部分成像分开的距离越大
Qy = Qx;
disp('方向余弦：')
disp(strcat(num2str(Qx),num2str(Qy),'1'));
x = -L/2:L/N:L/2-L/N;       
y = x;
% 网格化数据，这样的得到的是全息平面上的网格坐标
[X,Y] = meshgrid(x, y); 
% 参考光的振幅
Ar = max(max(abs(Uf))); %???????????????????
% 参考光的相位信息
Ur = Ar*exp(1j*k*(X.*Qx+y.*Qy));
% 在空域上进行相干叠加
Uh = Ur+Uz;
% 得到全息图
Wh = Uh.*conj(Uh);
Imax = max(max(Wh));
Ih = uint8(Wh./Imax*255);
imwrite(Ih, './result/Ih.bmp');
figstr = strcat('全息图宽度=',num2str(L),'mm');
figure(4),imshow(Ih,[]),colormap(gray);
xlabel(figstr);title('模拟形成全息图');
