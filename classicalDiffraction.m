function diffIm = ClassicalDiffraction()
    % ClassicalDiffraction 
    % 经典衍射积分的D-FFT运算
    %
    % Syntax: diffIm = ClassicalDiffraction()
    %
    % 经典衍射的D-FFT计算，正衍射，返回衍射图样U
    
    [fileName,filePath] = uigetfile(['./','*.*'],['输入图片'],100,100);
    [XRGB,MAP] = imread([filePath,fileName]);   % 读取图片
    
    if(length(size(XRGB)) >= 3)
        X0 = rgb2gray(XRGB);
    else
        X0 = XRGB;
    end
    imshow(X0,[]),colormap(gray);
    title('原始图样');
    z0 = input('衍射距离z0=？(mm)');
    
    h = 0.532e-3;   % 波长mm
    k = 2*pi/h;     % 波数
    [M0,N0] = size(X0);
    
    % N1 = max(M0,N0);
    % N难道说的是抽样点数？
    % N = 400;        % ?????????这一步是在干什么？，，，图片的实际宽度吗？
    % imresize的第二个参数是缩放！！！扩大或缩小到原来的N/N1倍
    % X1 = imresize(X0, N/N1);    %  这一步又是在干什么？原文：“修改为最大宽度为N点的图像”，你tm语文没学好吧！！
    
    % 显然这一步中max([M1,N1]) = N,这么做也就保证了下面的赋值操作不会出现下表是负数的情况了！！
    % 但是N代表的究竟是什么呢？真的是抽样点数吗？
    % [M1,N1] = size(X1); % --- M1是四倍的M0
    % X = zeros(N,N);
    % X(N/2 - M1/2+1:N/2+M1/2,N/2-N1/2+1:N/2+N1/2) = X1(1:M1,1:N1);  % 修改成N*N的点的图像
    
    U0 = double(X0); % 初始场复振幅
    % 不是很明白，不应该是抽样频率还有抽样间隔满足初始物光场的宽度吗？为什么这里要通过抽样点数N来确定初始物光场的宽度？！！
    Lx = sqrt(h*z0*N0);  % 计算同时满足振幅和相位抽样的物光场宽度
    Ly = sqrt(h*z0*M0);
    % L0 = input('初始场宽度L0=？(mm)');
    
    cal = input('角谱1，菲涅尔解析2，菲涅尔FFT 3，基尔霍夫4，瑞利-索末菲5？');
    
    % D-FFT计算
    Uf = fft2(U0,M0,N0);  % 傅里叶变换中的N,N,如果U0不满足NxN的大小的话，那么就在做傅里叶逆变换之前进行补0
    Uf = fftshift(Uf);  % 将低频成分移动到中间
    % 下面的switch就是用来计算每一种衍射的传递函数的
    switch cal
    case 1  % ----------------- 角谱衍射传递函数
        method = '角谱衍射传递函数计算衍射';
        n = 1:N0;
        m = 1:M0;
        % 看看下面的这个公式中的   “h”  ！！！这就说明了下面的trans中为什么没有波长了
        % 下面生成的频谱取值范围：yy和xx就是h*fy和h*fx！！！但是为什么确定了频谱的取值范围是怎么样的？
        x = h*(-N0/Lx/2+1/Lx*(n-1));
        y = h*(-M0/Ly/2+1/Ly*(m-1));
        [yy,xx] = meshgrid(y,x);
        trans = exp(1j*k*z0*sqrt(1-xx.^2-yy.^2));
        f2 = Uf.*trans;
        U = ifft2(f2,M0,N0); % 对f2做傅里叶逆变换得到衍射距离为z0的角谱衍射
        U = fftshift(U);
    case 2  % -------------- 菲涅尔解析传递函数
        method = '菲涅尔解析传递函数的计算';
        n = 1:N0;
        m = 1:M0;
        x = h*(-N0/Lx/2+1/Lx*(n-1));   
        y = h*(-M0/Ly/2+1/Ly*(m-1));
        [yy,xx] = meshgrid(y,x);    % 生成频率序列
    
        trans = exp(1j*k*z0*(1-(xx.^2+yy.^2)/2));
        f2 = Uf.*trans;
        U = ifft2(f2,M0,N0);
        U = fftshift(U);
    case 3  % ------------- 菲涅尔传递函数的FFT计算
        method = '菲涅尔传递函数的FFT计算';
        n = 1:N0;    % 这里的n是从1开始的！！！
        m = 1:M0;
    
        x = -Lx/2+Lx/N0*(n-1);
    
        y = -Ly/2+Ly/M0*(m-1);
    
        [yy,xx] = meshgrid(y,x);
        f2 = exp(1j*k*z0)/(1j*h*z0);
        f2 = f2.*exp(1j*k/2/z0*(xx.^2+yy.^2));
        f2 = fft2(f2,N,N);
        trans = fftshift(f2);
    
        f2 = Uf.*trans; % f2./trans;  trans - fftshift(f2) - fft2(f2) - f2 exp
        U = ifft2(f2,M0,N0);     %  这两步和课本上的不太一样
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
    
        n = 1:N0;
        m = 1:M0;
    
        x = -Lx/2+Lx/N0*(n-1);
        y = -Ly/2+Ly/M0*(m-1);
        [yy,xx] = meshgrid(y,x);
    
        f2 = exp(1j*k*sqrt(z0^2+xx.^2+yy.^2));
        f2 = f2./(1j*2*h*(z0^2 + xx.^2+yy.^2));
        f2 = f2.*(sqrt(z0^2+xx.^2+yy.^2));
        f2 = fft2(f2);
        trans = fftshift(f2);
    
        f2 = Uf.*trans;
        U = ifft2(f2,M0,N0);
        U = fftshift(U);
    
    case 5  % ------------------- 瑞利-索末菲传递函数的D-FFT计算
        method = '瑞利-索末菲传递函数计算衍射';
    
        n = 1:N;
        m = 1:M;
    
        x = -Lx/2+Lx/N0*(n-1);
        y = -Ly/2+Ly/M0*(m-1);
        [yy,xx] = meshgrid(y,x);
        f2 = exp(1j*k*sqrt(z0^2+xx.^2+yy.^2));
        f2 = f2./(1j*h*(z0^2+xx.^2+yy.^2));
        f2 = fft2(f2,M0,N0);
        trans = fftshift(f2);
    
        Uf = Uf.*trans;
    
        U = ifft2(Uf);
        U = fftshift(U);
    
    end
        diffIm = U;
    end