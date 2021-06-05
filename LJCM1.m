% ------------------ 矩形孔的夫琅禾费衍射 --------------------

clear all;close all;clc;

um = 1e-3;  % 单位

h = 0.532*um;  % 波长
k = 2*pi / h;
N = 500;
L0 = 1;         % 不明白为什么通光孔径都要比平面光波还要大
wx = 10;
wy = 10;
d = 1000;

% ------------------------------------------------------------

n = 1:N;  % 空域序列
x = -L0/2+L0/N*(n-1);   % x空域序列
y = x;                                  % y等于x
Lx = 2*pi*wx/h/d*x;     % 这里代码是对的，但是课本上是错的，课本上的参数中没有加上 \pi
Lx(N/2+1) = 1;                 %
sincx = sin(Lx)./Lx;
% sincx(N/2+1) = 1;
Ly = 2*pi*wy/h/d*y;
Ly(N/2+1) = 1;
sincy = sinc(Ly)./Ly;
% sincy(N/2+1) = 1;
fx = sincx.*sincx;
fy = sincy.*sincy;
I = zeros(N,N);
C = (4*wx*wy/h/d)^2;

for p=1:N
    for q = 1:N
        I(p,q) = C*fx(p)*fy(q);
    end
end

figstr = strcat('X向宽度=',num2str(2*wx),'mm,Y向宽度=',...
    num2str(2*wy),'重建物平面宽度=',num2str(L..0),'mm');
Imax = max(max(I));
p = 10;
while p
    strp = num2str(p);
    figure(1);imshow(I,[0,Imax/p]);
    title(strcat('Imax/',strp,'限幅显示的衍射斑强度图像'));
    xlabel(figstr);
    p = input('Imax/p,限幅显示，p="?(按Enter结束)');

end
Tx = h*d/wx/2;
strTx = num2str(Tx);
Ix = I(N/2,:);
figure(2),plot(x,Ix);title(strcat('X方向周期=',strTx,'mm'));
xlabel(figstr);
