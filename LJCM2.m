% ------------------------------近似条件下的夫琅禾费衍射-----------------------------
%   通过参考课本上的讲解，自己实现一下近似条件下的夫琅禾费衍射
%   实验说明：
%       平面光澜上具有中心在坐标原点的矩形孔，wx,wy分别是矩形孔沿着
%       坐标x，y方向的半宽度。光澜被垂直入射的平面波照射，紧贴在孔径
%       后方的物平面光场分布为，门函数相乘
%       通过一个透镜，将平面波变成了球面波
%       然后，在透镜后焦面上，计算光场分布（夫琅禾费衍射）
% ---------------------------------------------------------------------------------------------------
clear all;clc;close all;
% 初始条件
% 定义单位
um = 1e-3;
mm = 1;
nm = 1e-6;

h = 532*nm;     %     波长
L = 1*mm;         %     假设平面光的照射区域为1mm
wx = 10*mm;   %     x半宽度
wy = 10*mm;   %     y半宽度
d = 1000*mm;  %    焦距，
k = 2*pi/h;         %     波数
N = 500;        % 计算的总次数，就是x，y序列的个数
n = 1:N;
% ----------------------------------------------------------开始计算-----------------------------------------------
% 当n-1 = N/2的时候，就会等于0
% 也就是说当n = N/2 + 1的时候，会出现x=0的情况
x = -L/2 + L/N*(n-1);   % 生成x序列，这里就相当于在x轴上进行了抽样
y = x;
Lx = 2*x*wx*pi/h/d;
Lx(N/2+1) = 1;              % 这里的Lx还要放到sinc函数里面，所以这里就不能出现0，
% 经过计算得到N/2 + 1这一项会出现0
sincx = sinc(Lx);
Ly = 2*y*wy*pi/h/d;
Ly(N/2+1) = 1;
sincy = sinc(Ly);
% 光强分布
I = zeros(N,N);
% 这里有个坑
% U = 4*wx*wy*exp(1j*k*(x.*x+y.*y)/2/d).*sincx.*sincx.*sincy.*sincy;
C = (4*wx*wy)^2;
Cx = sincx.*sincx;
Cy = sincy.*sincy;
for p=1:N
    for q = 1:N
        I(p,q) = C*Cx(p)*Cy(q);
    end
end

% 展示在figure上的
figstr = strcat('X向宽度=',num2str(2*wx),'mm,Y向宽度=',...
    num2str(2*wy),'重建物平面宽度=',num2str(L),'mm');

Imax = max(max(I));
p =100;
while p

    strp = num2str(p);

    figure(1),imshow(I,[0,Imax/p]);
    
    title(strcat('Imax/',strp,'限幅显示的衍射斑强度图像'));
    xlabel(figstr);
    p = input('Imax/p,限幅显示，p="?(按Enter结束)');

end
figure(2)

Tx = h*d/wx/2;
strTx = num2str(Tx);
Ix = I(N/2,:);
figure(2),plot(x,Ix);title(strcat('X方向周期=',strTx,'mm'));
xlabel(figstr);


% plot3(x,y,I);
% xlabel('x'),ylabel('y'),zlabel('z');