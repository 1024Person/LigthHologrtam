% 罗曼三型编码

clc;clear all;
% 读取图片到程序中
U = imread('./res/E.bmp');
U = double(U(:,:,1));    % 将图片转换成double类型的图片
r = 100;c = 100;
U = imresize(U,[r,c]);
figure(1);
imshow(U,[]);
[r,c] = size(U);
ef1 =1;         % 这里的ef1越大，再现的图像就越清晰，这里可以深究一下，

temp = U.*exp(1i*rand(r,c)*2*pi*ef1);
imshow(temp),colormap(gray);title('随机相位图');

% 这里的思考方向可以是在时域上的乘上一个相位然后频域上对应着的操作，我记得好像是频移！！！
UF = fftshift(fft2(U.*exp(1i*rand(r,c)*2*pi*ef1)));
s =6;      % 设置编码单元的大小
w =s;      % 设置通光孔径的宽度,通光孔径设置成s，最合适！！，只是测试了一下整数，因为W控制着透过率，如果W越高，那么透过率就越大，透过的能量就越多，显示的就越清楚
wth = round(w/2);   % 通光孔经宽度的一半
Am = abs(UF);       %
Ph = angle(UF);     % 取相位
Ph = mod(Ph,2*pi);  % 取2pi的余数
Ph = Ph ./ 2 / pi;   % 对相位进行编码
maxAm = max(Am(:));    % 最大的Am
% 下面简单的对振幅进行归一化,通过实践发现，下面进行简单的归一化是非常不可取的！！！
% Am = Am / maxAm;
% 但是这里设置一个阈值又有什么道理呢？
% 如果只是简单的归一化的话，最后形成的图像会什么也不是！！！
ef2 =1.5;  % 滤波操作，这里对于每一个图片的阈值设置不一样，我读进来一个200x200规格的图片，对这个图片的滤波操作最好在4
th = maxAm / ef2;
Am(Am > th) = th;
CGH = zeros(r*s,c*s);      % 全息图矩阵
% -------------------- 罗曼三型编码 ----------------
Lmn = round(Am/th*s/2);       % 通光孔径高度的一半
Pmn = round(Ph*s/2);       % 编码中心距离通光孔径中心距离的一半
% ------------------ 对全息图进行赋值 ---------------

for i = 1:r
    
    for j = 1:c
        cgh = zeros(s,s);      % 编码单元
        if Lmn(i,j) == 0   % 如果振幅为0的话，透过率无论是多少，都没有用了，所以就不用管振幅为0的透过率了
        
        elseif Pmn(i,j) <= wth  % 如果没有发生模式溢出的话
            % 通光孔径的透过范围
            cgh(s/2 - Lmn(i,j) + 1:s/2+Lmn(i,j),s/2 + Pmn(i,j) - wth +1 :s - wth + Pmn(i,j)) = 1;
        elseif Pmn(i,j) > wth   % 发生模式溢出
            cgh(s/2 - Lmn(i,j) + 1: s/2 +Lmn(i,j),s/2 + Pmn(i,j) - wth:s) = 1;
            cgh(s/2- Lmn(i,j)  + 1: s/2 +Lmn(i,j),1:Pmn(i,j) - wth) = 1;   
        end
        CGH((i-1)*s+1:i*s,(j-1)*s+1:j*s) = cgh;
    end
end
disp('迂回相位全息图生成完成！');
figure(2);imshow(CGH);

% ------------------------ 再现部分 ----------------------------

% 这里不再将低频成分移动到中间了显然低频成分就是直流分量
% 但是这里操作是不可以的，在这里是再现部分了，
% 我们要做的就是改善全息图的生成质量
RU=fftshift(ifft2(CGH));
RI=RU.*conj(RU);
figure(3);imshow(RI,[]);
colormap(pink)
