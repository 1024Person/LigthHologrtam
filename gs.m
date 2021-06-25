% ------------- 对于GS算法的理解 --------------
%通过光源和目标图像的振幅信息获取要叠加相位，
%图片中没有相位信息，只有物体的振幅信息
%
%所以GS算法的目的是为了获取一幅图像的相位信息
%
%流程：
%1.因为不知道相位信息，所以先随即赋值一个相位信息
%2.对g(x,y)作一次傅里叶变换
%3.已知衍射光斑的强度分布（再现图像的振幅分布）|G'(fx,fy)|的情况下,
% 用已知振幅替换|F(fx,fy)|,保持相位不变
% --------------------------------------------
function O = gs()
clc;
close all;

[fileName,filePath] = uigetfile(['./','*.*'],['输入初始图像'],200,200);


%抽样
image0 = imread([filePath,fileName]);     % 向程序中读取图片
temp = size(image0);
if(length(temp) >= 3)
    % 转成灰度图
    image0 = rgb2gray(image0);
    
end


image0 = im2double(image0);     % 数据类型的转换,归一化
[M,N] = size(image0);           % 获取图像的维度信息

% 从图片中读取的数据，不包括相位信息，
% 加入随机相位信息
PI = 3.14159;
% image1 = image0 光源的振幅和预设图片的振幅相同
image1 = image0;                % 光源振幅
% 这里生成的相位取值范围是[0,2*pi];
phase = 2*PI*rand(M,N);         % 原始的gs算法，对这里的随机相位的生成有着很大的依赖，或者说是敏感
% 如果这里的初始相位生成的不好，那么下面的相位要收敛于预设相位的速度就会慢上很多，所以下面加上了一个负反馈
% image1是光源经过加上一个附加相位之后的值，
image1 = image1.*exp(1j*phase);    % 振幅加上相位信息,image1是复振幅，假设这是Target Intensity，然后通过不断调节相位信息，来获取真正的相位
% 刚刚测试发现，只有对image1取模的时候，才是image0，对image1取实部的时候，就不是image0了
% 也就是说，如果对一个图片加上一个相位的话，这个图片的光强是不会发生变换的，
%

% 计算编码
% 第一次傅里叶逆变换,变换到频域
image2 = ifft2(ifftshift(image1));
% 无论下面怎么对相位做一些变化，结果光强都是不变的
for t = 1:500
    %     disp(['第',num2str(t),'次循环开始']);
    % 迭代判据
    imgangle = angle(image2);        % 取相位
    image = exp(1i*imgangle);
    image = fftshift(fft2(image));   %  还原到空域
    imgabs = abs(image)/max(max(abs(image)));   % 振幅归一化
    % image0中不包含相位信息，这里需要生成相位带有相位信息的全息图
    sim = corrcoef(image0,imgabs);   % 取相关系数，相关系数就是说看看image0和imgabs的线性相关性
    if sim(1,2) >= 0.9995
        %满足条件，跳出循环
        break
    else
        % 开始迭代
        % 单位振幅
        imgangle = angle(image2);  % 取出image2的相位
        % 这里之所以将image2的振幅变成1，是因为我们设定的光源是均匀的，
        % 说通俗点就是处处相等，之所以这样设定有两个原因
        % 1.简单 2.好实现，现实中直接提供激光就好了
        image2 = exp(1i*imgangle); % 将image2的振幅变成1，得到新的image2
        % 傅里叶正变换变换到空域
        % 此时的image3 已经是经过了图中的圈二步骤了，
        image3 = fftshift(fft2(image2));
        % 负反馈振幅,从这里开始往下就看不懂了
        % 首先，计算出来image3，也就是image2的空域光场，
        % 对空域进行取相位，这时候的相位就是进行完圈二之后的相位了，
        % 但是这个相位不一定就是最好的相位，
        imgangle = angle(image3);
        % 丢弃掉image3的振幅信息，因为这里的振幅信息并不是对的
        % 下面会将image0也就是原图的image0的振幅信息重新赋值给image3
        % 注意这里是空域，可以直接将振幅进行赋值，结果也就是空域上的振幅
        image3 = exp(1i*imgangle);
        
        % 按照原来的gs算法，这一步应该是将image3的振幅赋值为预设图片的振幅(image0)
        % 但是博客上说这样收敛速度太慢，就进行了下面的处理
        % 这样得到的image3就是图中的圈三，但是不是最好的圈三
        %         image3 = image3.*(image0 + rand(1,1)*(image0 - imgabs));
        image3 = image3 .* (image0);
        % 进行圈四操作，使用逆傅里叶变换
        image2 = ifft2(ifftshift(image3));
    end
end
% 循环跳出来之后，的结果就是带着相位的image2的模和image0的相关系数达到了0.995及以上，那么就说名这个带着相位信息的image2等效于image0
% 因为最后成像的关键点就是image2的模
% disp('循环出来了！！！')
imgangle = angle(image2);   % 获取image2的相位，就是途中左下角的Hologram
image4 = exp(1i*imgangle);  % 频域图片

% 还原相位图片,将image4还原到空域上去
image4 = fftshift(fft2(image4));
% 归一化，因为是灰度图，而且归一化对整体也没有什么影响
imgabs = abs(image4)/max(max(abs(image4)));
angle_ = angle(image4); % 空域上的相位图片

%  存储位相全息图
% 这一部为什么要这样作？加上pi/(2*pi)
% imgangle = (imgangle + pi) / (2*pi);
% 等会再试试
imgangle = angle(image2);
O = image4;
end
