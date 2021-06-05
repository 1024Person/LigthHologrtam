close all;
clc;

% 使用200x200,然后使用resize方法进行更改，最后可以得到比较清晰的像
% 今天晚上可以去试试
U0=imread('./res/E.bmp');
U0=double(U0(:,:,1)); 
r=100; c=100;                   
U0=imresize(U0,[r,c]);         
figure,imshow(U0,[])                      
[r,c]=size(U0);                           
ef1=0.8;
% 这里之所以要乘上一个随机相位我猜测，是因为常数的傅里叶变换是脉冲函数，脉冲函数理论上在0处的函数值为无限大！
% 但是这里即使是加上随机相位也是有讲究的，这里随机相位的范围是  -ef1*pi ~ ef1 * pi
FU0=fftshift(fft2(U0.*exp(1i*rands(r,c)*pi*ef1)));         
Am=abs(FU0);                               
Ph=mod(angle(FU0),2*pi);                 
Ph=Ph./2/pi;                         % 这个不是对相位进行归一化，而是对相位进行编码                      
w=3;                                 % 通光孔径的宽度
wth=round(w/2);                      % 通光孔径宽度的一半
s = 6;
CGH=zeros(r*s,c*s);                     
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%罗曼三型编码
MaxAm=max(abs(FU0(:)));
ef2=1.5;
th=MaxAm/ef2;                               % 他这里没有设置归一化，而是设置了一个阈值！！
Am(Am>th)=th;                               % 阈值是th
lmn=round(Am/th*w);                         % 但是在对Amn进行编码的时候，这里相当于先进行了归一化，然后又开始了乘上通光孔经的宽度
pmn=round(Ph*w);
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% 全息图的赋值
% 在一个编码单元上，通光孔径覆盖的地方，透过率为1，通光孔径没有覆盖的地方为0
for m=1:r
    for n=1:c
        cgh=zeros(s,s);    % 编码单元               
        if lmn(m,n)==0     % 如果这一点的振幅本来就是0那么就不用管了，因为就算透过率为1，也没有，直接将0赋值给迂回相位全息图就好了
        elseif pmn(m,n)<=wth
			cgh(w+1-lmn(m,n):w+lmn(m,n),w-wth+pmn(m,n)+1:2*w-wth+pmn(m,n))=1;
		elseif pmn(m,n)>wth
			cgh(w+1-lmn(m,n):w+lmn(m,n),w-wth+pmn(m,n)+1:s)=1;
			cgh(w+1-lmn(m,n):w+lmn(m,n),1:pmn(m,n)-wth)=1;
        end
        CGH((m-1)*s+1:m*s,(n-1)*s+1:n*s)=cgh;
    end
end
figure;imshow(CGH,[]);                    
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% 我明白了，中间之所以会出现最亮的光点是因为那个0级衍射，
% 然后0及衍射是低频成分，我们在这里还将低频成分移动到了中心，所以中间非常的亮，
% 这里将下面的fftshift去掉就会发现最亮的光点分散到了屏幕的两边
RU=fftshift(ifft2(10000*CGH));                  
RI=RU.*conj(RU);                        
figure;imshow(RI,[]),colormap(pink);
