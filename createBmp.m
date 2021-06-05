close all;


E = zeros(100,100);

E(20:30,20:80) = 1;
E(20:90,20:30) = 1;
E(50:60,20:80) = 1;
E(80:90,20:80) = 1;
figure,imshow(E,[]);
imwrite(E,'./res/E.bmp');