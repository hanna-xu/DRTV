clear all;close all;
addpath(genpath('.'));

lambda=50;
imagepairs=1:20;

time=[];
iter=300;
for i=1:length(imagepairs)
    pic=num2str(i)
    ImageIR=imread(strcat('IR\',pic,'.bmp'));
    ImageVis=imread(strcat('VIS\',pic,'.bmp'));

    tic
    ImageIR = 255*im2double(ImageIR);
    ImageVis = 255*im2double(ImageVis);
    divK = 2;
    
    [m1,n1]=size(ImageVis);
    m=fix(m1/divK);
    n=fix(n1/divK);
    ImageIR=imresize(ImageIR,[m,n]);
    
    [ImageFus,Obj(i,:)]= SIRF(ImageIR,ImageVis,divK,lambda,iter);
    
    close all;
    ImageVis=uint8(ImageVis);
    ImageIR=uint8(ImageIR);
    ImageFus=uint8(ImageFus);
    
    toc
    time=[time toc];
    subplot(131),imshow(ImageVis);
    subplot(132),imshow(ImageIR);
    subplot(133),imshow(ImageFus);
    imwrite(ImageFus,strcat('results\',pic,'.bmp'));
end