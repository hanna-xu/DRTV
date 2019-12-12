function [Y,Obj]= SIRF(IR,Vis, divK, lambda, iter)
%%%% min |RX-M|_2^2/(2*lambda) + |X-T(P)|_VTV
% Infr : Infrared image
% Vis : Visible image
% divK: the resolution difference between the Vis and Infr 
% lambda: regularization parameter
% iter: the number of iterations one does to minimize the energy


% Output----
% Y: the fused image

%%% TV settings
parsin.MAXITER=8; 
parsin.tv='iso';

l = min(min(Vis)); u =max(max(Vis));
if((l==-Inf)&&(u==Inf))
    project=@(x)x;
elseif (isfinite(l)&&(u==Inf))
    project=@(x)(((l<x).*x)+(l*(x<=l)));
elseif (isfinite(u)&&(l==-Inf))
    project=@(x)(((x<u).*x)+((x>=u)*u));
elseif ((isfinite(u)&&isfinite(l))&&(l<u))
    project=@(x)(((l<x)&(x<u)).*x)+((x>=u)*u)+(l*(x<=l));
else
    error('lower and upper bound l,u should satisfy l<u');
end

IR = double(IR);
Vis = double(Vis);
[m,n] = size(Vis);
[c1]=find(isnan(Vis));
Vis(c1) = zeros(size(c1));
L=1;    %Lipschitz constant
x = zeros(m,n);
x = imresize(IR,[m,n]);
y = x;
tnew=1;
Obj=zeros(1,iter);

for k=1:3
    told=tnew;
    xp=x;
    dd = imresize(y,[m/divK,n/divK]) - IR;
    df2 = imresize(dd,[m,n]);
    yg = y - df2/L;
    [z, R, S]=denoise_TV_MT(yg-Vis, 2*lambda/L,-inf,inf,[],[], parsin);
    x = project(z+Vis);
    tnew=(1+sqrt(1+4*told^2))/2;
    y=x+((told-1)/tnew)*(x-xp);
    %imshow(uint8(x),[]);
    
    %%   formulate the Obj value
    fx(k) = 1/2*sum(sum((imresize(x,[m/divK,n/divK]) - IR).^2));
    [P,Q]=size(Vis);
    tmpx1=zeros(P,Q);     tmpx2=zeros(P,Q);
    tmpv1=zeros(P,Q);     tmpv2=zeros(P,Q);
    tmpx1(1:P-1,:)=x(1:P-1,:)-x(2:P,:);
    tmpx2(:,1:Q-1)=x(:,1:Q-1)-x(:,2:Q);
    tmpv1(1:P-1,:)=Vis(1:P-1,:)-Vis(2:P,:);
    tmpv2(:,1:Q-1)=Vis(:,1:Q-1)-Vis(:,2:Q);
    gradx=sqrt(tmpx1.^2 + tmpx2.^2);
    gradv=sqrt(tmpv1.^2 + tmpv2.^2);
    gx(k) = lambda*sum(sum(abs(gradx-gradv)));
    Obj(k)=fx(k)+gx(k);
end


while abs((Obj(k)-Obj(k-1)))/Obj(k-1)>=0.003
    k=k+1;
    told=tnew;
    xp=x;
    
    dd = imresize(y,[m/divK,n/divK]) - IR;
    df2 = imresize(dd,[m,n]);
    yg = y - df2/L;
   [z, R, S]=denoise_TV_MT(yg-Vis, 2*lambda/L,-inf,inf,R, S,parsin);
    
    x = project(z+Vis);
    tnew=(1+sqrt(1+4*told^2))/2;
    y=x+((told-1)/tnew)*(x-xp);
    %imshow(uint8(x),[]);
    
    %%   formulate the Obj value
    fx(k) = 1/2*sum(sum((imresize(x,[m/divK,n/divK]) - IR).^2));
    [P,Q]=size(Vis);
    tmpx1=zeros(P,Q);     tmpx2=zeros(P,Q);
    tmpv1=zeros(P,Q);     tmpv2=zeros(P,Q);
    tmpx1(1:P-1,:)=x(1:P-1,:)-x(2:P,:);
    tmpx2(:,1:Q-1)=x(:,1:Q-1)-x(:,2:Q);
    tmpv1(1:P-1,:)=Vis(1:P-1,:)-Vis(2:P,:);
    tmpv2(:,1:Q-1)=Vis(:,1:Q-1)-Vis(:,2:Q);
    gradx=sqrt(tmpx1.^2 + tmpx2.^2);
    gradv=sqrt(tmpv1.^2 + tmpv2.^2);
    gx(k) = lambda*sum(sum(abs(gradx-gradv)));
    Obj(k)=fx(k)+gx(k);
    % % % MSEaf=sum(sum((imresize(IR,[m,n])-x).^2))/(m*n);
    % % % MSEbf=sum(sum((Vis-x).^2))/(m*n);
    % % % MSE=(MSEaf+MSEbf)/2;
    % % %
    % % % %% PSNR
    % % % Obj(k)=10*log10((max(max(x)))^2/MSE);
end
k
Y = x;
