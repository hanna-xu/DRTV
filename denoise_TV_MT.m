function [X_den,R,S,iter]=denoise_TV_MT(Xobs,lambda,l,u,R_init, S_init,pars)
%This function implements the FISTA method for JTV denoising problems.

%Define the Projection onto the box
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

% Assigning parameres according to pars and/or default values
flag=exist('pars', 'var');
if (flag&&isfield(pars,'MAXITER'))
    MAXITER=pars.MAXITER;
else
    MAXITER=100;
end
% if (flag&&isfield(pars,'epsilon'))
%     epsilon=pars.epsilon;
% else
%     epsilon=1e-4;
% end
% if(flag&&isfield(pars,'print'))
%     prnt=pars.print;
% else
%     prnt=1;
% end
if(flag&&isfield(pars,'tv'))
    tv=pars.tv;
else
    tv='iso';
end

[m,n]=size(Xobs);
% clear P; clear R;
if(isempty(R_init))
    R=zeros(m-1,n);    S=zeros(m,n-1);
    U=zeros(m-1,n);    V=zeros(m,n-1);
else
    R=R_init;    S=S_init;
    U=R_init;    V=S_init;
end

tk=1;tkp1=1;
count=0;i=0;
D=zeros(m,n);%fval=inf;fun_all=[];
while((i<MAXITER)&&(count<5))
    %    fold=fval;
    i=i+1;
    %     Dold=D;
    Rold=R;Sold=S;
    tk=tkp1;
    
    temp=Xobs-lambda*Lforward(U, V, m, n);
    D=project(temp);
    
    %%%%%%%%%%
    % Taking a step towards minus of the gradient
    step0=1/(8*lambda);
    
    [tQ1, tQ2]=Ltrans(D, m, n);
    R=U+step0*tQ1;
    S=V+step0*tQ2;
    
    
    %%%%%%%%%%
    % Peforming the projection step
    switch tv
        case 'iso'
            A=[R.^2;zeros(1,n)]+[S.^2,zeros(m,1)];
            A=sqrt(max(A,1));
            R=R./A(1:m-1,:);
            S=S./A(:,1:n-1);
        case 'l1'
            %P{1}=P{1}./(max(abs(P{1}),1));
            %P{2}=P{2}./(max(abs(P{2}),1));
            R=R./(max(abs(R),1));
            S=S./(max(abs(S),1));
            %error('unknown type of total variation. should be iso or l1');
    end
    
    %%%%%%%%%%
    %Updating R and t
    tkp1=(1+sqrt(1+4*tk^2))/2;
    
    step=(tk-1)/tkp1;
    
    U=R+step*(R-Rold);
    V=S+step*(S-Sold);
   
    
    %     re=norm(D-Dold,'fro')/norm(D,'fro');
    %     if (re<epsilon)
    %         count=count+1;
    %     else
    %         count=0;
    %     end
    %     C=Xobs-lambda*Lforward(P, m, n, T);
    %     PC=project(C);
    %     fval=-norm(C-PC,'fro')^2+norm(C,'fro')^2;
    %     fun_all=[fun_all;fval];
end
X_den=D;iter=i;


    function X=Lforward(P1, P2, m, n)
        X=zeros(m,n);
        X(1:m-1,:)=P1;
        X(:,1:n-1)=X(:,1:n-1)+P2;
        X(2:m,:)=X(2:m,:)-P1;
        X(:,2:n)=X(:,2:n)-P2;
    end

    function [P1, P2]=Ltrans(X, m, n)
        %       [m,n]=size(X);
        
        P1=X(1:m-1,:)-X(2:m,:);
        P2=X(:,1:n-1)-X(:,2:n);
        
    end
end

