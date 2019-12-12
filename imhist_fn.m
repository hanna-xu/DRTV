function [hist_out]=imhist_fn(x)

[p,q]=size(x);
hist_out=zeros(1,256);
for ii=1:p
   for jj=1:q
      hist_out(x(ii,jj)+1)=hist_out(x(ii,jj)+1)+1;
   end
end