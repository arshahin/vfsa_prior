function [sample1,count]=sample1_reproducable(xmod,dx,nx,xmin,xmax,randf,count)
ix=2.*(randf(count)-0.5).*nx;
count=count+1;
sample1=xmod+ix.*dx;
while (sample1 < xmin || sample1 >xmax) ;
    ix=2.*(randf(count)-0.5).*nx;
    count=count+1;
    sample1=xmod+ix.*dx;
end
    







