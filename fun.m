function fun=fun(x)

%% x should be in interval [-10,10]
%% function has a global minimum at x=0 and it has a couple of local
%% minimum at x=-6 and x=6

fun=x.^2./100+sin(x-pi./2)+1;


%%%example
%  xx=linspace(-10,10,100);
%  fx=fun(xx);
%  plot(xx,fx)

