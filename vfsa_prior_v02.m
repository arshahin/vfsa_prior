clc
close all
clear

% randf=rand(10^7,1); %%% this random vector is already generated and sa
% save randdata randf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
randtype=1; %%%% random(non-reproducable results)
% randtype=1; %%%% reproducable results 

if randtype==1
    load randdata;
    count=1;
end

xmin=-10;  xmax=10;  dx=0.05;

%%% Interval
intv=xmin:dx:xmax;
N=max(size(intv));

%%%%%%%%%%%%%%%%%%%%%%%%%% objective function
xx=xmin:dx:xmax;
funx=fun(xx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Prior PDFs types 
prior_type=2; 
%%%%%%%% 1 for A unimodal Prior with low uncertainty :  sharp Gaussian N(0,1)
%%%%%%%% 2 for A unimodal Prior with high uncertainty: broad Gaussian N(0,3)
%%%%%%%% 3 for A deviated unimodal Prior simulating : broad Gaussian N(2,3)
%%%%%%%% 4 for A equiprobable bimodal Prior : N(0,1)+N(-4,1)
%%%%%%%% 5 for A bimodal Prior simulating high probability on target : N(0,1)+0.5*N(-4,1)
%%%%%%%% 6 for A bimodal Prior simulating low probability on target : 0.5*N(0,1)+N(-4,1)
%%%%%%%% 7 for An erroneous unimodal Prior  : N(6,1)

if prior_type==1
    prior=pdf('Normal',intv,0,1);    
elseif prior_type==2
    prior=pdf('Normal',intv,0,3);
elseif prior_type==3
    prior=pdf('Normal',intv,2,3);
elseif prior_type==4
    prior=pdf('Normal',intv,0,1)+pdf('Normal',intv,-4,1) ;
elseif prior_type==5
    prior=pdf('Normal',intv,0,1)+0.5.*pdf('Normal',intv,-4,1) ;
elseif prior_type==6
    prior=0.5.*pdf('Normal',intv,0,1)+pdf('Normal',intv,-4,1) ;
elseif prior_type==7
    prior=pdf('Normal',intv,6,1);
else
    print('prior_type is not the library')
end

cdf_prior=cumsum(prior);
cp=max(cdf_prior);
if cp~=1; cdf_prior=cdf_prior./max(cdf_prior); 
    prior(1)=cdf_prior(1);
    for i=2:N 
        prior(i)=cdf_prior(i)-cdf_prior(i-1); 
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cooling scheduale for VFSA  
nruns=1;
maxtemps=40;
nmov=3;
temp0=1;
decay=1.0;
t0=1;
nx=round(abs(xmax-xmin)./dx)+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% weigths on priors 
weight_type=2;
%%%%%%%% 1 for pure VFSA (no prior is applied)
%%%%%%%% 2 for Polynomial increase (with increasing the iteration number weights are increasing with a 4th order Polynomial)
%%%%%%%% 3 linear incease 
%%%%%%%% 4 for exponential
%%%%%%%% 5 fixed weigth 0.50 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mm=4.0;%%% power for polynomial weigths, it could be greater or smaller than one
a=0.0; %%% amount of weigth at iter==1
b=0.5; %%% amount of weigth at iter==N
c=0.65; %%% decay for exponential weigth
N=maxtemps; %%% will be used for polynomial and linear weigths


for k=1:maxtemps
    tmp1(k)=t0.*exp(-decay.*(k-1).^0.5);
    if weight_type==1
        wtemp=zeros(size(tmp1)); %%%% pure VFSA 
    elseif weight_type==2
        wtemp(k)=(b-a).*((k-1)./(N-1)).^mm+a; %% Polynomial
    elseif weight_type==3
        wtemp(k)=(b-a).*((k-1)./(N-1))+a; %%% Linear inceasing
    elseif weight_type==4
        wtemp(k)=(b-a).*(1-exp(-c.*(k-1).^0.5))+a; %%%exponential
    elseif weight_type==5
        wtemp(k)=0.50; %%% fixed weigth
    else
        print('weight_type is not the library')
    end       
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VFSA loop
temp=zeros(maxtemps,1);

for j=1:nruns
    ermod=0;
    if randtype==0
        xmod=sample1(xmin,dx,nx,xmin,xmax);
    elseif randtype==1
        [xmod,count]=sample1_reproducable(xmin,dx,nx,xmin,xmax,randf,count);
    else
        sprintf('randtype must be either 1 or zero')
        return
    end

    emod=fun(xmod);
    
    
    for jtemp=1:maxtemps
        temp(jtemp)=temp0.*exp(-decay.*(jtemp-1).^0.5);        
        tmp=t0.*exp(-decay.*(jtemp-1).^0.5);
        
        
        
        if weight_type==1
            w=0; % pure VFSA 
        elseif weight_type==2
            w=(b-a).*((jtemp-1)./(N-1)).^mm+a; %%% Polynomial
        elseif weight_type==3
            w=(b-a).*((jtemp-1)./(N-1))+a; %%% Linear increasing
        elseif weight_type==4
            w=(b-a).*(1-exp(-c.*(jtemp-1).^0.5))+a; %%%exponential
        elseif weight_type==5
            w=0.50; %%% fixed weigth 0.50
        else
            print('weight_type is not in the library')
        end
        
        for jmov=1:nmov
           methodmc=1;
           if randtype==0               
               [xtrial,pdf(:,jtemp),cdf(:,jtemp)]=pdf_combine(prior,nx,xmin,xmax,xmod,tmp,w,methodmc);
           else
               [xtrial,pdf(:,jtemp),cdf(:,jtemp),count]=pdf_combine_reproducable(prior,nx,xmin,xmax,xmod,tmp,w,randf,count,methodmc);
           end
            xtrial=round((xtrial-xmin)./dx).*dx+xmin;
            etrial=fun(xtrial);
            
            if etrial< emod
                emod=etrial;
                xmod=xtrial;                
            else
                arg=(etrial-emod)./temp(jtemp);
                if arg>1.e6
                    pde=0.001;
                else
                    pde=exp(-arg);
                end
                if randtype==0
                    if pde>rand
                        emod=etrial;
                        xmod=xtrial;
                    end
                else
                    if pde>randf(count)
                        emod=etrial;
                        xmod=xtrial;
                    end
                    count=count+1;
                end
            end              
            
        end %%%% end move
        ermod1(jtemp,j)=emod;
%         a(jtemp)=jtemp;
%         b(jtemp)=temp(jtemp);
        m1(jtemp,j)=xmod;       
%         e(jtemp)=emod;
    end
    qq=min(ermod1(:,j));
    best=find(qq==ermod1(:,j));
    if max(size(best))>1
        best=round(mean(best));
    end
    m2(j)=m1(best,j);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures
ss1=sum(sum(ermod1));
merror=mean(ermod1');

fs=16;

figure(1);
plot(xx,funx,'b','LineWidth',2);
title('Objective function','FontSize',fs);xlabel('X','FontSize',fs);ylabel('F(X)','FontSize',fs);set(gca,'FontSize',fs)

figure(2);
plot(intv,prior,'b','LineWidth',2); 
ylabel('PDF','FontSize',fs);
xlabel('X','FontSize',fs)
title('Prior PDF','FontSize',fs)
% ylim([0 0.02])
grid
set(gca,'FontSize',fs)


figure(3);
plot(intv,cdf_prior,'b','LineWidth',2); 
ylabel('CDF','FontSize',fs);
xlabel('X','FontSize',fs)
title('Prior CDF','FontSize',fs)
grid
set(gca,'FontSize',fs)


figure(4)
plot(tmp1,'b','LineWidth',2);
hold on;
plot(wtemp,'r','LineWidth',2);
legend('Temp','Weigth');
xlim([1 maxtemps]);
title('Temp and weigths vs. iteration','FontSize',fs);
xlabel('Iteration No.','FontSize',fs);
ylabel('Temp&Weigths','FontSize',fs);
set(gca,'FontSize',fs)


figure(5);
plot(ermod1,'b','LineWidth',2);
hold on
xlim([1 maxtemps]);
xlabel('Iteration No.','FontSize',fs);
ylabel('Absolute error','FontSize',fs);
set(gca,'FontSize',fs)




pdf1=pdf(:,1);
pdf10=pdf(:,10);
pdf20=pdf(:,20);
pdf30=pdf(:,30);



cdf1=cdf(:,1);
cdf10=cdf(:,10);
cdf20=cdf(:,20);
cdf30=cdf(:,30);


figure(6);
plot(intv,pdf1,'b','LineWidth',2);hold on;plot(intv,pdf10,'g','LineWidth',2); hold on; plot(intv,pdf20,'c','LineWidth',2); hold on; plot(intv,pdf30,'r','LineWidth',2);
hold on; plot(intv,prior,'--k','LineWidth',2);
ylabel('PDF','FontSize',fs);
ylim([0 0.08])
xlabel('X','FontSize',fs)
title('Dynamic PDF','FontSize',fs)
legend('Iter1','Iter10','Iter20','Iter30','Prior');
set(gca,'FontSize',fs)

figure(7);
plot(intv,cdf1,'b','LineWidth',2);hold on;plot(intv,cdf10,'g','LineWidth',2); hold on; plot(intv,cdf20,'c','LineWidth',2); hold on; plot(intv,cdf30,'r','LineWidth',2);
hold on; plot(intv,cdf_prior,'--k','LineWidth',2);
ylabel('CDF','FontSize',fs);
xlabel('X','FontSize',fs)
legend('Iter1','Iter10','Iter20','Iter30','Prior');
title('Dynamic CDF','FontSize',fs)
set(gca,'FontSize',fs)













