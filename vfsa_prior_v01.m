clc
close all
clear

% randf=rand(10^7,1); %%% this random vector is already generated and sa
% save randdata randf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
randtype=0; %%%% random(non-reproducable results)
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
prior_type=1; 
%%%%%%%% 1 for A unimodal Prior with low uncertainty :  sharp Gaussian N(0,1)
%%%%%%%% 2 for A unimodal Prior with high uncertainty: broad Gaussian N(0,3)
%%%%%%%% 3 for A deviated unimodal Prior simulating : broad Gaussian N(2,3)
%%%%%%%% 4 for A equiprobable bimodal Prior : N(0,1)+N(-4,1)
%%%%%%%% 5 for A bimodal Prior simulating high probability on target : N(0,1)+0.5*N(-4,1)
%%%%%%%% 6 for A bimodal Prior simulating low probability on target : 0.5*N(0,1)+N(-4,1)
%%%%%%%% 7 for An erroneous unimodal Prior  : N(6,1)

if prior_type==1
    prior1=pdf('Normal',intv,0,1);
    cdf_prior1=cumsum(prior1);
elseif prior_type==2
    prior2=pdf('Normal',intv,0,3);
elseif prior_type==3
    prior3=pdf('Normal',intv,2,3);
elseif prior_type==4
    prior4=pdf('Normal',intv,0,1)+pdf('Normal',intv,-4,1) ;
elseif prior_type==5
    prior5=pdf('Normal',intv,0,1)+0.5.*pdf('Normal',intv,-4,1) ;
elseif prior_type==6
    prior6=0.5.*pdf('Normal',intv,0,1)+pdf('Normal',intv,-4,1) ;
elseif prior_type==7
    prior7=pdf('Normal',intv,6,1);
else
    print('prior_type is not the library')
end

cdf_prior1=cumsum(prior1);
cdf_prior2=cumsum(prior2);
cdf_prior3=cumsum(prior3);
cdf_prior4=cumsum(prior4);
cdf_prior5=cumsum(prior5);
cdf_prior6=cumsum(prior6);
cdf_prior7=cumsum(prior7);

cp1=max(cdf_prior1);
cp2=max(cdf_prior2);
cp3=max(cdf_prior3);
cp4=max(cdf_prior4);
cp5=max(cdf_prior5);
cp6=max(cdf_prior6);
cp7=max(cdf_prior7);%%% should be equal to 1

if cp1~=1; cdf_prior1=cdf_prior1./max(cdf_prior1);  prior1(1)=cdf_prior1(1);for i=2:N; prior1(i)=cdf_prior1(i)-cdf_prior1(i-1); end; end;
if cp2~=1; cdf_prior2=cdf_prior2./max(cdf_prior2);  prior2(1)=cdf_prior2(1);for i=2:N; prior2(i)=cdf_prior2(i)-cdf_prior2(i-1); end; end;
if cp3~=1; cdf_prior3=cdf_prior3./max(cdf_prior3);  prior3(1)=cdf_prior3(1);for i=2:N; prior3(i)=cdf_prior3(i)-cdf_prior3(i-1); end; end;
if cp4~=1; cdf_prior4=cdf_prior4./max(cdf_prior4);  prior4(1)=cdf_prior4(1);for i=2:N; prior4(i)=cdf_prior4(i)-cdf_prior4(i-1); end; end;
if cp5~=1; cdf_prior5=cdf_prior5./max(cdf_prior5);  prior5(1)=cdf_prior5(1);for i=2:N; prior5(i)=cdf_prior5(i)-cdf_prior5(i-1); end; end;
if cp6~=1; cdf_prior6=cdf_prior6./max(cdf_prior6);  prior6(1)=cdf_prior6(1);for i=2:N; prior6(i)=cdf_prior6(i)-cdf_prior6(i-1); end; end;
if cp7~=1; cdf_prior7=cdf_prior7./max(cdf_prior7);  prior7(1)=cdf_prior7(1);for i=2:N; prior7(i)=cdf_prior7(i)-cdf_prior7(i-1); end; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cooling scheduale for VFSA  
nruns=1;
maxtemps=40;
nmov=3;
temp0=1;
decay=1.0;
t0=1;
nx=round(abs(xmax-xmin)./dx)+1;
%%%%%%%%%%%%%%%%Prior type
prior=prior1;
cdf_prior=cdf_prior1;
% prior_add=0; %%%% 0 for pure VFSA, 1 for constant 0.25, 2 for constant 0.5; 3 for linear increase, ....
% %                 4 for linear decrease, 5 for poly, 6 for exp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% weigths on priors 
mm=4.0;%%% power for polynomial weigths, it could be greater or smaller than one
a=0.0; %%% amount of weigth at iter==1
b=0.5; %%% amount of weigth at iter==N
c=0.65; %%% decay for exponential weigth
N=maxtemps; %%% will be used for polynomial and linear weigths


for k=1:maxtemps
    tmp1(k)=t0.*exp(-decay.*(k-1).^0.5);
    w1(k)=(b-a).*(1-exp(-c.*(k-1).^0.5))+a; %%%exponential
    w2(k)=(b-a).*((k-1)./(N-1))+a; %%% Linear inceasing
    w22(k)=(a-b).*((k-1)./(N-1))+b; %%% Linear deceasing
    w3(k)=(b-a).*((k-1)./(N-1)).^mm+a; %% Polynomial
    w4(k)=0; %%% pure VFSA
    w5(k)=0.25; %%% fixed weigth
    w6(k)=0.50; %%% fixed weigth    
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
%     xmod=-9.0;
    emod=fun(xmod);
    
    
    for jtemp=1:maxtemps
        temp(jtemp)=temp0.*exp(-decay.*(jtemp-1).^0.5);        
        tmp=t0.*exp(-decay.*(jtemp-1).^0.5);
        
        
%       w=0.0; %%% pure VFSA
        w=(b-a).*((jtemp-1)./(N-1)).^mm+a; %%% Polynomial  
%           w=(b-a).*((jtemp-1)./(N-1))+a; %%% Linear increasing
%         w=(a-b).*((jtemp-1)./(N-1))+b; %%% Linear decreasing
%         w=0.25; %%% fixed weigth 0.25
%         w=(b-a).*(1-exp(-c.*(jtemp-1).^0.5))+a; %%%exponential
%         w=0.50; %%% fixed weigth 0.50

%         
        for jmov=1:nmov
%             xtrial=walk(xmod,dx,xmin,xmax,tmp);
%            [xtrial,pdf(:,jtemp),cdf(:,jtemp)]=walk_model(xmod,dx,xmin,xmax,tmp);
            methodmc=1;
           if randtype==0
               
               [xtrial,pdf(:,jtemp),cdf(:,jtemp)]=pdf_combine(prior,nx,xmin,xmax,xmod,tmp,w,methodmc);
%             xtrial=pdf_combine(prior2,nx,xmin,xmax,xmod,tmp,w);
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
% ss2=sum(merror)


fs=16;

figure(1);
plot(xx,funx,'LineWidth',2);
title('Objective function','FontSize',fs);xlabel('X','FontSize',fs);ylabel('F(X)','FontSize',fs);set(gca,'FontSize',fs)

figure(2);
plot(intv,prior1,'b','LineWidth',2); hold on;plot(intv,prior2,'g','LineWidth',2); hold on; plot(intv,prior3,'c','LineWidth',2)
hold on;plot(intv,prior4,'r','LineWidth',2);hold on;plot(intv,prior5,'y','LineWidth',2);hold on;plot(intv,prior6,'k','LineWidth',2);
ylabel('PDF','FontSize',fs);
xlabel('X','FontSize',fs)
title('Prior PDF','FontSize',fs)
% ylim([0 0.02])
% grid
legend('N(0,1)','N(0,3)','N(2,3)','N(0,1)+N(-4,1)','N(0,1)+0.5*N(-4,1)','0.5*N(0,1)+N(-4,1)')
set(gca,'FontSize',fs)


figure(3);
plot(intv,cdf_prior1,'b','LineWidth',2); hold on;plot(intv,cdf_prior2,'g','LineWidth',2); hold on; plot(intv,cdf_prior3,'c','LineWidth',2); 
hold on;plot(intv,cdf_prior4,'r','LineWidth',2);hold on;plot(intv,cdf_prior5,'y','LineWidth',2);hold on;plot(intv,cdf_prior6,'k','LineWidth',2);
ylabel('CDF','FontSize',fs);
xlabel('X','FontSize',fs)
legend('N(0,1)','N(0,3)','N(2,3)','N(0,1)+N(-4,1)','N(0,1)+0.5*N(-4,1)','0.5*N(0,1)+N(-4,1)')
title('Prior CDF','FontSize',fs)
% grid
set(gca,'FontSize',fs)


figure(33)
plot(tmp1,'--k','LineWidth',2);
hold on;
plot(w1,'r','LineWidth',2);hold on;plot(w2,'b','LineWidth',2);hold on;plot(w22,'--b','LineWidth',2);hold on;plot(w3,'-oc','LineWidth',2);hold on;plot(w4,'k','LineWidth',2);hold on;
plot(w5,'m','LineWidth',2);hold on;plot(w6,'g','LineWidth',2)
legend('Temp','Weigth');
xlim([1 maxtemps]);
ylim([-0.10 1]);
legend('Temp','W(Exponential)','W(Linear increase)','W(Linear decrease)','W(Polynomial[4])','Pure VFSA(W=0)','Fixed W=0.25','Fixed W=0.50')
title('Temp and weigths vs. iteration','FontSize',fs);xlabel('Iteration No.','FontSize',fs);ylabel('Temp&Weigths','FontSize',fs);set(gca,'FontSize',fs)


figure(4);
plot(ermod1,'b','LineWidth',2);
hold on
% plot(merror,'-or','LineWidth',2);
% text(10,0.9,'o-o Mean','color','r','FontSize',fs) 
% text(10,1,'--- Absolute','color','b','FontSize',fs) 
xlim([1 maxtemps]);
% ylim([0 max(fun(xx))]);
% title('Absolute error"200 runs"','FontSize',fs);
xlabel('Iteration No.','FontSize',fs);
ylabel('Absolute error','FontSize',fs);
set(gca,'FontSize',fs)

% 
% save purevfsa ermod1

% savefile='merror_purevfsa_sharpgaussian.mat';
% save(savefile, 'merror');
% load merror_purevfsa_sharpgaussian;

% figure(6);

% xlim([1 maxtemps]);
% ylim([0 1]);
% title('Mean of absolute errors','FontSize',fs);xlabel('Iteration No.','FontSize',fs);ylabel('Absolute error','FontSize',fs);set(gca,'FontSize',fs)



% ermod2=reshape(ermod1,nruns*maxtemps,1); 
% figure(7);
% hist(ermod2,20); ylim([0 nruns*maxtemps]); xlim([0 3]); 
% title('Absolute error histogram','FontSize',fs);xlabel('Absolute error','FontSize',fs);ylabel('Frequency','FontSize',fs);set(gca,'FontSize',fs) 

figure(5);
hist(m2); %ylim([0 160]); xlim([-0.40 0.40]); 
title('Optimum model','FontSize',fs);xlabel('Optimum model','FontSize',fs);ylabel('Frequency','FontSize',fs);set(gca,'FontSize',fs) 


% pdf1=permute(pdf(:,1,1:nruns),[1 3 2]);
% pdf10=permute(pdf(:,10,1:nruns),[1 3 2]);
% pdf20=permute(pdf(:,20,1:nruns),[1 3 2]);
% pdf30=permute(pdf(:,30,1:nruns),[1 3 2]);
% 
% 
% cdf1=permute(cdf(:,1,1:nruns),[1 3 2]);
% cdf10=permute(cdf(:,10,1:nruns),[1 3 2]);
% cdf20=permute(cdf(:,20,1:nruns),[1 3 2]);
% cdf30=permute(cdf(:,30,1:nruns),[1 3 2]);


pdf1=pdf(:,1);
pdf10=pdf(:,10);
pdf20=pdf(:,20);
pdf30=pdf(:,30);
% pdf40=pdf(:,40);
% pdf50=pdf(:,50);


cdf1=cdf(:,1);
cdf10=cdf(:,10);
cdf20=cdf(:,20);
cdf30=cdf(:,30);
% cdf40=cdf(:,40);
% cdf50=cdf(:,50);

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
% legend('Iter1','Iter10','Iter20','Iter30','Prior');
title('Dynamic CDF','FontSize',fs)
set(gca,'FontSize',fs)


%%%%%%%%%%%%% Summarizing sum(errors) table

% sum_error =
%        Prior1 Prior2 Prior3 Prior4 Prior5 Prior6 Prior7
% w=0         743   817   716   744   722   669     734
% w(poly)     687   756   714   721   707   695     715
% w(linear)   587   695   741   681   637   663    1492
% w(exp)      382   598   605   466   434   543    3805
% w=0.25      325   561   562   418   412   442    3963
% w=0.5       246   455   513   388   313   434    4123

% sum_error=[743 817 716 744 722 669 734; 687 756 714 721 707 695 715; 587 695 741 681 637 663 1492; 382 598 605 466 434 543 3963;325 561 562 418 412 442 3805; 246 455 513 388 313 434 4123];
% sum_error=[743 817 716 744 722 669; 687 756 714 721 707 695; 587 695 741 681 637 663; 382 598 605 466 434 543;325 561 562 418 412 442; 246 455 513 388 313 434];
% sum_error=[743 743 743 743 743 743; 687 756 714 721 707 695; 587 695 741 681 637 663; 382 598 605 466 434 543;325 561 562 418 412 442; 246 455 513 388 313 434];



% figure(8);
% plot(sum_error,'LineWidth',2);
% ylabel('Error summation','FontSize',fs);
% xlabel('W=0       W(polynomial)       W(linear)            W(exp)           W=0.25          W=0.50','FontSize',fs)
% % legend('N(0,1)','N(0,3)','N(2,3)','N(0,1)+N(-4,1)','N(0,1)+0.5*N(-4,1)','0.5*N(0,1)+N(-4,1)','N(6,1)')
% legend('N(0,1)','N(0,3)','N(2,3)','N(0,1)+N(-4,1)','N(0,1)+0.5*N(-4,1)','0.5*N(0,1)+N(-4,1)')
% title('Summation of errors','FontSize',fs)
% grid
% set(gca,'FontSize',fs)
%  

% mean_poly=mean(687,756,714,721,707,695,715);

    











