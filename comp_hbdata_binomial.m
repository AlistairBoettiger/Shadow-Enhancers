clear all;

load hbSD_010611

figure(1); clf; 

%     cc11  12  13  14 %
Ncc = [125,250,500,1000]; % expected number of on cells.  
cc = {cc11; cc12; cc13; cc14}; 

for c = 1:4
    mu = zeros(1,G);
    sigma = zeros(1,G);
    bi_sig = zeros(1,G); 
    plot_miss = cell(1,G); 
     for k=1:G;     
         plot_miss{k} = foff{k}(cc{c}{k}) ;%foff{k}(cc14{k});
         mu(k) = mean(plot_miss{k})*Ncc(c);
         sigma(k) = std(plot_miss{k})*Ncc(c);
         bi_sig(k) = mu(k)*( 1- nanmean(plot_miss{k}) );
     end
     rels = logical(1-isnan(sigma));

     corc = corrcoef(bi_sig(rels), sigma(rels)   )  ;
     
     D = nanmedian((sigma - bi_sig)./bi_sig);
     
     disp(sigma);
     disp(bi_sig);
     
     figure(1); subplot(2,2,c); scatter(sigma,bi_sig); 
     lin = [0,bi_sig,1.1*max(bi_sig)];
     hold on; plot(lin,lin,'k');
     xlim([0,1.1*max(sigma)]);
     ylim([0,1.1*max(bi_sig)]);
     title(['hb nucs: ', num2str(Ncc(c)), '  corr= ',num2str(corc(1,2),3), '  dist = ',num2str(D)]); 
 
end

set(gcf,'color','w');


%% Simulate binomial Response across different times. 

V = 200; % number of concentration points to check
time = 20 ; % min in cell cylce
Ts = time*60/10; % number of time points to check


D = 4.5;% diffusion rate of bcd (according to Dostatni)
cG = 4.8; % c Gregor

  N = 500;

C = linspace(2,7,V);  % range of concentrations to explore; 
      Th = 1:1000; % range of thresholds number of molecules to explore
      
    low = cG-cG*.1; % 10% less than threshold bcd concentration
    high = cG + cG*.1; % 10% more than threshold bcd concentration

    miss_rate = zeros(1,Ts); 

 figure(4); clf;   
ploop = [.22, .13, .22*.13]; col = {'blue','green','red'}; 
a_effs = [.03, .05, .03+.05]/1.9; 

for e = 1:3; 
    a = a_effs(e); % effective enhancer size
    
    for t = 1:Ts
        T = t*10;% 7*60;
         phalf = 1-gammainc(a*cG*D*T,floor(Th),'upper');
         [err,thresh] = min( abs(phalf - .5) );
         pint =  1-gammainc(a*C*D*T,floor(thresh),'upper');
        [err,cut] = min( abs(C - low));  % find where lines intersect
        miss_rate(t) = pint(cut); 
       
        
%         figure(1); clf; plot(C,pint); hold on; plot([low,low],[0,1],'c');
%         hold on; plot([high,high],[0,1],'r');
%         figure(1); hold on; plot([C(cut)], [pint(cut)],'r.')
    end

%     figure(2); clf; plot(miss_rate, 'k.'); 
%     bvar = N*miss_rate.*(1-miss_rate);
%     figure(2); clf; plot(bvar); 
%     figure(3); clf; hist(miss_rate,0:.01:1); xlim([0,1]);

  
    bf = zeros(Ts,N+1);



    
    for t=1:Ts
        bf(t,:) = binopdf(0:1:N,N,1-(1-miss_rate(t))*(1-ploop(e) ) );
    end

    xbf = linspace(0,1,N+1);
    bfD = sum(bf);
    bfD = bfD/sum(bfD)*N;

    figure(4); hold on; plot(xbf,bfD,'color',col{e});
end
legend('no primary','no shadow','control');
title(['cc13 T = ',num2str(time), 'min  N = ',num2str(N), ' cells']); 
set(gcf,'color','w');


%sum(bfD)*1/N


%% Dostatni numbers

V = 100;


%770/(.003*1*4.8)  % time to count all 770 molecules (takes 891 min, error 3%)
% 1/sqrt(X) = .1 -> X = 100 molecules need to count to have 10% error.
% 100/(.003*1*4.8) --> 115 minutes (Thomas numbers)
% 100/(.003 * 4.5 *4.8) --> still 25 minutes, Dostatni numbers
% 100/(.003*10 * 4.5 *4.8) my enhancer size estimate 2.5 minutes


a = .003;
T = 25*60;% 7*60;
D = 4.5;
c = 4.8;
BP =    1/sqrt(D*a*c*T);

    theta = 93;% 25; % 

C = linspace(2, 7, V); 

intp = zeros(1,V);
for k=1:V;
    c= C(k);  % c = 2 
    lambda = a*c*D*T;
    n = 0:100;
    p = lambda.^n./factorial(n)*exp(-lambda);
%      figure(2); clf; 
%      plot(p,'linewidth',4); set(gcf,'color','k'); colordef black;
%      set(gca,'color','k');
%     sum( p(logical(1-isnan(p)))  );


   %  sum(p(1:theta))  % the manual integration method
   %  gammainc(a*c*D*T,floor(theta),'upper') % the cdf method
    intp(k) = 1-gammainc(a*c*D*T,floor(theta),'upper');
end

h1 = figure(3); clf; set(gca,'color','w'); colordef white; 
plot(C,intp); hold on; plot([4.8,4.8],[0,1],'k'); 
plot([4.3,4.3],[0,1],'r'); 
legend('prob detecting > \theta','boundary conc.','10% less than boundary', 'Location','NorthWest' ); 
xlabel('c, molecules / um^3');
ylabel(['probability of seeing >', num2str(theta),' molecules in time T']);
title(['a = ',num2str(a,3), 'um  T = ',num2str(T,3),'s  D = ',num2str(D,2), 'um^2/s', '  1/(DacT)^{1/2} = ' num2str(BP,3)]);
set(gcf,'color','w');



 print('-depsc', '-tiff', [fout,'dostatni_calc_25min.eps']);

%%  multiple binding sites new effective receptor size of primary enhancer

% factor of pi from binding model version of Berg-Purcell

a = .003;
T = 7*60;
c = 4.8;
b = a*24;
m=6;
D = 4.5;

dcM = 1/sqrt(pi*D*c*T).*sqrt(1/(m*a) + 1/(2*b))

%dcD = 1/sqrt(D*c*T*a_eff) = dcM
a_eff = 1/(dcM*sqrt(D*c*T))^2


%%
V = 100;

a = a_eff;
T = 7*60;
D = 4.5;
c = 4.8;

BP =    1/sqrt(D*a*c*T);
C = linspace(2, 7, V); 
intp = zeros(1,V);
for k=1:V;
    c= C(k); 
    lambda = a*c*D*T;
    theta = 430; 
    intp(k) =  1- gammainc(lambda,floor(theta),'upper');
end

h2 = figure(3); clf; plot(C,intp); hold on; plot([4.8,4.8],[0,1],'k'); 
plot([4.3,4.3],[0,1],'r'); 
legend('prob detecting > \theta','boundary conc.','10% less than boundary', 'Location','NorthWest' ); 
xlabel('c, molecules / um^3');
ylabel(['probability of seeing >', num2str(theta),' molecules in time T']);
title(['a = ',num2str(a,3), 'um  T = ',num2str(T,3),'s  D = ',num2str(D,2), 'um^2/s', '  1/(DacT)^{1/2} = ' num2str(BP,3)]);
set(gcf,'color','w');

% saveas(h2,[fout,'aeff_calc.eps'],'eps');

print('-depsc', '-tiff', [fout,'aeff_calc2.eps']);

%% Effective size of shadow

a = .003;
T = 7*60;
c = 4.8;
b = a*20;  % the 3 bcd sites are only 200 bp apart
m=4;
D = 4.5;

dcM = 1/sqrt(pi*D*c*T).*sqrt(1/(m*a) + 1/(2*b))

%dcD = 1/sqrt(D*c*T*a_eff) = dcM
a_eff = 1/(dcM*sqrt(D*c*T))^2


a = a_eff;
T = 7*60;
D = 4.5;
c = 4.8;

BP =    1/sqrt(D*a*c*T);
C = linspace(2, 7, V); 
intp = zeros(1,V);
for k=1:V;
    c= C(k); 
    lambda = a*c*D*T;
    theta = 295; 
    intp(k) =  1- gammainc(lambda,floor(theta),'upper');
end

h2 = figure(3); clf; plot(C,intp); hold on; plot([4.8,4.8],[0,1],'k'); 
plot([4.3,4.3],[0,1],'r'); 
legend('prob detecting > \theta','boundary conc.','10% less than boundary', 'Location','NorthWest' ); 
xlabel('c, molecules / um^3');
ylabel(['probability of seeing >', num2str(theta),' molecules in time T']);
title(['a = ',num2str(a,3), 'um  T = ',num2str(T,3),'s  D = ',num2str(D,2), 'um^2/s', '  1/(DacT)^{1/2} = ' num2str(BP,3)]);
set(gcf,'color','w');

% saveas(h2,[fout,'aeff_calc.eps'],'eps');

print('-depsc', '-tiff', [fout,'shadow_calc.eps']);




a = .003;
T = 7*60;
c = 4.8;
b = a*4;  % the 3 bcd sites are almost back to back
m=3;
D = 4.5;

dcM = 1/sqrt(pi*D*c*T).*sqrt(1/(m*a) + 1/(2*b))

%dcD = 1/sqrt(D*c*T*a_eff) = dcM
a_eff_Driever = 1/(dcM*sqrt(D*c*T))^2

