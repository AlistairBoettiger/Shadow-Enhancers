%%                          Plot_DotComps2_hb_regs.m                             %%
%
% Analyzing New Shadow data
%
%
% Alistair Boettiger                                   Date Begun: 03/05/10
% Levine Lab                                     Functional Since: 03/06/10
%                                                   Last Modified: 11/19/10

%% Description
% comparison
%
%
%% Updates
% Modified 10/18/10 to incldue repression of ectopic expression. 
% Modified to hack failed combine data code
% Modified 11/02 for figures
% Modified 11/16 to split hb into 3 regions. 
% New Code adapted from Plot_DotComps2_hb.m on 11/17 to segment by regions.
%  uses additional code divde_reds to split hb domain into 3rds and analyze
%  by region.  
% Modified 11/19 added more statistical signficance testing lines.  
%% Source Code
clear all;

folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed';

emb_roots = {'MP09_22C_y_hb'; % 1
            'MP02_22C_y_hb'; % 2
            'MP02_30C_y_hb';% 3
            'MP02_30C_LacZ_hb';% 4   % yes it's actually yellow
            'MP01_22C_y_hb'; % 5 
            'BAC09_22C_y_hb'; % 6 
            'BAC09_30C_y_hb';% 7 
            'BAC01_30C_y_hb';% 8
            'BAC02_22C_y_hb' 
            };  
          %1 +6, 2+9, 3+4, 

    

     % all names
names = {'2 enhancers, 22C';
         'no shadow, 22C';
         'no shadow, 30C';
         'no shadow b, 30C'; 
         'no primary, 22C';
         '2 enh 22C';
         '2 enh 30C';
         'no primary 30C';
         'no shadow b, 22C'
         };




N = 100;
K = length(emb_roots); 
G= length(names);



miss_rate1 = cell(1,K); 
miss_rate2 = cell(1,K); 
miss_rate3 = cell(1,K); 
nd = cell(1,K); 
age_table = cell(1,K);

for z=1:K
    miss_rate1{z} = zeros(N,1);
    miss_rate2{z} = zeros(N,1);
    miss_rate3{z} = zeros(N,1);
    nd{z} = zeros(N,1);
end

xmin = .2; xmax = .9; ymin = .15; ymax = .4;
% as fractions of the original image dimensions.  

%%
for z=1:K % k=2;
    for n=   1:N
        if n<10
            emb = ['0',num2str(n)];
        else
            emb = num2str(n);
        end

        try
        load([folder,'/',emb_roots{z},emb,'_data.mat']);   
             [miss_rate1{z}(n),miss_rate2{z}(n),miss_rate3{z}(n)] =  divide_regs(L2,H,pts1,pts2,ptr_nucin2,handles.In);
           if length(H) > 2000
               im_dim = 2048;
           else
               im_dim = 1024;
           end
               
           lims =  round([xmin,xmax,ymin,ymax]*im_dim); 
           nd{z}(n) = NucDensity(cent,lims,0);
            age_table{z}{n,1} = [folder,'/',emb_roots{z},emb,'_data.mat']; %  
            age_table{z}{n,2} = nd{z}(n); 
        catch ME
            disp(ME.message); 
            %disp(['can not find file' folder,'/',emb_roots{z},emb,'_data.mat']);
        end
    end
end

close all; 

% save hb_SD-region_data; 
% save hb_SD-10-21-10;% includes ectopic repression analysis
% save hb_SD-10-18-10;% includes ectopic repression analysis
% save hb_SD-10-13-10;
%  save hb_SD-9-13-10
 %  save hb_shadow_yellow_data;   % load hb_shadow_yellow_data;

 % clear all; load hb_SD-10-21-10;
  %1 +6, 2+9, 3+4, 
  

%    [miss_cnt,miss_rate,nd,ectop_rate] = merge_data(2,9,N,miss_cnt,miss_rate,nd, ectop_rate);
%       [miss_cnt,miss_rate,nd,ectop_rate] = merge_data(1,7,N,miss_cnt,miss_rate,nd, ectop_rate);
%   [miss_cnt,miss_rate,nd,ectop_rate] = merge_data(3,4,N,miss_cnt,miss_rate,nd, ectop_rate);


%%
 % clear all; load hb_SD-region_data; 

ND = cell2mat(nd); 
age_offset = 4.8;

emb_cycle = age_offset + log2( nonzeros( sort(ND(:)) ) );
figure(2); clf; plot( emb_cycle ,'r.');


title(['hb embryos, N = ',num2str(length(nonzeros(ND(:))) )  ],'FontSize',15);
set(gca,'FontSize',15); grid on;
set(gcf,'color','w'); ylabel('log_2(nuc density)');  xlabel('embryo number'); 
ylim([10,14.99]);


%%
cc14 =cell(1,G); cc13 = cell(1,G); cc12 = cell(1,G); cc11 = cell(1,G); cc10 = cell(1,G); cc9 = cell(1,G); 
for z=1:G
    logage =   age_offset + log2( ND(:,z) );    
    cc14{z} = logage >14;
    cc13{z} = logage <14  & logage> 13;
    cc12{z}  = logage <13 & logage > 12;
    cc11{z} = logage <12 ;
end

%% Plot Fraction of missing nuclei distributions

xlab = 'fraction of missed nuclei';
F = 12; % FontSize; 
ymax = .02; pts = 1000;

%~~~~ Plot Fraction of missing nuclei distributions  ~~~~
plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = miss_rate1{k}(cc14{k}); end
figure(31); clf;  subplot(3,1,1);
colordef white; set(gcf,'color','w');
x = linspace(0,1,10);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .05;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,.5*ymax]); title('Anterior Early cc14');


%~~~~ Plot Fraction of missing nuclei distributions  ~~~~
plot_miss = cell(1,G); 
for k=1:G;     plot_miss{k} = miss_rate1{k}(cc13{k}); end
figure(31); subplot(3,1,2);
colordef white; set(gcf,'color','w');
x = linspace(0,1,20);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,ymax]); title('Anterior cc13');

%~~~~ Plot Fraction of missing nuclei distributions ~~~~~
plot_miss = cell(1,G); 
for k=1:G;     plot_miss{k} = miss_rate1{k}(cc12{k}|cc11{k} ); end
 figure(31); subplot(3,1,3);
colordef white; set(gcf,'color','w');
x = linspace(0,1,15);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .15;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,.5*ymax]);  title('Anterior cc11 & 12');


%% Plot Fraction of missing nuclei distributions

xlab = 'fraction of missed nuclei';
F = 12; % FontSize; 
ymax = .02; pts = 1000;

%~~~~ Plot Fraction of missing nuclei distributions  ~~~~
plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = miss_rate2{k}(cc14{k}); end
figure(32); clf;  subplot(3,1,1);
colordef white; set(gcf,'color','w');
x = linspace(0,1,10);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .05;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,.5*ymax]); title('central Early cc14');


%~~~~ Plot Fraction of missing nuclei distributions  ~~~~
plot_miss = cell(1,G); 
for k=1:G;     plot_miss{k} = miss_rate2{k}(cc13{k}); end
figure(32); subplot(3,1,2);
colordef white; set(gcf,'color','w');
x = linspace(0,1,20);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,ymax]); title('central cc13');

%~~~~ Plot Fraction of missing nuclei distributions ~~~~~
plot_miss = cell(1,G); 
for k=1:G;     plot_miss{k} = miss_rate2{k}(cc12{k}|cc11{k} ); end
 figure(32); subplot(3,1,3);
colordef white; set(gcf,'color','w');
x = linspace(0,1,15);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .15;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,.5*ymax]);  title('central cc11 & 12');





%% Plot Fraction of missing nuclei distributions

xlab = 'fraction of missed nuclei';
F = 12; % FontSize; 
ymax = .02; pts = 1000;

%~~~~ Plot Fraction of missing nuclei distributions  ~~~~
plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = miss_rate3{k}(cc14{k}); end
figure(33); clf;  subplot(3,1,1);
colordef white; set(gcf,'color','w');
x = linspace(0,1,10);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .05;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,.5*ymax]); title('boundary, early cc14');


%~~~~ Plot Fraction of missing nuclei distributions  ~~~~
plot_miss = cell(1,G); 
for k=1:G;     plot_miss{k} = miss_rate3{k}(cc13{k}); end
figure(33); subplot(3,1,2);
colordef white; set(gcf,'color','w');
x = linspace(0,1,20);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,ymax]); title('boundary cc13');

%~~~~ Plot Fraction of missing nuclei distributions ~~~~~
plot_miss = cell(1,G); 
for k=1:G;     plot_miss{k} = miss_rate3{k}(cc12{k}|cc11{k} ); end
 figure(33); subplot(3,1,3);
colordef white; set(gcf,'color','w');
x = linspace(0,1,15);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .15;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,.5*ymax]);  title('boundary cc11 & 12');




%%








%% Plot Fraction of missing nuclei distributions and combine duplicates

% Hack 
G = 9;
     % all names
names = {'2 enhancers, 22C';
         'no shadow, 22C';
         'no shadow, 30C';
         'no shadow B only, 30C'; 
         'no primary, 22C';
         '2 enh B only 22C';
         '2 enh 30C';
         'no primary 30C';
         'no shadow B only, 22C'
         };
     



xlab = 'fraction of missed nuclei';
F = 12; % FontSize; 
ymax = .02; pts = 1000;
plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = miss_rate1{k}(cc14{k}); end
plot_miss{1} = [plot_miss{1}; plot_miss{6}];  
plot_miss{2} = [plot_miss{2}; plot_miss{9}]; 
plot_miss{3} = [plot_miss{3}; plot_miss{4}]; 
plot_miss{6} = [];plot_miss{9} = [];plot_miss{4} = [];

figure(34); clf;  subplot(3,1,1);
colordef white; set(gcf,'color','w');

x = linspace(0,1,20);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .05;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,.5*ymax]); title('Anterior Early cc14');


%~~~~ Plot Fraction of missing nuclei distributions  ~~~~

plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = miss_rate1{k}(cc13{k}); end
plot_miss{1} = [plot_miss{1}; plot_miss{6}];  
plot_miss{2} = [plot_miss{2}; plot_miss{9}]; 
plot_miss{3} = [plot_miss{3}; plot_miss{4}]; 
plot_miss{6} = [];plot_miss{9} = [];plot_miss{4} = [];
figure(34); subplot(3,1,2);
colordef white; set(gcf,'color','w');
x = linspace(0,1,20);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,ymax]); title('Anterior cc13');



%~~~~ Plot Fraction of missing nuclei distributions ~~~~~

plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = miss_rate1{k}(cc12{k}|cc11{k} ); end
plot_miss{1} = [plot_miss{1}; plot_miss{6}];  
plot_miss{2} = [plot_miss{2}; plot_miss{9}]; 
plot_miss{3} = [plot_miss{3}; plot_miss{4}]; 
plot_miss{6} = [];plot_miss{9} = [];plot_miss{4} = [];
 figure(34); subplot(3,1,3);
colordef white; set(gcf,'color','w');
x = linspace(0,1,20);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .15;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,.5*ymax]);  title('Anterior cc11 & 12');



%%

xlab = 'fraction of missed nuclei';
F = 12; % FontSize; 
ymax = .02; pts = 1000;
plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = miss_rate3{k}(cc14{k}); end
plot_miss{1} = [plot_miss{1}; plot_miss{6}];  
plot_miss{2} = [plot_miss{2}; plot_miss{9}]; 
plot_miss{3} = [plot_miss{3}; plot_miss{4}]; 
plot_miss{6} = [];plot_miss{9} = [];plot_miss{4} = [];

figure(35); clf;  subplot(3,1,1);
colordef white; set(gcf,'color','w');

x = linspace(0,1,20);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .05;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,.5*ymax]); title('Boundary Early cc14');

bdata14 = plot_miss;

%~~~~ Plot Fraction of missing nuclei distributions  ~~~~

plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = miss_rate3{k}(cc13{k}); end
plot_miss{1} = [plot_miss{1}; plot_miss{6}];  
plot_miss{2} = [plot_miss{2}; plot_miss{9}]; 
plot_miss{3} = [plot_miss{3}; plot_miss{4}]; 
plot_miss{6} = [];plot_miss{9} = [];plot_miss{4} = [];
figure(35); subplot(3,1,2);
colordef white; set(gcf,'color','w');
x = linspace(0,1,20);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,ymax]); title('Boundary cc13');

p2enh_v_noS = ranksum(plot_miss{1},plot_miss{2})

p_noshadow_22v30 = ranksum(plot_miss{3},plot_miss{2})
p_noprimary_22v30 = ranksum(plot_miss{5},plot_miss{8})
p_2enhancer_22v30 = ranksum(plot_miss{1},plot_miss{7})

p_22_shadow_v_primary = ranksum(plot_miss{2},plot_miss{5})
p_30_shadow_v_primary = ranksum(plot_miss{3},plot_miss{8})

bdata13 = plot_miss;

%~~~~ Plot Fraction of missing nuclei distributions ~~~~~

plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = miss_rate3{k}(cc12{k}|cc11{k} ); end
plot_miss{1} = [plot_miss{1}; plot_miss{6}];  
plot_miss{2} = [plot_miss{2}; plot_miss{9}]; 
plot_miss{3} = [plot_miss{3}; plot_miss{4}]; 
plot_miss{6} = [];plot_miss{9} = [];plot_miss{4} = [];
 figure(35); subplot(3,1,3);
colordef white; set(gcf,'color','w');
x = linspace(0,1,20);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,.5*ymax]);  title('Boundary cc11 & 12');

bdata12 = plot_miss;

%%

xlab = 'fraction of missed nuclei';
F = 12; % FontSize; 
ymax = .02; pts = 1000;
plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = miss_rate2{k}(cc14{k}); end
plot_miss{1} = [plot_miss{1}; plot_miss{6}];  
plot_miss{2} = [plot_miss{2}; plot_miss{9}]; 
plot_miss{3} = [plot_miss{3}; plot_miss{4}]; 
plot_miss{6} = [];plot_miss{9} = [];plot_miss{4} = [];

figure(36); clf;  subplot(3,1,1);
colordef white; set(gcf,'color','w');
x = linspace(0,1,20);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,.5*ymax]); title('central Early cc14');


p_noshadow22_cvb = ranksum(bdata14{2},plot_miss{2})

p_noshadow30_cvb = ranksum(bdata14{3},plot_miss{3})

p_noprimary22_cvb = ranksum(bdata14{5},plot_miss{5})

p_noprimary30_cvb = ranksum(bdata14{8},plot_miss{8})


%~~~~ Plot Fraction of missing nuclei distributions  ~~~~

plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = miss_rate2{k}(cc13{k}); end
plot_miss{1} = [plot_miss{1}; plot_miss{6}];  
plot_miss{2} = [plot_miss{2}; plot_miss{9}]; 
plot_miss{3} = [plot_miss{3}; plot_miss{4}]; 
plot_miss{6} = [];plot_miss{9} = [];plot_miss{4} = [];
figure(36); subplot(3,1,2);
colordef white; set(gcf,'color','w');
x = linspace(0,1,20);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,ymax]); title('central cc13');

p_noshadow_22v30 = ranksum(plot_miss{3},plot_miss{2})
p_noprimary_22v30 = ranksum(plot_miss{5},plot_miss{8})
p_2enhancer_22v30 = ranksum(plot_miss{1},plot_miss{7})

p_22_shadow_v_primary = ranksum(plot_miss{2},plot_miss{5})
p_30_shadow_v_primary = ranksum(plot_miss{3},plot_miss{8})


p_noshadow22_cvb = ranksum(bdata13{2},plot_miss{2})

p_noshadow30_cvb = ranksum(bdata13{3},plot_miss{3})

p_noprimary22_cvb = ranksum(bdata13{5},plot_miss{5})

p_noprimary30_cvb = ranksum(bdata13{8},plot_miss{8})

%~~~~ Plot Fraction of missing nuclei distributions ~~~~~

plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = miss_rate2{k}(cc12{k}|cc11{k} ); end
plot_miss{1} = [plot_miss{1}; plot_miss{6}];  
plot_miss{2} = [plot_miss{2}; plot_miss{9}]; 
plot_miss{3} = [plot_miss{3}; plot_miss{4}]; 
plot_miss{6} = [];plot_miss{9} = [];plot_miss{4} = [];
 figure(36); subplot(3,1,3);
colordef white; set(gcf,'color','w');
x = linspace(0,1,20);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,.5*ymax]);  title('central cc11 & 12');


p_noshadow22_cvb = ranksum(bdata12{2},plot_miss{2})

p_noshadow30_cvb = ranksum(bdata12{3},plot_miss{3})

p_noprimary22_cvb = ranksum(bdata12{5},plot_miss{5})

p_noprimary30_cvb = ranksum(bdata12{8},plot_miss{8})

