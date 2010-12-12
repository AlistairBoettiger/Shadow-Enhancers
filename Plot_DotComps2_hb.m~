%%                          Plot_DotComps2_hb.m                             %%
%
% Analyzing New Shadow data
%
%
% Alistair Boettiger                                   Date Begun: 03/05/10
% Levine Lab                                     Functional Since: 03/06/10
%                                                   Last Modified: 12/09/10

%% Description
% comparison
%
%
%% Updates
% Modified 10/18/10 to incldue repression of ectopic expression. 
% Modified to hack failed combine data code
% Modified 12/09/10 to add new datasets


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
            'BAC02_22C_y_hb'; % 9
            'BAC01b_22C_y_hb';
            'BAC01b_30C_y_hb';
            'BAC02b_30C_y_hb';
            };  
          %1 +6, 2+9, 3+4, 
     

%      % all names
% names = {'2 enhancers, 22C';
%          'no shadow, 22C';
%          'no shadow, 30C';
%          'no shadow, 30C, b'; 
%          'no primary, 22C';
%          '2 enh 22C';
%          '2 enh 30C';
%          'no primary 30C';
%          'no shadow, 22C, b';
%          'no primary 22C, b';
%          'no primary 30C, b';
%          'no shadow 30C, c'
%          };

     

N = 100;
K = length(emb_roots); 
% G= length(names);


     age_table = cell(1,K);
miss_cnt = cell(1,K); 
miss_rate = cell(1,K); 
nd = cell(1,K); 
lowon = cell(1,K); 
ectop_cnt = cell(1,K);
ectop_rate = cell(1,K);

for z=1:K
    miss_cnt{z} = zeros(N,1);
    miss_rate{z} = zeros(N,1); 
    lowon{z} = zeros(N,1); 
    nd{z} = zeros(N,1);
    age_table{z} = cell(N,2);
    ectop_cnt{z} = zeros(N,1);
    ectop_rate{z} = zeros(N,1);
end

xmin = .2; xmax = .9; ymin = .15; ymax = .4;
% as fractions of the original image dimensions.  

for z=1:K % k=2;
    for n=   1:N
        if n<10
            emb = ['0',num2str(n)];
        else
            emb = num2str(n);
        end

        try
        load([folder,'/',emb_roots{z},emb,'_data.mat']);   
        % get the indices of all nuclei in green that are not also red.  
        % require these nuclei also fall in the 'region' for red nuclei.  
       % s29_miss_cnt(n) =  length(intersect(setdiff(pts2,pts1), ptr_nucin2));
          % miss_cnt{z}(n) =  length(intersect(setdiff(pts2,pts1), ptr_nucin2));  
             [miss_rate{z}(n), ptr_nucin2] = anlz_major_reg2(folder,emb_roots{z},emb );
           % miss_rate{z}(n) = miss_cnt{z}(n)/length(pts2); 
           % lowon{z}(n) = lowon_fxn(H,handles,nin2,ptr_nucin2,emb); 
          
           
           ectop_cnt{z}(n) = length(intersect(setdiff(pts1,pts2),setdiff(ptr_nucin1,ptr_nucin2)')); 
           if length(H) > 2000
               im_dim = 2048;
           else
               im_dim = 1024;
           end
               
           lims =  round([xmin,xmax,ymin,ymax]*im_dim); 
           nd{z}(n) = NucDensity(cent,lims,0);
            age_table{z}{n,1} = [folder,'/',emb_roots{z},emb,'_data.mat']; %  
            age_table{z}{n,2} = nd{z}(n); 
            ectop_rate{z}(n) = ectop_cnt{z}(n)/nd{z}(n);
        catch ME
            disp(ME.message); 
            %disp(['can not find file' folder,'/',emb_roots{z},emb,'_data.mat']);
        end
    end
end

close all; 

save hb_SD-12-09-10
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
  clear all; load hb_SD-12-09-10

 % Merge data
 
%  %
%       % all names
% names = {'2 enhancers, 22C';  % 1
%          'no shadow, 22C';    % 2
%          'no shadow, 30C';    % 3
%          'no shadow, 30C, b'; % 4
%          'no primary, 22C';   % 5
%          '2 enh 22C';         % 6
%          '2 enh 30C';         % 7
%          'no primary 30C';    % 8
%          'no shadow, 22C, b'; % 9
%          'no primary 22C, b'; % 10 
%          'no primary 30C, b'; % 11
%          'no shadow 30C, c'   % 12
%          };


 Nnuc = cell(1,6); foff = cell(1,6);  

 foff{1} = [miss_rate{1}; miss_rate{6}] ;  Nnuc{1} = [nd{1}; nd{6}] ;  % 2 enhancer 22C
 
 foff{2} = miss_rate{7};  Nnuc{2} = nd{7}; % 2 enhancer 30C
 
 foff{3} = [miss_rate{2}; miss_rate{9}];  Nnuc{3} = [nd{2}; nd{9}]; % no shadow 22C
  
 foff{4} = [miss_rate{3}; miss_rate{4}; miss_rate{12}];  Nnuc{4} = [nd{3}; nd{4}; nd{12}]; % no shadow 30C
 
 foff{5} = [miss_rate{5}; miss_rate{10}];  Nnuc{5} = [nd{5}; nd{10}]; % no primary 22C

 foff{6} = [miss_rate{8}; miss_rate{11}];  Nnuc{6} = [nd{8}; nd{11}]; % no primary 30C
 
 G = length(foff);
 
 for k=1:G
    data = nonzeros(foff{k});
    foff{k} = [data; zeros(200-length(data),1)];
    data = nonzeros(Nnuc{k});
    Nnuc{k} = [data; zeros(200-length(data),1)];
 end
     
 names = {'control 22C';
          'control 30C';
          'no shadow 22C';
          'no shadow 30C';
          'no primary 22C';
          'no primary 30C'
          };
      
 
 %%
 
 
%ND = cell2mat(nd); 
ND = cell2mat(Nnuc); 
age_offset = 4.8;

emb_cycle = age_offset + log2( nonzeros( sort(ND(:)) ) );
figure(2); clf; plot( emb_cycle ,'r.');


title(['hb embryos, N = ',num2str(length(nonzeros(ND(:))) )  ],'FontSize',15);
set(gca,'FontSize',15); grid on;
set(gcf,'color','w'); ylabel('log_2(nuc density)');  xlabel('embryo number'); 
ylim([10,14.99]);


%%
G= length(names);
cc14 =cell(1,G); cc13 = cell(1,G); cc12 = cell(1,G); cc11 = cell(1,G); cc10 = cell(1,G); cc9 = cell(1,G); 
for z=1:G
    logage =   age_offset + log2( ND(:,z) );
    
    cc14{z} = logage >14;
    cc13{z} = logage <14  & logage> 13;
    cc12{z}  = logage <13 & logage > 12;
    cc11{z} = logage <12 & logage > 0 ;
    foff{z}(foff{z}==Inf) = 0; 
end

%% Plot Fraction of missing nuclei distributions

xlab = 'fraction of missed nuclei';
F = 12; % FontSize; 
ymax = .02; pts = 1000;
xmax = 1;

plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = foff{k}(cc14{k}); end


figure(33); clf;  subplot(3,1,1);
colordef white; set(gcf,'color','w');

x = linspace(0,1,15);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .01;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,.5*ymax]); title('Early cc14');
xlim([0,xmax]);

figure(30); clf;  subplot(3,1,1); 
BoxDist(plot_miss,names,'fraction missed');
xlim([0,1]);


%~~~~ Plot Fraction of missing nuclei distributions  ~~~~

plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = foff{k}(cc13{k}); end

figure(33); subplot(3,1,2);
colordef white; set(gcf,'color','w');

x = linspace(0,1,40);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .05;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,ymax]); title('cc13');
xlim([0,xmax]);


figure(30); subplot(3,1,2); 
BoxDist(plot_miss,names,'fraction missed');
xlim([0,1]);

%~~~~ Plot Fraction of missing nuclei distributions ~~~~~

plot_miss = cell(1,G); % k = 10;
 for k=1:G;     plot_miss{k} = foff{k}(cc12{k}|cc11{k} ); end


 figure(33); subplot(3,1,3); % figure(2); clf; 
colordef white; set(gcf,'color','w');

x = linspace(0,1,12);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .05;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,.5*ymax]);  title('cc11 & 12');
xlim([0,xmax]);

%


figure(30); subplot(3,1,3); 
BoxDist(plot_miss,names,'fraction missed');
xlim([0,1]);



%%  Ectopic expression rate

plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = ectop_rate{k}(cc14{k}); end
% for k=1:G;     plot_miss{k} = miss_rate{k}; end
pts = 50; 

 figure(5); clf;
% colordef black; set(gcf,'color','k');
colordef white; set(gcf,'color','w');

x = linspace(0,1,14);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,14);

title('cc14 embryos');

xlab = 'fraction of ectopic on nuclei';


plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = ectop_rate{k}(cc13{k}); end

 figure(6); clf;
% colordef black; set(gcf,'color','k');
colordef white; set(gcf,'color','w');

x = linspace(0,1,15);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation
y = CompDist(plot_miss,x,xx,method,sigma,names,xlab,14);
title('cc13 embryos'); 

xlab = 'fraction of ectopic on nuclei';


plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = ectop_rate{k}(cc12{k}|cc11{k}); end
% for k=1:G;     plot_miss{k} = miss_rate{k}; end


 figure(7); clf;
% colordef black; set(gcf,'color','k');
colordef white; set(gcf,'color','w');

x = linspace(0,1,25);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,14);

title('cc12 and 11 embryos');

xlab = 'fraction of ectopic on nuclei';




