%%                          Plot_DotComps2_hbLacZ.m                      %%
%
% Analyzing New Shadow data
%
%
% Alistair Boettiger                                   Date Begun: 03/05/10
% Levine Lab                                     Functional Since: 03/06/10
%                                                   Last Modified: 03/07/11

%% Description
% comparison
%
%
%% Updates
% Changed age ID code, 10/26/10
% updated 03/08/11 to count total nuclei and to save data and use new
% plotting methods.

%% Source Code
clear all;

folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed/'; % upload data 
data_folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Code_Data/'; % export folder

emb_roots ={ 
 'hb2enh_22C_LacZ_hb';
 'hbP_22C_LacZ_hb';
% 'C33_22C_LacZ_hb';
 'hb2enh_29C_LacZ_hb';
 'C55_29C_LacZ_hb';
 'C33_29C_LacZ_hb';
 };


names = {
    '2 enhancer 22C';
    'primary 22C';
 %   'shadow 22C';
    '2 enhancer 30C';
    'primary 30C';
    'shadow 30C';
    };
    
    

N = 60;
K = length(emb_roots); 
G= length(names);

     age_table = cell(1,K);

miss_cnt = cell(1,K); 
miss_rate = cell(1,K); 
nd = cell(1,K); 
lowon = cell(1,K); 
cell_var = cell(1,K); 
ectop_cnt = cell(1,K);
ectop_rate = cell(1,K);
endog_cnt = cell(1,K);
rept_cnt = cell(1,K); 

for z=1:K
    miss_cnt{z} = zeros(N,1);
    miss_rate{z} = zeros(N,1); 
    lowon{z} = zeros(N,1); 
    nd{z} = zeros(N,1);
    age_table{z} = cell(N,2);
    ectop_cnt{z} = zeros(N,1);
    ectop_rate{z} = zeros(N,1);
    endog_cnt{z} = zeros(N,1);
    rept_cnt{z} = zeros(N,1);    
end


xmin = .2; xmax = .9; ymin = .15; ymax = .4;
% as fractions of the original image dimensions.  


for z=1:K % K % k=2;
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
          miss_cnt{z}(n) =  length(intersect(setdiff(pts2,pts1), ptr_nucin2));  
         %  miss_cnt{z}(n) = anlz_major_reg(folder,emb_roots{z},emb );
           miss_rate{z}(n) = miss_cnt{z}(n)/length(pts2); 
          %    [lowon{z}(n)] = lowon_fxn(H,handles,nin2,ptr_nucin2,[emb_roots{z},emb],1); 
          endog_cnt{z}(n) = length(pts2); 
          rept_cnt{z}(n) = length(pts1);
           %lowon{z}(n) = lowon_fxn(H,handles,all_nucs,pts2,nin2,Cell_bnd);  
          
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

save([data_folder,'hb_LacZ_data_030811']);
disp('finished');



%%
 clear all; 
 
data_folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Code_Data/'; % export folder
load([data_folder,'hb_LacZ_data_030811']);
ND = cell2mat(nd); 

age_offset = 4.8;
emb_cycle = age_offset + log2( nonzeros( sort(ND(:)) ) );
figure(10); clf; plot( emb_cycle ,'r.');

T_embs = length(nonzeros(ND(:))) ;
title(['hb embryos, N = ',num2str(T_embs)  ],'FontSize',15);
set(gca,'FontSize',15); grid on;
set(gcf,'color','w'); ylabel('log_2(nuc density)');  xlabel('embryo number'); 
ylim([10,14.99]); xlim([0,T_embs + 10]);


%%
cc14 =cell(1,G); cc13 = cell(1,G); cc12 = cell(1,G); cc11 = cell(1,G); cc10 = cell(1,G); cc9 = cell(1,G); 
for z=1:G
    logage =   age_offset + log2( ND(:,z) );
    
    cc14{z} = logage >14;
    cc13{z} = logage <14  & logage> 13;
    cc12{z}  = logage <13 & logage > 12;
    cc11{z} = logage <12 ;
end


%% Compare endog vs rept
F = 12;
xlab = 'expressing fraction of nuclei';

names = {'hb 2 enhancers, hb 22C';      % 1
         'hb P, hb 22C';                % 2
         %'hb S, hb 22C';
         'hb 2 enhancers, hb 30C';      % 3
         'hb P, hb 30C';                % 4
         'hb S, hb 30C';                % 5
         'hb 2 enhancers, LacZ 22C';    % 6
         'hb P, LacZ 22C';              % 7
         %'hb S, LacZ 22C';
         'hb 2 enhancers, LacZ 30C';    % 8
         'hb P, LacZ 30C';              % 9
         'hb S, LacZ 30C';              % 10
         };

colordef white; 
endog = cell(1,G); 
rept = cell(1,G);
 for k=1:G;     endog{k} = endog_cnt{k}(cc13{k})./ND(cc13{k},k)/5; end
 for k=1:G;     rept{k} = rept_cnt{k}(cc13{k})./ND(cc13{k},k)/5; end

  data = cat(2,endog,rept);    
  Ts = length(data);% number of tracks
  pW = zeros(Ts);
  pA = zeros(Ts); 
  for i=1:Ts
    for j = 1:Ts
     pW(i,j) = ranksum(data{i},data{j});   % Wilcox Rank Sum
     pA(i,j)=anovan([data{i}',data{j}'],{[zeros(1,length(data{i})),ones(1,length(data{j}))]},'display','off'); % 2-way ANOVA
    end
  end
 Wpvals = ['p_{16} = ',num2str(pW(1,6),2), '   p_{27} = ',num2str(pW(2,7),2) ];
  disp(['22C pairwise Wilcoxon rank sum:  ', Wpvals]);
 
 Wpvals = ['p_{38} = ',num2str(pW(3,8),2), '   p_{49} = ',num2str(pW(4,9),2) , '   p_{27} = ',num2str(pW(5,10),2) ];
 Apvals = ['p_{38} = ',num2str(pA(3,8),2), '   p_{49} = ',num2str(pA(4,9),2) , '    p_{38} = ',num2str(pA(5,10),2)  ];
 disp(['30C pairwise Wilcoxon rank sum:  ', Wpvals]);
 disp(['2-way ANOVA:  ',Apvals]);
%  
%  figure(1); clf;
%  cityscape(data,names,xlab,F);
 
 figure(3); clf;
  cumhist(data,names,xlab,F);
  title(['pairwise Wilcoxon:  ' Wpvals]);
  set(gcf,'color','w');

%%
F = 12;
xlab = 'fraction of missed nuclei';
names = {
    '2 enhancer 22C';
    'primary 22C';
  %  'shadow 22C';
    '2 enhancer 30C';
    'primary 30C';
    'shadow 30C';
    };
    

plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = miss_rate{k}(cc13{k}); end


  data = plot_miss;    
  
  
  Ts = length(data);% number of tracks
  pW = zeros(Ts);
  pA = zeros(Ts); 
  for i=1:Ts
    for j = 1:Ts
     pW(i,j) = ranksum(data{i},data{j});   % Wilcox Rank Sum
     pA(i,j)=anovan([data{i}',data{j}'],{[zeros(1,length(data{i})),ones(1,length(data{j}))]},'display','off'); % 2-way ANOVA
    end
  end
 Wpvals = ['p_{12} = ',num2str(pW(3,4),2), '   p_{13} = ',num2str(pW(3,5),2) , '    p_{23} = ',num2str(pW(4,5),2)  ];
 Apvals = ['p_{34} = ',num2str(pA(3,4),2), '   p_{35} = ',num2str(pA(3,5),2) , '    p_{23} = ',num2str(pA(4,5),2)  ];
 disp(['pairwise Wilcoxon rank sum:  ', Wpvals]);
 disp(['2-way ANOVA:  ',Apvals]);
 
%  figure(1); clf;
%  cityscape(data,names,xlab,F);
 
 figure(3); clf;
  cumhist(data([2:5]),names([2:5]),xlab,F);
  title(['pairwise Wilcoxon:  ' Wpvals]);
  set(gcf,'color','w');


    
disp([names{3}, ': ' ,num2str(median([data{3}])),'+/-',num2str(std([data{3}])),  ' missing']);
disp([names{4}, ': ' ,num2str(median([data{4}; data{2}])),'+/-',num2str(std([data{4}; data{2}])),  ' missing']);
disp([names{5}, ': ' ,num2str(median([data{5}])),'+/-',num2str(std([data{5}])),  ' missing']);

%%  Ectopic expression rate

xlab = 'ectopic expression rate';

plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = ectop_rate{k}(cc14{k}); end
  
 data = plot_miss;
  Ts = length(data);% number of tracks
  pW = zeros(Ts);
  pA = zeros(Ts); 
  for i=1:Ts
    for j = 1:Ts
     pW(i,j) = ranksum(data{i},data{j});   % Wilcox Rank Sum
     pA(i,j)=anovan([data{i}',data{j}'],{[zeros(1,length(data{i})),ones(1,length(data{j}))]},'display','off'); % 2-way ANOVA
    end
  end
 Wpvals = ['p_{12} = ',num2str(pW(1,2),2), '   p_{13} = ',num2str(pW(1,3),2) , '    p_{23} = ',num2str(pW(2,3),2)  ];
 Apvals = ['p_{12} = ',num2str(pA(1,2),2), '   p_{13} = ',num2str(pA(1,3),2) , '    p_{23} = ',num2str(pA(2,3),2)  ];
 disp(['pairwise Wilcoxon rank sum:  ', Wpvals]);
 disp(['2-way ANOVA:  ',Apvals]);
 
 figure(2); clf;
 cityscape(data,names,xlab,F);
 
 figure(4); clf;
  cumhist(data,names,xlab,F);
  title(['pairwise Wilcoxon:  ' Wpvals]);
  set(gcf,'color','w');



























% 
% 
% 
% 
% 
% ND = cell2mat(nd); 
% age_offset = 4.8;
% 
% emb_cycle = age_offset + log2( nonzeros( sort(ND(:)) ) );
% figure(2); clf; plot( emb_cycle ,'r.');
% 
% 
% title(['hb embryos, N = ',num2str(length(nonzeros(ND(:))) )  ],'FontSize',15);
% set(gca,'FontSize',15); grid on;
% set(gcf,'color','w'); ylabel('log_2(nuc density)');  xlabel('embryo number'); 
% ylim([10,14.99]);
% 
% 
% %%
% G= length(names);
% cc14 =cell(1,G); cc13 = cell(1,G); cc12 = cell(1,G); cc11 = cell(1,G); cc10 = cell(1,G); cc9 = cell(1,G); 
% for z=1:G
%     logage =   age_offset + log2( ND(:,z) );
%     
%     cc14{z} = logage >14;
%     cc13{z} = logage <14  & logage> 13;
%     cc12{z}  = logage <13 & logage > 12;
%     cc11{z} = logage <12 ;
% end
% 
% %% Plot Fraction of missing nuclei distributions
% 
% xlab = 'fraction of missed nuclei';
% F = 12; % FontSize; 
% ymax = .02; pts = 1000;
% 
% plot_miss = cell(1,G); 
%  for k=1:G;     plot_miss{k} = miss_rate{k}(cc14{k}); end
% 
% 
% figure(33); clf;  subplot(3,1,1);
% colordef white; set(gcf,'color','w');
% 
% x = linspace(0,1,8);  % range and number of bins for histogram
% xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
%  method = 'pcubic'; % method for interpolation
% sigma = .05;  % smoothing factor for interpolation
% CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
% ylim([0,.5*ymax]); title('Early cc14');
% 
% 
% %~~~~ Plot Fraction of missing nuclei distributions  ~~~~
% 
% plot_miss = cell(1,G); 
%  for k=1:G;     plot_miss{k} = miss_rate{k}(cc13{k}); end
% 
% figure(33); subplot(3,1,2);
% colordef white; set(gcf,'color','w');
% 
% x = linspace(0,1,15);  % range and number of bins for histogram
% xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
%  method = 'pcubic'; % method for interpolation
% sigma = .1;  % smoothing factor for interpolation
% CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
% ylim([0,ymax]); title('cc13');
% 
% 
% 
% %~~~~ Plot Fraction of missing nuclei distributions ~~~~~
% 
% plot_miss = cell(1,G); 
%  for k=1:G;     plot_miss{k} = miss_rate{k}(cc12{k}|cc11{k} ); end
% 
% 
%  figure(33); subplot(3,1,3);
% colordef white; set(gcf,'color','w');
% 
% x = linspace(0,1,15);  % range and number of bins for histogram
% xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
%  method = 'pcubic'; % method for interpolation
% sigma = .15;  % smoothing factor for interpolation
% CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
% ylim([0,.5*ymax]);  title('cc11 & 12');
% 
% 
% %%
% 
% 
% 
% 
% %% Total Expression Variability
% 
% %%
% xlab = '\sigma/\mu';
% 
% 
% plot_miss = cell(1,G); 
%  for k=1:G;     plot_miss{k} = lowon{k}(cc14{k}); end
% % for k=1:G;     plot_miss{k} = miss_rate{k}; end
% 
% 
%  figure(1); clf;
% % colordef black; set(gcf,'color','k');
% colordef white; set(gcf,'color','w');
% 
% x = linspace(0,1,14);  % range and number of bins for histogram
% xx = linspace(0,1,1000); % range a number of bins for interpolated distribution
%  method = 'pcubic'; % method for interpolation
% sigma = .1;  % smoothing factor for interpolation
% y = CompDist(plot_miss,x,xx,method,sigma,names,xlab,14);
% 
% title('cc13 embryos');
% 
% 
