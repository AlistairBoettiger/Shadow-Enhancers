%%                          Plot_DotComps2_sna.m                         %%
%
% Analyzing New Shadow data
%
%
% Alistair Boettiger                                   Date Begun: 03/05/10
% Levine Lab                                     Functional Since: 03/06/10
%                                                   Last Modified: 01/12/11

%% Description
% comparison
%
%
%% Updates
% Revised 01/12/11 to use most recent  formulation of age structure and
% plotting tools.  
% Really should make this into a function

%% Source Code
clear all;

folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed';

emb_roots = {'MP05_29C_y_sna'; % 1
    'MP05_22C_y_sna';          % 2
    'MP06xYW_30C_y_sna';       % 3
    'MP06xYW_22C_y_sna';       % 4
    'MP10_29C_y_sna';          % 5
    'MP10_22C_y_sna';          % 6
    'MP05xYW_30C_sna_y-full';  % 7
    'MP10xYW_30C_sna_y-full';  % 8
    'MP05xdl6_25C_pt1';        % 9
    'MP05xdl6_25C_pt2';        % 10
    'MP10xdl6_25C_pt1';        % 11 
    'MP10xdl6_25C_pt2'};       % 12
          

 names = {'2 enh 22C';
          '2 enh 30C';
          'no primary 22C';
          'no primary 30C';
          'no shadow 22C';
          'no shadow 30C';
          '2 enh dl6';
          'no primary dl6'};
 
N = 70;
K = length(emb_roots); 
G= length(names);

miss_cnt = cell(1,K); 
miss_rate = cell(1,K); 
nd = cell(1,K); 
lowon = cell(1,K); 

for z=1:K
    miss_cnt{z} = zeros(N,1);
    miss_rate{z} = zeros(N,1); 
    lowon{z} = zeros(N,1); 
    nd{z} = zeros(N,1);
end


xmin = .2; xmax = .8; ymin = .25; ymax = .75;
% as fractions of the original image dimensions.  

for z=1:K % k=2;
    for n=1:N
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
           miss_rate{z}(n) = miss_cnt{z}(n)/length(pts2); 
           lowon{z}(n) = lowon_fxn(H,handles,nin2,ptr_nucin2,emb,0); 
           %lowon{z}(n) = lowon_fxn(H,handles,all_nucs,pts2,nin2,Cell_bnd);  
          
           if length(H) > 2000
               im_dim = 2048;
           else
               im_dim = 1024;
           end
               
           lims =  round([xmin,xmax,ymin,ymax]*im_dim); 
           nd{z}(n) = NucDensity(cent,lims,0);

        catch ME
            disp(ME.message); 
            %disp(['can not find file' folder,'/',emb_roots{z},emb,'_data.mat']);
        end
    end
end

 % save snail_SD_011211





%% 
% clear all; 
load('/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Code_Data/snail_SD_011211.mat');

foff{1} = miss_rate{6};                 Nnuc{1} = nd{6}; % 2 enh 22 C, MP10
foff{2} = [miss_rate{5},miss_rate{8}];  Nnuc{2} = [nd{5},nd{8}]; % 2 enh 30C MP10
foff{3} = miss_rate{2};                 Nnuc{3} = nd{2}; % MP05 22C 
foff{4} = [miss_rate{1}; miss_rate{7}]; Nnuc{4} = [nd{1}, nd{7}]; % MP05 30C
foff{5} = miss_rate{4};                 Nnuc{5} = nd{4};  % MP06 22C
foff{6} = miss_rate{3};                 Nnuc{6} = nd{3};  % MP06 30C
foff{7} = [miss_rate{11},miss_rate{12}]; Nnuc{7} = [nd{11}, nd{12}];% MP10 dl6 25C
foff{8} = [miss_rate{9},miss_rate{10}]; Nnuc{8} = [nd{9},nd{10}]; % MP05 dl6 25C

 
 for k=1:G
    data = nonzeros(foff{k});
    foff{k} = [data; zeros(200-length(data),1)];
    data = nonzeros(Nnuc{k});
    Nnuc{k} = [data; zeros(200-length(data),1)];
 end

 names = {'2 enh 22C';
          '2 enh 30C';
          'no primary 22C';
          'no primary 30C';
          'no shadow 22C';
          'no shadow 30C';
          '2 enh dl6';
          'no primary dl6'};
 
% 
% 
% emb_roots = {'MP05_29C_y_sna'; % 1
%     'MP05_22C_y_sna';          % 2
%     'MP06xYW_30C_y_sna';       % 3
%     'MP06xYW_22C_y_sna';       % 4
%     'MP10_29C_y_sna';          % 5
%     'MP10_22C_y_sna';          % 6
%     'MP05xYW_30C_sna_y-full';  % 7
%     'MP10xYW_30C_sna_y-full';  % 8
%     'MP05xdl6_25C_pt1';        % 9
%     'MP05xdl6_25C_pt2';        % 10
%     'MP10xdl6_25C_pt1';        % 11 
%     'MP10xdl6_25C_pt2'};       % 12

%%

ND = cell2mat(Nnuc); 
age_offset = 5.3;

emb_cycle = age_offset + log2( nonzeros( sort(ND(:)) ) );
figure(2); clf; plot( emb_cycle ,'r.');


title(['sna embryos, N = ',num2str(length(nonzeros(ND(:))) )  ],'FontSize',15);
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

 names = {'2 enh 22C';
          '2 enh 30C';
          'no primary 22C';
          'no primary 30C';
          'no shadow 22C';
          'no shadow 30C';
          '2 enh dl6';
          'no primary dl6'};


plot_miss = cell(1,G); 
for k=1:G;     plot_miss{k} = foff{k}(cc14{k}); end
%for k=1:G;     plot_miss{k} = miss_rate{k}; end


 figure(1); clf;
 colordef black; set(gcf,'color','k');
%colordef white; set(gcf,'color','w');

x = linspace(0,1,30);  % range and number of bins for histogram
xx = linspace(0,1,100); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation
subplot(2,1,1);
CompDist(plot_miss([1,3,5]),x,xx,method,sigma,names([1,3,5]),xlab,12)

subplot(2,1,2);
sigma = .1; x = linspace(0,1,20);
CompDist(plot_miss([2,4,6]),x,xx,method,sigma,names([2,4,6]),xlab,12)



labs = {'30C','22C'};
figure(30); clf;  
BoxDist(plot_miss,names,'fraction missed',labs);
xlim([0,1]);

%%



%%  expression
F = 14;
xlab = 'missed expression';

plot_miss = cell(1,G); 
for k=1:G;     plot_miss{k} = foff{k}(cc14{k}); end
  
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

  Wpvals = ['   p_{24} = ',num2str(pW(2,4),2) , '    p_{26} = ',num2str(pW(2,6),2) , '    p_{78} = ',num2str(pW(7,8),2) ];
 figure(4); clf;
  cumhist(data,names,xlab,F);
  title(['pairwise Wilcoxon:  ' Wpvals]);
  set(gcf,'color','w');

disp([names{1},': ' ,num2str(median([data{1}])),'+/-',num2str(std([data{1}])),  ' missed']);
disp([names{2},': ' ,num2str(median([data{2}])),'+/-',num2str(std([data{2}])),  ' missed']);
disp([names{3},': ' ,num2str(median([data{3}])),'+/-',num2str(std([data{3}])),  ' missed']);




 figure(3); clf;
  cumhist(data([1:6]),names([1:6]),xlab,F);
  title(['pairwise Wilcoxon:  ' ['   p_{24} = ',num2str(pW(2,4),2) , '    p_{26} = ',num2str(pW(2,6),2)  ];]);
  set(gcf,'color','w');

 figure(5); clf;
  cumhist(data([1,3,5]),names([1,3,5]),xlab,F);
  title(['pairwise Wilcoxon:  ' ['   p_{12} = ',num2str(pW(1,3),2) , '    p_{13} = ',num2str(pW(1,5),2) , '    p_{23} = ',num2str(pW(3,5),2) ];]);
  set(gcf,'color','w');


   figure(6); clf;
  cumhist(data([2,4,6]),names([2,4,6]),xlab,F);
  title(['pairwise Wilcoxon:  ' ['   p_{12} = ',num2str(pW(2,4),2) , '    p_{13} = ',num2str(pW(2,6),2) , '    p_{23} = ',num2str(pW(4,6),2) ];]);
  set(gcf,'color','w');
  
   figure(7); clf;
  cumhist(data([1,3,7,8]),names([1,3,7,8]),xlab,F);
  title(['pairwise Wilcoxon:  ' ['   p_{13} = ',num2str(pW(2,7),2) , '    p_{24} = ',num2str(pW(3,8),2)  ];]);
  set(gcf,'color','w');

%%  Compare to bionmial

N = 700; % estimate of number of cells

for k=1:G
    mu(k) = mean(plot_miss{k});
    sig(k) = std(plot_miss{k});
    bisig(k) = sqrt( mu(k)*N*(1-mu(k)) )/N;
end

figure(3); clf;
scatter(sig,bisig);



% %% Plot Fraction of missing nuclei distributions
% 
% xlab = 'fraction of missed nuclei';
% 
% 
% plot_miss = cell(1,G); 
% for k=1:G;     plot_miss{k} = miss_rate{k}(cc14{k}); end
% 
%  figure(1); clf;
% % colordef black; set(gcf,'color','k');
% colordef white; set(gcf,'color','w');
% 
% % x = linspace(0,1,8);  % range and number of bins for histogram
% % xx = linspace(0,1,100); % range a number of bins for interpolated distribution
% %  method = 'pcubic'; % method for interpolation
% % sigma = .1;  % smoothing factor for interpolation
% % CompDist(plot_miss,x,xx,method,sigma,names,xlab)
% 
% BoxDist(plot_miss,names,xlab);
% set(gcf,'color','k');
% 
% V= plot_miss;
%     P_var = zeros(G,G); 
%     for i=1:G
%         for j=1:G   
%    P_var(i,j) =  log10(ranksum(V{i},V{j}));
%         end
%     end
%     
%     figure(6); clf; imagesc(P_var); colorbar;  colormap('gray');
%     set(gca,'YtickLabel', str2mat(names{:}),'YTick',1:6,'fontsize',15,...
%     'YMinorTick','on'); title(xlab);
%     set(gcf,'color','k');

% %% Plot Total mRNA variability distribuitons 
% 
% xlab = 'variability in total transcript (\sigma/\mu)';
% 
% 
% plot_lowon = cell(1,G); 
% for k=1:G; plot_lowon{k} = lowon{k}(cc14{k}); end
% figure(2); clf;
%  colordef black; set(gcf,'color','k');
% %colordef white; set(gcf,'color','w');
% 
% 
% %  x = linspace(0,1,20);
% %   xx = linspace(0,1,100); 
% %  method = 'pcubic';
% % sigma = .1; 
% % CompDist(plot_lowon,x,xx,method,sigma,names,xlab)
% 
% BoxDist(plot_lowon,names,xlab);
% set(gcf,'color','k');
% 
% V= plot_lowon;
%     P_var = zeros(G,G); 
%     for i=1:G
%         for j=1:G   
%    P_var(i,j) =  log10(ranksum(V{i},V{j}));
%         end
%     end
%     
%     figure(6); clf; imagesc(P_var); colorbar;  colormap('gray');
%     set(gca,'YtickLabel', str2mat(names{:}),'YTick',1:6,'fontsize',15,...
%     'YMinorTick','on'); title(xlab);
%     set(gcf,'color','k');
% 
% %% Variability between immidiate neighbors 
%  
% xlab = 'variability in total transcript among neighbors \sigma/mu';
% % x = linspace(0,1,33);
% %   xx = linspace(0,1,100); 
% %  method = 'pcubic';
% % sigma = .1; 
% 
% plot_cell_var = cell(1,G); 
% for k=1:G;  plot_cell_var{k} = cell_var{k}(cc14{k}); end
% plot_cell_var{4} = [plot_cell_var{4}; cell_var{7}(cc14{7})];
% 
%  figure(3); clf;
% % colordef black; set(gcf,'color','k');
% colordef white; set(gcf,'color','w');
% BoxDist(plot_cell_var,names,xlab);
% set(gcf,'color','k');
% % CompDist(plot_cell_var,x,xx,method,sigma,names,xlab);