%%                          Plot_DotComps2_Kr.m                          %%
%
% Analyzing New Shadow data
%
%
% Alistair Boettiger                                   Date Begun: 03/05/10
% Levine Lab                                     Functional Since: 09/13/10
%                                                   Last Modified: 03/07/11

%% Description
% comparison
%
%
%% Updates
% Modified 10/18 to also count ectopic nuclei 
% Modified 02/28/11 to use cityscape and cumulative sum.  
% Modified 03/07/11 to compare total reporter cells to total endogneous.  
%% Source Code
clear all;



folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed';

emb_roots = {'kr2enh_22C_LacZ_kr';  
             'krCD1_22C_LacZ_kr';
             'krCD2_22C_LacZ_kr'
            };  
          

names = {'Kr 2 enhancers, 22C';  % 
         'Kr distal, 22C';  % CD1
         'Kr proximal, 22C'  % CD2
         };


N = 40;
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
          miss_cnt{z}(n) =  length(setdiff(pts2,pts1));
      %     miss_cnt{z}(n) =  length(intersect(setdiff(pts2,pts1), ptr_nucin2));  
          % miss_cnt{z}(n) = anlz_major_reg(folder,emb_roots{z},emb );
          endog_cnt{z}(n) = length(pts2); 
          rept_cnt{z}(n) = length(pts1);
           miss_rate{z}(n) = miss_cnt{z}(n)/length(pts2); 
          %  [lowon{z}(n),cell_var{z}(n)] = lowon_fxn(H,handles,nin2,ptr_nucin2,[emb_roots{z},emb],0);     
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
        end
    end
end

close all; 

data_folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Code_Data/';
save([data_folder,'kr_LacZ_data_030711']);

% save kr_LacZ_data2;
%save kr_LacZ_data_ect;



%%




%%

    clear all; 
data_folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Code_Data/';
load([data_folder,'kr_LacZ_data_ect']);
%%
  
 % load kr_LacZ_data_ect;
ND = cell2mat(nd); 

emb_cycle = 4.8 + log2( nonzeros( sort(ND(:)) ) );
figure(10); clf; plot( emb_cycle ,'r.');

T_embs = length(nonzeros(ND(:))) ;
title(['kr embryos, N = ',num2str(T_embs)  ],'FontSize',15);
set(gca,'FontSize',15); grid on;
set(gcf,'color','w'); ylabel('log_2(nuc density)');  xlabel('embryo number'); 
ylim([10,14.99]); xlim([0,T_embs + 10]);


%%
cc14 =cell(1,G); cc13 = cell(1,G); cc12 = cell(1,G); cc11 = cell(1,G); cc10 = cell(1,G); cc9 = cell(1,G); 
for z=1:G
    logage =   4.8 + log2( ND(:,z) );
    
    cc14{z} = logage >14;
    cc13{z} = logage <14  & logage> 13;
    cc12{z}  = logage <13 & logage > 12;
    cc11{z} = logage <12 ;
end




%% Fraction of missing nuclei
F = 14;
xlab = 'fraction of missed nuclei';

names = {'Kr 2 enhancers, 22C';  % 
         'Kr distal, 22C';  % CD1
         'Kr proximal, 22C'  % CD2
         };

colordef white; 
plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = miss_rate{k}(cc14{k}); end


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
 Wpvals = ['p_{CA} = ',num2str(pW(1,2),2), '   p_{CB} = ',num2str(pW(1,3),2) , '    p_{AB} = ',num2str(pW(2,3),2)  ];
 Apvals = ['p_{12} = ',num2str(pA(1,2),2), '   p_{13} = ',num2str(pA(1,3),2) , '    p_{23} = ',num2str(pA(2,3),2)  ];
 disp(['pairwise Wilcoxon rank sum:  ', Wpvals]);
 disp(['2-way ANOVA:  ',Apvals]);
 
 figure(1); clf;
 cityscape(data,names,xlab,F);
 
 figure(3); clf;
  cumhist(data,names,xlab,F);
  title(['pairwise Wilcoxon:  ' Wpvals]);
  set(gcf,'color','w');


    
disp([names{1}, ': ' ,num2str(median([data{1}])),'+/-',num2str(std([data{1}])),  ' missing']);
disp([names{2}, ': ' ,num2str(median([data{2}])),'+/-',num2str(std([data{2}])),  ' missing']);
disp([names{3}, ': ' ,num2str(median([data{3}])),'+/-',num2str(std([data{3}])),  ' missing']);


  

%% Compare endog vs rept
F = 12;
xlab = 'expressing nuclei';

names = {'Kr 2 enhancers, Kr';
         'Kr CD1, Kr';  % CD1
         'Kr CD2, Kr' % CD2
         'Kr 2 enhancers, LacZ';
         'Kr CD1, LacZ'; % CD1
         'Kr CD2, LacZ'  % CD2
         };

colordef white; 
endog = cell(1,G); 
rept = cell(1,G);
 for k=1:G;     endog{k} = endog_cnt{k}(cc14{k})./ND(cc14{k},k)/5; end
 for k=1:G;     rept{k} = rept_cnt{k}(cc14{k})./ND(cc14{k},k)/5; end

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
 Wpvals = ['p_{14} = ',num2str(pW(1,4),2), '   p_{25} = ',num2str(pW(2,5),2) , '    p_{36} = ',num2str(pW(3,6),2)  ];
 Apvals = ['p_{14} = ',num2str(pA(1,4),2), '   p_{25} = ',num2str(pA(2,5),2) , '    p_{36} = ',num2str(pA(3,6),2)  ];
 disp(['pairwise Wilcoxon rank sum:  ', Wpvals]);
 disp(['2-way ANOVA:  ',Apvals]);
 
 figure(1); clf;
 cityscape(data,names,xlab,F);
 
 figure(3); clf;
  cumhist(data,names,xlab,F);
  title(['pairwise Wilcoxon:  ' Wpvals]);
  set(gcf,'color','w');


%%  Ectopic expression rate

xlab = 'ectopic expression rate';

names = {'Kr 2 enhancers, 22C';
         'Kr distal, 22C';
         'Kr proximal, 22C'
         };

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
 Wpvals = ['p_{KI} = ',num2str(pW(1,2),2), '   p_{KJ} = ',num2str(pW(1,3),2) , '    p_{IJ} = ',num2str(pW(2,3),2)  ];
 Apvals = ['p_{12} = ',num2str(pA(1,2),2), '   p_{13} = ',num2str(pA(1,3),2) , '    p_{23} = ',num2str(pA(2,3),2)  ];
 disp(['pairwise Wilcoxon rank sum:  ', Wpvals]);
 disp(['2-way ANOVA:  ',Apvals]);
 
 figure(2); clf;
 cityscape(data,names,xlab,F);
 
 figure(4); clf;
  cumhist(data,names,xlab,F);
  title(['pairwise Wilcoxon:  ' Wpvals]);
  set(gcf,'color','w');

      
disp([names{1}, ': ' ,num2str(median([data{1}])),'+/-',num2str(std([data{1}])),  ' missing']);
disp([names{2}, ': ' ,num2str(median([data{2}])),'+/-',num2str(std([data{2}])),  ' missing']);
disp([names{3}, ': ' ,num2str(median([data{3}])),'+/-',num2str(std([data{3}])),  ' missing']);

