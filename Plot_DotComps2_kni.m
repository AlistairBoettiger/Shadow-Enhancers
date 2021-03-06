%%                          Plot_DotComps2_kni.m                          %%
%
% Analyzing New Shadow data
%
%
% Alistair Boettiger                                   Date Begun: 03/05/10
% Levine Lab                                     Functional Since: 10/14/10
%                                                   Last Modified: 02/28/11

%% Description
% comparison
%
%
%% Updates
% Modified 10/18/10 to also count ectopically active nuclei
% Modified 02/28/11 to use cityscape and cumulative sum.  


%% Source Code
clear all;



folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed';

emb_roots = {'kni_2enh_22C_LacZ_kni';  
             'kni_int_22C_LacZ_kni';
             'kni_5p_22C_LacZ_kni'
            };  
          

names = {'kni 2 enhancers, 22C';
         'kni proximal, 22C'; % intronic
         'kni distal, 22C'% 5p
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
endog_frac  = cell(1,K); 
rept_frac  = cell(1,K); 

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
    endog_frac{z}  =  zeros(N,1);   
    rept_frac{z}  =  zeros(N,1);   
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
       %  miss_cnt{z}(n) =  length(setdiff(pts2,pts1));
      %   miss_cnt{z}(n) =  length(intersect(setdiff(pts2,pts1), ptr_nucin2));  
          [miss_cnt{z}(n), ptr_nucin2] = anlz_major_reg(folder,emb_roots{z},emb );
          miss_rate{z}(n) = miss_cnt{z}(n)/length(pts2); 
          endog_cnt{z}(n) = length(pts2); 
          rept_cnt{z}(n) = length(pts1);
           
        %   ectop_cnt{z}(n) = length(setdiff(pts1,pts2));
           ectop_cnt{z}(n) = length(intersect(setdiff(pts1,pts2),setdiff(ptr_nucin1,ptr_nucin2)')); 
        %   [lowon{z}(n), cell_var{z}(n)] =   lowon_fxn(H,handles,nin2,ptr_nucin2,[emb_roots{z},emb],1); 
          
           endog_frac{z}(n) = length(intersect(ptr_nucin2,pts2))/length(ptr_nucin2); 
           rept_frac{z}(n) = length(intersect(ptr_nucin1,pts1))/length(ptr_nucin2); 
        
            [lowon{z}(n)] = lowon_fxn(H,handles,nin2,ptr_nucin2,[emb_roots{z},emb],0);  
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

save([data_folder,'kni_LacZ_data_031611']);
% save([data_folder,'kni_LacZ_data_030711']);
disp('finished');

% Avoid dependence on endogenous region finding for computing ectopic
% expression.  

% save kni_LacZ_data_2;
 %save kni_LacZ_data;
%  save kni_LacZ_data_ect; % also record ectopic expression rate.  

% % Concatinate data sets from different sessions
%[miss_cnt,miss_rate,nd,lowon] = merge_data(3,4,N,miss_cnt,miss_rate,nd,lowon);





%%
% clear all; load  kni_LacZ_data_ect;
 %clear all; load  kni_LacZ_data;
% clear all; load kni_LacZ_data_2;  % used 03/01/11


     clear all; 
data_folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Code_Data/';
%load([data_folder,'kni_LacZ_data_2']);
load([data_folder,'kni_LacZ_data_031611']);% ([data_folder,'kni_LacZ_data_030711']);
 
%%
 
ND = cell2mat(nd); 


emb_cycle = 4.9 + log2( nonzeros( sort(ND(:)) ) );
figure(10); clf; plot( emb_cycle ,'r.');

T_embs = length(nonzeros(ND(:))) ;
title(['kni embryos, N = ',num2str(T_embs)  ],'FontSize',15);
set(gca,'FontSize',15); grid on;
set(gcf,'color','w'); ylabel('log_2(nuc density)');  xlabel('embryo number'); 
ylim([10,14.99]); xlim([0,T_embs + 10]);


%%
cc14 =cell(1,G); cc13 = cell(1,G); cc12 = cell(1,G); cc11 = cell(1,G); cc10 = cell(1,G); cc9 = cell(1,G); 
for z=1:G
    logage =   4.9 + log2( ND(:,z) );
    
    cc14{z} = logage >14;
    cc13{z} = logage <14  & logage> 13;
    cc12{z}  = logage <13 & logage > 12;
    cc11{z} = logage <12 ;
end


%% Compare endog vs rept
F = 12;
xlab = 'expressing fraction of nuclei';

names = {'kni 2 enhancers, kni';
         'kni int, kni';
         'kni 5p, kni';
         'kni 2 enhancers, LacZ';
         'kni int, LacZ';
         'kni 5p, LacZ'
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

%% Fraction missed 
F = 14;
xlab = 'fraction of missed nuclei';
names = {'kni 2 enhancers, 22C';
         'kni distal, 22C';
         'kni proximal, 22C'
         };

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
 Wpvals = ['p_{GF} = ',num2str(pW(1,2),2), '   p_{GE} = ',num2str(pW(1,3),2) , '    p_{EF} = ',num2str(pW(2,3),2)  ];
 Apvals = ['p_{12} = ',num2str(pA(1,2),2), '   p_{13} = ',num2str(pA(1,3),2) , '    p_{23} = ',num2str(pA(2,3),2)  ];
 disp(['pairwise Wilcoxon rank sum:  ', Wpvals]);
 disp(['2-way ANOVA:  ',Apvals]);
 
 figure(1); clf;
 cityscape(data,names,xlab,F);
 
 figure(3); clf;
  cumhist(data,names,xlab,F);
  title(['pairwise Wilcoxon:  ' Wpvals]);
  set(gcf,'color','w');



  
disp([names{1},': ' ,num2str(median([data{1}])),'+/-',num2str(std([data{1}])),  ' missing']);
disp([names{2},': ' ,num2str(median([data{2}])),'+/-',num2str(std([data{2}])),  ' missing']);
disp([names{3},': ' ,num2str(median([data{3}])),'+/-',num2str(std([data{3}])),  ' missing']);


  

%%  Ectopic expression rate
F = 14;
xlab = 'ectopic expression rate';

plot_miss = cell(1,G); 
% for k=1:G;     plot_miss{k} = ectop_rate{k}(cc14{k}); end
 for k=1:G;     plot_miss{k} = ectop_rate{k}(cc14{k}).*ND(cc14{k},k)./rept_cnt{k}(cc14{k}) ; end
  
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
 Wpvals = ['p_{GF} = ',num2str(pW(1,2),2), '   p_{GE} = ',num2str(pW(1,3),2) , '    p_{EF} = ',num2str(pW(2,3),2)  ];
 Apvals = ['p_{12} = ',num2str(pA(1,2),2), '   p_{13} = ',num2str(pA(1,3),2) , '    p_{23} = ',num2str(pA(2,3),2)  ];
 disp(['pairwise Wilcoxon rank sum:  ', Wpvals]);
 disp(['2-way ANOVA:  ',Apvals]);
 
 figure(2); clf;
 cityscape(data,names,xlab,F);
 
 figure(4); clf;
  cumhist(data,names,xlab,F);
  title(['pairwise Wilcoxon:  ' Wpvals]);
  set(gcf,'color','w');

disp([names{1},': ' ,num2str(median([data{1}])),'+/-',num2str(std([data{1}])),  ' ectopic']);
disp([names{2},': ' ,num2str(median([data{2}])),'+/-',num2str(std([data{2}])),  ' ectopic']);
disp([names{3},': ' ,num2str(median([data{3}])),'+/-',num2str(std([data{3}])),  ' ectopic']);

%%   'Completeness' of Pattern
F=14;
xlab = 'expressing nuclei';

names = {'kni 2 enhancers, kni';
         'kni int, kni';
         'kni 5p, kni';
         'kni 2 enhancers, LacZ';
         'kni int, LacZ';
         'kni 5p, LacZ'
         };

colordef white; 
endog = cell(1,G); 
rept = cell(1,G);
 for k=1:G;     endog{k} = endog_frac{k}(cc14{k}); end
 for k=1:G;     rept{k} = rept_frac{k}(cc14{k}); end

  data = cat(2,endog,rept);    
  Ts = length(data);% number of tracks
  pW = zeros(Ts);
  pA = zeros(Ts); 
  for i=1:Ts
    for j = 1:Ts
     pW(i,j) = ranksum(data{i},data{j});   % Wilcox Rank Sum
    % pA(i,j)=anovan([data{i}',data{j}'],{[zeros(1,length(data{i})),ones(1,length(data{j}))]},'display','off'); % 2-way ANOVA
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

  
  disp([names{1}, ': ' ,num2str(median([data{1}])),'+/-',num2str(std([data{1}])),  ' missing']);
disp([names{2}, ': ' ,num2str(median([data{2}])),'+/-',num2str(std([data{2}])),  ' missing']);
disp([names{3}, ': ' ,num2str(median([data{3}])),'+/-',num2str(std([data{3}])),  ' missing']);
disp([names{4}, ': ' ,num2str(median([data{4}])),'+/-',num2str(std([data{4}])),  ' missing']);
disp([names{5}, ': ' ,num2str(median([data{5}])),'+/-',num2str(std([data{5}])),  ' missing']);
disp([names{6}, ': ' ,num2str(median([data{6}])),'+/-',num2str(std([data{6}])),  ' missing']);




