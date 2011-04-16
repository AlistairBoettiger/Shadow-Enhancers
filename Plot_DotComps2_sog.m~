%%                          Plot_DotComps2.m                             %%
%
% Analyzing New Shadow data
%
%
% Alistair Boettiger                                   Date Begun: 03/05/10
% Levine Lab                                     Functional Since: 03/06/10
%                                                   Last Modified: 04/02/10

%% Description
% comparison
%
%
%% Updates
% rewritten on 04/02/10 to streamline addition of new tracks using cell
% notation to replace repetitive coding.  

%% Source Code
clear all;

folder =  '/Volumes/Data/Lab Data/Shadow_data';

emb_roots = { 'sog2enh_29C_LacZ_sog';
              'sog2enh_22C_LacZ_sog';
              'SogP30C_LacZ_sog';
              'sogPxYW_22C_LacZ_sog'; 
              'sogSxYW_29C_LacZ_sog';
              'sogSxYW_22C_LacZ_sog' ;
              'sogP_C81_22C_LacZ_sog';};
   

names = {
            '2 enhancers, 30C';
            '2 enhancers, 22C';
            'no shadow, 30C';
            'no shadow, 22C';
            'no primary, 30C';
            'no primary, 22C';};

poor_images = { ['02','07'];
                ['09','07','06','03']; 
                '06'; 
                [];
                [];
                [];
                ['15','13','12','10','09']};
                
        
N = 40;
K = length(emb_roots); 
G = length(names); 


miss_cnt = cell(1,K); 
miss_rate = cell(1,K); 
lowon = cell(1,K); 
cell_var = cell(1,K); 
nd = cell(1,K); 

for z=1:K
    miss_cnt{z} = zeros(N,1);
    miss_rate{z} = zeros(N,1); 
    lowon{z} = zeros(N,1); 
    cell_var{z} = zeros(N,1); 
    nd{z} = zeros(N,1);
end


xmin =200; xmax =800; ymin=150; ymax =300;
for z=1:K % z=7;
    for n=1:N
        if n<10
            emb = ['0',num2str(n)];
        else
            emb = num2str(n);
        end

        if isempty( strfind(poor_images{z},emb) ) == 1;
             
            try                     
            load([folder,'/',emb_roots{z},emb,'_data.mat']);   
            % get the indices of all nuclei in green that are not also red.  
            % require these nuclei also fall in the 'region' for red nuclei.  
           % s29_miss_cnt(n) =  length(intersect(setdiff(pts2,pts1), ptr_nucin2));
               miss_cnt{z}(n) =  length(unique( intersect(setdiff(pts2,pts1), ptr_nucin2 )))  ;
               miss_rate{z}(n) = miss_cnt{z}(n)/length( unique(pts2)); 
               [lowon{z}(n), cell_var{z}(n)] = lowon_fxn(H,handles,nin2,ptr_nucin2,[emb_roots{z},emb]);  
               nd{z}(n) = NucDensity(cent,[xmin,xmax,ymin,ymax]);
            catch ME
                disp(ME.message); 
                %disp(['can not find file' folder,'/',emb_roots{z},emb,'_data.mat']);
            end
        end
    end
end

% save sog_shadow_data

% load sog_shadow_data
%[miss_cnt,miss_rate,nd,lowon] = merge_data(4,7,N,miss_cnt,miss_rate,nd,lowon);
%% visualize distribution of nuclear densities

 clear all;
load('/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Code_Data/sog_shadow_data'); 
 
ND = cell2mat(nd); 
age_offset = 5.9;

emb_cycle = age_offset + log2( nonzeros( sort(ND(:)) ) );  % density data computed using smaller area
figure(10); clf; plot( emb_cycle ,'r.');

T_embs = length(nonzeros(ND(:))) ;
title(['kni embryos, N = ',num2str(T_embs)  ],'FontSize',15);
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


colordef white; 
endog = cell(1,G); 
rept = cell(1,G);
 for k=1:G;     data{k} = plot_miss{k}(cc14{k}); end


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

 
 figure(3); clf;
  cumhist(data,names,xlab,F);
  title(['pairwise Wilcoxon:  ' Wpvals]);
  set(gcf,'color','w');

%% Fraction missed 
F = 14;
xlab = 'fraction of missed nuclei';


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
 
 
 figure(5); clf; imagesc(log10(pW)); colormap(bone);
   set(gcf,'color','w');
   set(gca,'YTickLabel',names); colorbar; 
   title('log p-values');

 
 figure(3); clf;
  cumhist(data,names,xlab,F);
 % title(['pairwise Wilcoxon:  ' Wpvals]);
  set(gcf,'color','w');



  disp([names{1}, ': ' ,num2str(median([data{1}])),'+/-',num2str(std([data{1}])),  ' missing']);
disp([names{2}, ': ' ,num2str(median([data{2}])),'+/-',num2str(std([data{2}])),  ' missing']);
disp([names{3}, ': ' ,num2str(median([data{3}])),'+/-',num2str(std([data{3}])),  ' missing']);
disp([names{4}, ': ' ,num2str(median([data{4}])),'+/-',num2str(std([data{4}])),  ' missing']);
disp([names{5}, ': ' ,num2str(median([data{5}])),'+/-',num2str(std([data{5}])),  ' missing']);
disp([names{6}, ': ' ,num2str(median([data{6}])),'+/-',num2str(std([data{6}])),  ' missing']);

%%
xlab = 'variability in total transcript (\sigma/\mu)';


plot_lowon = cell(1,G); 
for k=1:G; plot_lowon{k} = (lowon{k}(cc14{k})); end
 data = plot_lowon;    

%for k=1:G;  data{k} = cell_var{k}(cc14{k}); end



  Ts = length(data);% number of tracks
  pW = zeros(Ts);
  pA = zeros(Ts); 
  for i=1:Ts
    for j = 1:Ts
     pW(i,j) = ranksum(data{i},data{j});   % Wilcox Rank Sum
     pA(i,j)=anovan([data{i}',data{j}'],{[zeros(1,length(data{i})),ones(1,length(data{j}))]},'display','off'); % 2-way ANOVA
    end
  end
 Wpvals = ['p_{} = ',num2str(pW(1,2),2), '   p_{} = ',num2str(pW(1,3),2) , '    p_{} = ',num2str(pW(2,3),2)  ];
 Apvals = ['p_{12} = ',num2str(pA(1,2),2), '   p_{13} = ',num2str(pA(1,3),2) , '    p_{23} = ',num2str(pA(2,3),2)  ];
 disp(['pairwise Wilcoxon rank sum:  ', Wpvals]);
 disp(['2-way ANOVA:  ',Apvals]);
 
 
 figure(5); clf; imagesc(log10(pW)); colormap(bone);
   set(gcf,'color','w');
   set(gca,'YTickLabel',names); colorbar; 
   title('log p-values');

 
 figure(3); clf;
  cumhist(data([1,3,5]),names([1,3,5]),xlab,F);
  title(['pairwise Wilcoxon:  ' ['p_{12} = ',num2str(pW(1,3),2), '   p_{13} = ',num2str(pW(1,5),2) , '    p_{23} = ',num2str(pW(3,5),2)  ];]);
  set(gcf,'color','w');


   figure(4); clf;
  cumhist(data([2,4,6]),names([2,4,6]),xlab,F);

  title(['pairwise Wilcoxon:  ' ['p_{12} = ',num2str(pW(2,4),2), '   p_{13} = ',num2str(pW(2,6),2) , '    p_{23} = ',num2str(pW(4,6),2)  ];]);
  set(gcf,'color','w');

  

  disp([names{1}, ': ' ,num2str(median([data{1}])),'+/-',num2str(std([data{1}])),  ' missing']);
disp([names{2}, ': ' ,num2str(median([data{2}])),'+/-',num2str(std([data{2}])),  ' missing']);
disp([names{3}, ': ' ,num2str(median([data{3}])),'+/-',num2str(std([data{3}])),  ' missing']);
disp([names{4}, ': ' ,num2str(median([data{4}])),'+/-',num2str(std([data{4}])),  ' missing']);
disp([names{5}, ': ' ,num2str(median([data{5}])),'+/-',num2str(std([data{5}])),  ' missing']);
disp([names{6}, ': ' ,num2str(median([data{6}])),'+/-',num2str(std([data{6}])),  ' missing']);

