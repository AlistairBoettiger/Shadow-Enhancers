%%                          Plot_DotComps2_ConcDept_hb.m                 %%
%
% Analyzing New Shadow data
%
%
% Alistair Boettiger                                   Date Begun: 03/05/10
% Levine Lab                                     Functional Since: 03/06/10
%                                                   Last Modified: 03/28/11

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
% Modified 12/09 to fix merge data error and to streamline region
% processing.  
% Modified to save data 12/14
% Modified 12/15 to add additional output: compare center vs boundary
% Modified 03/21/11 to compare cell-by-cell
% Modified 03/28/11 to update figure plotting
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
            'BAC01b_22C_y_hb'; % 10
            'BAC01b_30C_y_hb'; % 11
            'BAC02b_30C_y_hb'; % 12
            };  
          %1 +6, 2+9, 3+4, 
     
 names = {'control 22C';
          'control 30C';
          'no shadow 22C';
          'no shadow 30C';
          'no primary 22C';
          'no primary 30C'
          };
      

N = 100;
Zs = length(emb_roots); 
G= length(names);


Layer_misses = cell(N,Zs); 
nd = cell(1,Zs); 
age_table = cell(1,Zs);


xmin = .2; xmax = .9; ymin = .15; ymax = .4;
% as fractions of the original image dimensions.  

%%
for z=1:Zs % k=2;   z =9
    for n=   1:N % n = 12   n= 8
        if n<10
            emb = ['0',num2str(n)];
        else
            emb = num2str(n);
        end

        try
            fname = [folder,'/',emb_roots{z},emb,'_data.mat'];
            load(fname);   
       %       [miss_rate{z,1}(n),miss_rate{z,2}(n),miss_rate{z,3}(n)] =  divide_regs(L2,H,pts1,pts2,ptr_nucin2,handles.In,0);
               
           if length(H) > 2000
               im_dim = 2048;
           else
               im_dim = 1024;
           end    
           lims =  round([xmin,xmax,ymin,ymax]*im_dim); 
            nd{z}(n) = NucDensity(cent,lims,0);
            age_table{z}{n,1} = fname; % for labeling purposes  
            age_table{z}{n,2} = nd{z}(n); 
            
            Layer_misses{n,z} = divide_regs(Reg1,Reg2,L2n1,H,cellbords,1); 
            
        catch ME
            disp(ME.message); 
        end
    end
end

close all; 

% save hb_SDreg_12-15-10

%% Merge Data
foff = cell(6,1); 

 foff{1} = [Layer_misses(1,:), Layer_misses(6,:)] ;   % 2 enhancer 22C
 foff{2} =  Layer_misses(7,:) ;  % 2 enhancer 30C
 foff{3} = [Layer_misses(2,:), Layer_misses(9,:)];   % no shadow 22C
 foff{4} = [Layer_misses(3,:), Layer_misses(4,:), Layer_misses(12,:)]; % no shadow 30C
 foff{5} = [Layer_misses(5,:), Layer_misses(10,:)];  % no primary 22C
 foff{6} = [Layer_misses(8,:), Layer_misses(11,:)]; % no primary 30C

 Ls = size(Layer_misses,2)
 
 Nnuc = cell(1,6); 
  Nnuc{1} = [nd{1}, zeros(1,Ls-length(nd{1})), nd{6}, zeros(1,Ls-length(nd{6}))] ;  % 2 enhancer 22C
  Nnuc{2} = [nd{7}, zeros(1,Ls-length(nd{7}))]; % 2 enhancer 30C
  Nnuc{3} = [nd{2}, zeros(1,Ls-length(nd{2})), nd{9}, zeros(1,Ls-length(nd{9}))]; % no shadow 22C
   Nnuc{4} = [nd{3}, zeros(1,Ls-length(nd{3})), nd{4}, zeros(1,Ls-length(nd{4})), nd{12}, zeros(1,Ls-length(nd{12}))]; % no shadow 30C
   Nnuc{5} = [nd{5}, zeros(1,Ls-length(nd{5})), nd{10}, zeros(1,Ls-length(nd{10}))]; % no primary 22C
    Nnuc{6} = [nd{8}, zeros(1,Ls-length(nd{8})), nd{11}, zeros(1,Ls-length(nd{11}))]; % no primary 30C
    
% square off entries
for k=1:G
        foff{k} = [foff{k},cell(1,300-length(foff{k}))];
        data = Nnuc{k}';
        Nnuc{k} = [data; zeros(300-length(data),1)];
 end
     

%%

ND = cell2mat(Nnuc); 
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
    cc11{z} = logage <12 & logage > 0 ;
end

%%


z = 6
emb_set = foff{z}(cc14{z});

N = length(emb_set); 
miss_freq = zeros(N,40); 
emb_norm = zeros(N,40); 
for n=1:N;
    L = length(emb_set{n});
    miss_freq(n,1:L) = emb_set{n};
    emb_norm(n,1:L) = emb_set{n}>0; 
end

F = 16;

miss_dist_hist = sum(miss_freq)./sum(emb_norm); 
figure(1); clf; 
colordef black;
bar(miss_dist_hist);
set(gcf,'color','k');
ylabel('fraction of missing nulcei','FontSize',F); 
xlabel('distance from hb-boundary (cells)','FontSize',F); 
title([names{z},'  cc14   N = ',num2str(N)],'FontSize',F);
set(gca,'FontSize',F);
 xlim([0,30.5]); ylim([0,.6]); 
 
 
 %% new plotting cc14
 
 clear all;
 
data_folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Code_Data/';

  load([data_folder,'hb_ConcDept']);
  
  
  F = 16;

    figure(1); clf; 
    colordef black;
    set(gcf,'color','k');

      
    data = zeros(G,30); 
    C = flipud(hsv(G));
    leg_lab = cell(G,1);  
  
for z = 1:6
    emb_set = foff{z}(cc14{z});

    N = length(emb_set); 
    miss_freq = zeros(N,40); 
    emb_norm = zeros(N,40); 
    for n=1:N;
        L = length(emb_set{n});
        miss_freq(n,1:L) = emb_set{n};
        emb_norm(n,1:L) = emb_set{n}>0; 
    end

    miss_dist_hist = sum(miss_freq)./sum(emb_norm); 
    data(z,1:30) = miss_dist_hist(1:30); 
    leg_lab{z} = ['hb ', names{z},'  cc14   N = ',num2str(N)];
end

plot(data([1,3],:)','.','MarkerSize',10); legend(leg_lab([1,3]));

figure(2); clf; set(gcf,'color','k');
bar(data([1,3,5],:)'); legend(leg_lab([1,3,5]));
set(gca,'FontSize',F);
xlabel('Distance from edge of expression (cells)');
ylabel('mean fraction of inactive nuclei'); xlim([0,30.5]); ylim([0,.25]); 

%% cc13


    figure(1); clf; 
    colordef black;
    set(gcf,'color','k');

      
    data = zeros(G,20); 
    C = flipud(jet(G));
    leg_lab = cell(G,1);  
  
for z = 1:6
    emb_set = foff{z}(cc13{z});

    N = length(emb_set); 
    miss_freq = zeros(N,40); 
    emb_norm = zeros(N,40); 
    for n=1:N;
        L = length(emb_set{n});
        miss_freq(n,1:L) = emb_set{n};
        emb_norm(n,1:L) = emb_set{n}>0; 
    end

    miss_dist_hist = sum(miss_freq)./sum(emb_norm); 
    data(z,1:20) = miss_dist_hist(1:20); 
    leg_lab{z} = ['hb ', names{z},'  cc13   N = ',num2str(N)];
end

plot(data([1,3,5],:)','.','MarkerSize',10); legend(leg_lab([1,3,5]));

%%
plot(data([1,3],:)','.-','MarkerSize',10); legend(leg_lab([1,3]));

figure(2); clf;
bar(data([1,3,5],:)'); legend(leg_lab([1,3,5]));
set(gca,'FontSize',F);
xlabel('Distance from edge of expression (cells)');
ylabel('mean fraction of inactive nuclei'); xlim([0,16.5]);

%%
figure(3); clf; % set(gcf,'color','k');
set(gcf,'color','w'); colordef white; 
Pdata = data([1,3,5],1:16);
%Pdata(1,:) = NaN*Pdata(1,:);
%Pdata(2,:) = NaN*Pdata(2,:);
len = length(Pdata);



bar(flipud(1:len)',fliplr(Pdata)'); legend(leg_lab([1,3,5]),'Location','NorthWest');
set(gca,'FontSize',F);
set(gca,'Xtick',1:len,'XtickLabel',num2str(fliplr(1:len)'));
xlabel('Distance from edge of expression (cells)');
ylabel('mean fraction of inactive nuclei'); xlim([0,17]); ylim([0,.35]); 



