%%                          Plot_DotComps2_kni.m                          %%
%
% Analyzing New Shadow data
%
%
% Alistair Boettiger                                   Date Begun: 03/05/10
% Levine Lab                                     Functional Since: 10/14/10
%                                                   Last Modified: 10/18/10

%% Description
% comparison
%
%
%% Updates
% Modified 10/18/10 to also count ectopically active nuclei

%% Source Code
clear all;



folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed';

emb_roots = {'kni_2enh_22C_LacZ_kni';  
             'kni_int_22C_LacZ_kni';
             'kni_5p_22C_LacZ_kni'
            };  
          

names = {'kni 2 enhancers, 22C';
         'kni int, 22C';
         'kni 5p, 22C'
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
       %   miss_cnt{z}(n) =  length(setdiff(pts2,pts1));
      %     miss_cnt{z}(n) =  length(intersect(setdiff(pts2,pts1), ptr_nucin2));  
           [miss_cnt{z}(n), ptr_nucin2] = anlz_major_reg(folder,emb_roots{z},emb );
           miss_rate{z}(n) = miss_cnt{z}(n)/length(pts2); 
        
           ectop_cnt{z}(n) = length(intersect(setdiff(pts1,pts2),setdiff(ptr_nucin1,ptr_nucin2)')); 
        %   [lowon{z}(n), cell_var{z}(n)] =   lowon_fxn(H,handles,nin2,ptr_nucin2,[emb_roots{z},emb],1); 
          
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


% save kni_LacZ_data_2;
 %save kni_LacZ_data;
%  save kni_LacZ_data_ect; % also record ectopic expression rate.  

% % Concatinate data sets from different sessions
%[miss_cnt,miss_rate,nd,lowon] = merge_data(3,4,N,miss_cnt,miss_rate,nd,lowon);





%%
% clear all; load  kni_LacZ_data_ect;
 %clear all; load  kni_LacZ_data;

clear all; load kni_LacZ_data_2;
 
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



%%
xlab = 'fraction of missed nuclei';


plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = miss_rate{k}(cc14{k}); end
% for k=1:G;     plot_miss{k} = miss_rate{k}; end


 figure(1); clf;
 colordef black; set(gcf,'color','k');
%colordef white; set(gcf,'color','w');

x = linspace(0,1,8);  % range and number of bins for histogram
xx = linspace(0,1,100); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation
y = CompDist(plot_miss,x,xx,method,sigma,names,xlab,14);

title('cc14 embryos');

%%


plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = miss_rate{k}(cc13{k}); end
% for k=1:G;     plot_miss{k} = miss_rate{k}; end


 figure(2); clf;
% colordef black; set(gcf,'color','k');
colordef white; set(gcf,'color','w');

x = linspace(0,1,7);  % range and number of bins for histogram
xx = linspace(0,1,100); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation
y = CompDist(plot_miss,x,xx,method,sigma,names,xlab,14);

title('cc13 embryos');

%%
plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = miss_rate{k}(cc12{k}|cc11{k}); end
% for k=1:G;     plot_miss{k} = miss_rate{k}; end


 figure(3); clf;
% colordef black; set(gcf,'color','k');
colordef white; set(gcf,'color','w');

x = linspace(0,1,7);  % range and number of bins for histogram
xx = linspace(0,1,100); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation
y = CompDist(plot_miss,x,xx,method,sigma,names,xlab,14);

title('cc11 and 12 embryos');




%%  Ectopic expression rate
xlab = 'fraction of ectopic on nuclei';
plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = ectop_rate{k}(cc14{k}); end
% for k=1:G;     plot_miss{k} = miss_rate{k}; end


 figure(5); clf;
% colordef black; set(gcf,'color','k');
colordef white; set(gcf,'color','w');

x = linspace(0,1,18);  % range and number of bins for histogram
xx = linspace(0,1,100); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation
y = CompDist(plot_miss,x,xx,method,sigma,names,xlab,14);

title('cc14 embryos');




plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = ectop_rate{k}(cc13{k}); end
% for k=1:G;     plot_miss{k} = miss_rate{k}; end


 figure(6); clf;
% colordef black; set(gcf,'color','k');
colordef white; set(gcf,'color','w');

x = linspace(0,1,15);  % range and number of bins for histogram
xx = linspace(0,1,100); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation
y = CompDist(plot_miss,x,xx,method,sigma,names,xlab,14);

title('cc13 embryos');

xlab = 'fraction of ectopic on nuclei';


%% Total Expression Variability

%%
xlab = '\sigma/\mu';


plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = lowon{k}(cc13{k}); end
% for k=1:G;     plot_miss{k} = miss_rate{k}; end


 figure(1); clf;
% colordef black; set(gcf,'color','k');
colordef white; set(gcf,'color','w');

x = linspace(0,1,14);  % range and number of bins for histogram
xx = linspace(0,1,1000); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation
y = CompDist(plot_miss,x,xx,method,sigma,names,xlab,14);

title('cc13 embryos');
