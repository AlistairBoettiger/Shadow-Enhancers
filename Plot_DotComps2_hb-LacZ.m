%%                          Plot_DotComps.m                             %%
%
% Analyzing New Shadow data
%
%
% Alistair Boettiger                                   Date Begun: 03/05/10
% Levine Lab                                     Functional Since: 03/06/10
%                                                   Last Modified: 10/26/10

%% Description
% comparison
%
%
%% Updates
% Changed age ID code, 10/26/10

%% Source Code
clear all;

folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed';




emb_roots ={ 'C33_29C_LacZ_hb';
 'hb2enh_29C_LacZ_hb';
 'hb2enh_22C_LacZ_hb';
'C33_22C_LacZ_hb';
 'hbP_22C_LacZ_hb';
 'C55_29C_LacZ_hb'};

names = {'shadow 30C';
    '2 enhancer 30C';
    '2 enhancer 22C';
    'shadow 22C';
   ' primary 22C';
    'primary 30C'};
    
    

N = 60;
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
          miss_cnt{z}(n) =  length(intersect(setdiff(pts2,pts1), ptr_nucin2));  
         %  miss_cnt{z}(n) = anlz_major_reg(folder,emb_roots{z},emb );
           miss_rate{z}(n) = miss_cnt{z}(n)/length(pts2); 
              [lowon{z}(n)] = lowon_fxn(H,handles,nin2,ptr_nucin2,[emb_roots{z},emb],0); 
           %lowon{z}(n) = lowon_fxn(H,handles,all_nucs,pts2,nin2,Cell_bnd);  
          
           if length(H) > 2000
               im_dim = 2048;
           else
               im_dim = 1024;
           end
               
           lims =  round([xmin,xmax,ymin,ymax]*im_dim); 
           nd{z}(n) = NucDensity(cent,lims,1);

        catch ME
            disp(ME.message); 
            %disp(['can not find file' folder,'/',emb_roots{z},emb,'_data.mat']);
        end
    end
end





ND = cell2mat(nd); 
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
    cc11{z} = logage <12 ;
end

%% Plot Fraction of missing nuclei distributions

xlab = 'fraction of missed nuclei';
F = 12; % FontSize; 
ymax = .02; pts = 1000;

plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = miss_rate{k}(cc14{k}); end


figure(33); clf;  subplot(3,1,1);
colordef white; set(gcf,'color','w');

x = linspace(0,1,8);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .05;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,.5*ymax]); title('Early cc14');


%~~~~ Plot Fraction of missing nuclei distributions  ~~~~

plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = miss_rate{k}(cc13{k}); end

figure(33); subplot(3,1,2);
colordef white; set(gcf,'color','w');

x = linspace(0,1,15);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,ymax]); title('cc13');



%~~~~ Plot Fraction of missing nuclei distributions ~~~~~

plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = miss_rate{k}(cc12{k}|cc11{k} ); end


 figure(33); subplot(3,1,3);
colordef white; set(gcf,'color','w');

x = linspace(0,1,15);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .15;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,.5*ymax]);  title('cc11 & 12');


%%




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


