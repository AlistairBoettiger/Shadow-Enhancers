%%                          Plot_DotComps.m                             %%
%
% Analyzing New Shadow data
%
%
% Alistair Boettiger                                   Date Begun: 03/05/10
% Levine Lab                                     Functional Since: 03/06/10
%                                                   Last Modified: 03/17/10

%% Description
% comparison
%
%
%% Updates


%% Source Code
clear all;

folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed';




emb_roots = { 'C33_29C_LacZ_hb';
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


xmin = .2; xmax = .8; ymin = .25; ymax = .75;
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
           lowon{z}(n) = lowon_fxn(H,handles,nin2,ptr_nucin2,emb); 
           %lowon{z}(n) = lowon_fxn(H,handles,all_nucs,pts2,nin2,Cell_bnd);  
          
           if length(H) > 2000
               im_dim = 2048;
           else
               im_dim = 1024;
           end
               
           lims =  round([xmin,xmax,ymin,ymax]*im_dim); 
           nd{z}(n) = NucDensity(cent,lims);

        catch ME
            disp(ME.message); 
            %disp(['can not find file' folder,'/',emb_roots{z},emb,'_data.mat']);
        end
    end
end




%%
 % clear all; load snail_shadow_data;

ND = cell2mat(nd); 
figure(1); clf;
plot( log2(sort(ND(:))) ,'r.');
%%
cc14 =cell(1,G); cc13 = cell(1,G); cc12 = cell(1,G); 
for z=1:G
    cc14{z} = log2(ND(:,z))>9;
    cc13{z} = log2(ND(:,z))<9  & log2(ND(:,z))>8;
    cc12{z} = log2(ND(:,z))<8;
end


%% Plot Fraction of missing nuclei distributions

xlab = 'fraction of missed nuclei';


plot_miss = cell(1,G); 
% for k=1:G;     plot_miss{k} = miss_rate{k}(cc13{k}); end
for k=1:G;     plot_miss{k} = miss_rate{k}; end


 figure(2); clf;
% colordef black; set(gcf,'color','k');
colordef white; set(gcf,'color','w');

x = linspace(0,1,15);  % range and number of bins for histogram
xx = linspace(0,1,100); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .15;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab)



%% Plot Fraction of missing nuclei distributions

xlab = 'fraction of missed nuclei';


plot_miss = cell(1,G); 
for k=1:G;     plot_miss{k} = miss_rate{k}(cc14{k}); end

 figure(1); clf;
% colordef black; set(gcf,'color','k');
colordef white; set(gcf,'color','w');

% x = linspace(0,1,8);  % range and number of bins for histogram
% xx = linspace(0,1,100); % range a number of bins for interpolated distribution
%  method = 'pcubic'; % method for interpolation
% sigma = .1;  % smoothing factor for interpolation
% CompDist(plot_miss,x,xx,method,sigma,names,xlab)

BoxDist(plot_miss,names,xlab);
set(gcf,'color','k');

V= plot_miss;
    P_var = zeros(G,G); 
    for i=1:G
        for j=1:G   
   P_var(i,j) =  log10(ranksum(V{i},V{j}));
        end
    end
    
    figure(6); clf; imagesc(P_var); colorbar;  colormap('gray');
    set(gca,'YtickLabel', str2mat(names{:}),'YTick',1:6,'fontsize',15,...
    'YMinorTick','on'); title(xlab);
    set(gcf,'color','k');

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
% %