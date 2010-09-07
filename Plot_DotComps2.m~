%%                          Plot_DotComps.m                             %%
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
            '2 enhancers, 29C';
            '2 enhancers, 22C';
            'no shadow, 29C';
            'no shadow, 22C';
            'no primary, 29C';
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

save sog_shadow_data

%% visualize distribution of nuclear densities
ND = cell2mat(nd); 
figure(5); clf;
plot( log2(sort(ND(:))) ,'r.');

%% Separate out age structure by division cycle
cc14 =cell(1,K); cc13 = cell(1,K); cc12 = cell(1,K); 
for z=1:K
    cc14{z} = log2(ND(:,z))>8;
    cc13{z} = log2(ND(:,z))<8 & log2(ND(:,z))>7;
    cc12{z} = log2(ND(:,z))<7;
end

%% Plot Fraction of missing nuclei distributions

xlab = 'fraction of missed nuclei';


plot_miss = cell(1,G); 
for k=1:G;     plot_miss{k} = miss_rate{k}(cc14{k}); end
 plot_miss{4} = [plot_miss{4}; miss_rate{7}(cc14{7})  ]; 

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

%% Plot Total mRNA variability distribuitons 

xlab = 'variability in total transcript (\sigma/\mu)';


plot_lowon = cell(1,G); 
for k=1:G; plot_lowon{k} = lowon{k}(cc14{k}); end
plot_lowon{4} = [plot_lowon{4}; lowon{7}(cc14{7})]; 

figure(2); clf;
% colordef black; set(gcf,'color','k');
colordef white; set(gcf,'color','w');

%  x = linspace(0,1,20);
%   xx = linspace(0,1,100); 
%  method = 'pcubic';
% sigma = .1; 
% CompDist(plot_lowon,x,xx,method,sigma,names,xlab)

BoxDist(plot_lowon,names,xlab); xlim([.3,.75]);
set(gcf,'color','k');

V= plot_lowon;
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

%% Variability between immidiate neighbors 
 
xlab = 'variability in total transcript among neighbors \sigma/mu';
% x = linspace(0,1,33);
%   xx = linspace(0,1,100); 
%  method = 'pcubic';
% sigma = .1; 

plot_cell_var = cell(1,G); 
for k=1:G;  plot_cell_var{k} = cell_var{k}(cc14{k}); end
plot_cell_var{4} = [plot_cell_var{4}; cell_var{7}(cc14{7})];

 figure(3); clf;
colordef black; set(gcf,'color','k');
BoxDist(plot_cell_var,names,xlab);
set(gcf,'color','k');
% CompDist(plot_cell_var,x,xx,method,sigma,names,xlab);

