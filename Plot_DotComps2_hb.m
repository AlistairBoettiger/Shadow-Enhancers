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

emb_roots = {'MP09_22C_y_hb';
            'MP02_22C_y_hb';
            'MP02_30C_y_hb';
            'MP02_30C_LacZ_hb'   % yes it's actually yellow
            };  
          

names = {'2 enhancers, 22C';
         'no shadow, 22C';
         'no shadow, 30C'
         'no shadow, 30C'
         };



s29_root = 'C33_29C_LacZ_hb';
sp29_root = 'hb2enh_29C_LacZ_hb';
sp22_root = 'hb2enh_22C_LacZ_hb';
s22_root = 'C33_22C_LacZ_hb';
p22_root = 'hbP_22C_LacZ_hb';
p29_root = 'C55_29C_LacZ_hb';


N = 90;
K = length(emb_roots); 
G= length(names);


     age_table = cell(1,K);

miss_cnt = cell(1,K); 
miss_rate = cell(1,K); 
nd = cell(1,K); 
lowon = cell(1,K); 

for z=1:K
    miss_cnt{z} = zeros(N,1);
    miss_rate{z} = zeros(N,1); 
    lowon{z} = zeros(N,1); 
    nd{z} = zeros(N,1);
    age_table{z} = cell(N,2);
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
          % miss_cnt{z}(n) =  length(intersect(setdiff(pts2,pts1), ptr_nucin2));  
           miss_cnt{z}(n) = anlz_major_reg(folder,emb_roots{z},emb );
           miss_rate{z}(n) = miss_cnt{z}(n)/length(pts2); 
           lowon{z}(n) = lowon_fxn(H,handles,nin2,ptr_nucin2,emb); 
          
          
           if length(H) > 2000
               im_dim = 2048;
           else
               im_dim = 1024;
           end
               
           lims =  round([xmin,xmax,ymin,ymax]*im_dim); 
           nd{z}(n) = NucDensity(cent,lims);
            age_table{z}{n,1} = [folder,'/',emb_roots{z},emb,'_data.mat']; %  
            age_table{z}{n,2} = nd{z}(n); 

        catch ME
            disp(ME.message); 
            %disp(['can not find file' folder,'/',emb_roots{z},emb,'_data.mat']);
        end
    end
end

close all; 



 %  save hb_shadow_yellow_data;   % load hb_shadow_yellow_data;
[miss_cnt,miss_rate,nd,lowon] = merge_data(3,4,N,miss_cnt,miss_rate,nd,lowon);


G=3;
%%
 % clear all; load snail_shadow_data;

ND = cell2mat(nd); 
figure(1); clf;
plot( log2(sort(ND(:))) ,'r.');

figure(2); clf;
plot(sort(ND(:)),'r.');

%%
cc14 =cell(1,G); cc13 = cell(1,G); cc12 = cell(1,G); cc11 = cell(1,G); cc10 = cell(1,G); cc9 = cell(1,G); 
for z=1:G
    logage =  log2(ND(:,z)) ;
    
    cc14{z} = logage >9;
    cc13{z} = logage <9  & logage>8;
    cc12{z}  = logage <8 & logage > 7.5;
    cc11{z} = logage <7.5 & logage > 7;
    cc10{z} =  logage < 7 & logage > 6.5;  
    cc9{z} = logage < 6.5;
end

%% 

% disp all embryos 

for z = 1:K
    for n=1:N
        xmin = .2; xmax = .8; ymin = .1; ymax = .75;
            if length(H) > 1500
               im_dim = 2048;
           else
               im_dim = 1024;
           end
        
         lims =  round([xmin,xmax,ymin,ymax/2]*im_dim); 
     

        
      %  logage = log2(age_table{z}{n,2});
        if cc9{z}(n) ==1;  % logage < 8
            try
            load(age_table{z}{n,1});
           figure(3);  clf; imshow(handles.It);   pause(2);    
        
           figure(4); clf; 
            set(gcf,'color','k'); 
 subplot(2,1,1); 
     Io = uint8(zeros(h,w,3));
     Io(:,:,1) = uint8(255*L1);
     Io(:,:,3) = 30*handles.In;
     Ired = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
     imshow(Ired); hold on;
      %plot(rna_x1,rna_y1,'r.');  
 subplot(2,1,2); 
     Io = uint8(zeros(h,w,3));
     Io(:,:,2) = uint8(255*L2);
     Io(:,:,3) = 30*handles.In;
     Igreen = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
     imshow(Igreen); hold on;
           
           
            catch 
            end
             % title([emb_roots{z}, '   cell cycle ', num2str( round(logage + 3 ))],'FontSize',12);
          
          %  nd{z}(n) = NucDensity(cent,lims,1);
           %  age_table{z}{n,2} = nd{z}(n);  
         %   pause(.02); 
    
        end
    end
end





%% Plot Fraction of missing nuclei distributions

xlab = 'fraction of missed nuclei';


plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = miss_rate{k}(cc9{k}); end
%for k=1:G;     plot_miss{k} = miss_rate{k}; end


 figure(2); clf;
% colordef black; set(gcf,'color','k');
colordef white; set(gcf,'color','w');

x = linspace(0,1,25);  % range and number of bins for histogram
xx = linspace(0,1,100); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab)


%%
xlab = 'fraction of missed nuclei';


plot_miss = cell(1,G); 
 for k=1:G;     plot_miss{k} = miss_rate{k}(cc11{k}|cc12{k}); end
%for k=1:G;     plot_miss{k} = miss_rate{k}; end


 figure(1); clf;
% colordef black; set(gcf,'color','k');
colordef white; set(gcf,'color','w');

x = linspace(0,1,25);  % range and number of bins for histogram
xx = linspace(0,1,100); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab)

%% Plot Fraction of missing nuclei distributions

xlab = 'fraction of missed nuclei';


plot_miss = cell(1,G); 
for k=1:G;     plot_miss{k} = miss_rate{k}(cc11{k}); end

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