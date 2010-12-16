%%                          Plot_DotComps2_hb_regs.m                             %%
%
% Analyzing New Shadow data
%
%
% Alistair Boettiger                                   Date Begun: 03/05/10
% Levine Lab                                     Functional Since: 03/06/10
%                                                   Last Modified: 12/09/10

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
% Modified to save data
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
     

%      % all names
% names = {'2 enhancers, 22C';
%          'no shadow, 22C';
%          'no shadow, 30C';
%          'no shadow, 30C, b'; 
%          'no primary, 22C';
%          '2 enh 22C';
%          '2 enh 30C';
%          'no primary 30C';
%          'no shadow, 22C, b';
%          'no primary 22C, b';
%          'no primary 30C, b';
%          'no shadow 30C, c'
%          };


 names = {'control 22C';
          'control 30C';
          'no shadow 22C';
          'no shadow 30C';
          'no primary 22C';
          'no primary 30C'
          };
      

N = 100;
K = length(emb_roots); 
G= length(names);



% miss_rate1 = cell(1,K); 
% miss_rate2 = cell(1,K); 
% miss_rate3 = cell(1,K); 
% nd = cell(1,K); 
% age_table = cell(1,K);
% 
% for z=1:K
%     miss_rate1{z} = zeros(N,1);
%     miss_rate2{z} = zeros(N,1);
%     miss_rate3{z} = zeros(N,1);
%     nd{z} = zeros(N,1);
% end


miss_rate = cell(K,3); 
nd = cell(1,K); 
age_table = cell(1,K);
for z=1:K
    
    for r=1:3
        miss_rate{z,r} = zeros(N,1);
    end
    nd{z} = zeros(N,1);
end


xmin = .2; xmax = .9; ymin = .15; ymax = .4;
% as fractions of the original image dimensions.  

%%
for z=1:K % k=2;
    for n=   1:N
        if n<10
            emb = ['0',num2str(n)];
        else
            emb = num2str(n);
        end

        try
        load([folder,'/',emb_roots{z},emb,'_data.mat']);   
             %[miss_rate1{z}(n),miss_rate2{z}(n),miss_rate3{z}(n)] =  divide_regs(L2,H,pts1,pts2,ptr_nucin2,handles.In);
              [miss_rate{z,1}(n),miss_rate{z,2}(n),miss_rate{z,3}(n)] =  divide_regs(L2,H,pts1,pts2,ptr_nucin2,handles.In,0);
           if length(H) > 2000
               im_dim = 2048;
           else
               im_dim = 1024;
           end
               
           lims =  round([xmin,xmax,ymin,ymax]*im_dim); 
           nd{z}(n) = NucDensity(cent,lims,0);
            age_table{z}{n,1} = [folder,'/',emb_roots{z},emb,'_data.mat']; %  
            age_table{z}{n,2} = nd{z}(n); 
        catch ME
            disp(ME.message); 
            %disp(['can not find file' folder,'/',emb_roots{z},emb,'_data.mat']);
        end
    end
end

close all; 

save hb_SDreg_12-15-10


%%

 % clear all; load hb_SDreg_12-15-10
 % Merge data
 
%  %
%       % all names
% names = {'2 enhancers, 22C';  % 1
%          'no shadow, 22C';    % 2
%          'no shadow, 30C';    % 3
%          'no shadow, 30C, b'; % 4
%          'no primary, 22C';   % 5
%          '2 enh 22C';         % 6
%          '2 enh 30C';         % 7
%          'no primary 30C';    % 8
%          'no shadow, 22C, b'; % 9
%          'no primary 22C, b'; % 10 
%          'no primary 30C, b'; % 11
%          'no shadow 30C, c'   % 12
%          };

foff = cell(6,3);
 for r = 1:3
 foff{1,r} = [miss_rate{1,r}; miss_rate{6,r}] ;   % 2 enhancer 22C
 foff{2,r} = miss_rate{7,r};  % 2 enhancer 30C
 foff{3,r} = [miss_rate{2,r}; miss_rate{9,r}];   % no shadow 22C
 foff{4,r} = [miss_rate{3,r}; miss_rate{4,r}; miss_rate{12,r}]; % no shadow 30C
 foff{5,r} = [miss_rate{5,r}; miss_rate{10,r}];  % no primary 22C
 foff{6,r} = [miss_rate{8,r}; miss_rate{11,r}]; % no primary 30C
 end
 
 Nnuc = cell(1,6); 
  Nnuc{1} = [nd{1}; nd{6}] ;  % 2 enhancer 22C
  Nnuc{2} = nd{7}; % 2 enhancer 30C
  Nnuc{3} = [nd{2}; nd{9}]; % no shadow 22C
   Nnuc{4} = [nd{3}; nd{4}; nd{12}]; % no shadow 30C
   Nnuc{5} = [nd{5}; nd{10}]; % no primary 22C
    Nnuc{6} = [nd{8}; nd{11}]; % no primary 30C


 

 
 for k=1:G
     for r=1:3
        data = nonzeros(foff{k,r});
        foff{k,r} = [data; zeros(200-length(data),1)];
     end
    data = nonzeros(Nnuc{k});
    Nnuc{k} = [data; zeros(200-length(data),1)];
 end
     

 

%%
 % clear all; load hb_SD-region_data; 

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
    for r=1:3
        foff{z,r}(foff{z,r}==Inf) = 0;  % remove infs as well 
    end
end

%% Plot Fraction of missing nuclei distributions
close all;


fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Results/';
reg = {'Anterior '; 'Central'; 'Boundary'};
xlab = 'fraction of missed nuclei';
F = 12; % FontSize; 
ymax = .02; pts = 1000;

for r = 1:3

    f1=figure(30+r); clf;  subplot(3,1,1);
    set(f1,'Position',[50,550,400,1000]);
    %~~~~ Plot Fraction of missing nuclei distributions  ~~~~
plot_miss14 = cell(1,G); 
 for k=1:G;     plot_miss14{k} = foff{k,r}(cc14{k}); end

colordef white; set(gcf,'color','w');
x = linspace(0,1,15);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .01;  % smoothing factor for interpolation
CompDist(plot_miss14,x,xx,method,sigma,names,xlab,F)
ylim([0,.5*ymax]); title([reg{r}, ' early cc14']);

f2 = figure(20+r); subplot(3,1,1); 
   set(f2,'Position',[450,550,400,1000]);
BoxDist(plot_miss14,names,[reg{r}, ' early cc14'],{'30C','22C'});
xlim([0,1]);

  P = zeros(G);
    for j = 1:G
        for k=1:G
        P(j,k) = ranksum(plot_miss14{j},plot_miss14{k});
        end
    end
       f3=figure(10+r); subplot(3,1,1); imagesc(P); colormap(jet); colorbar; colordef white; set(gca,'color','w');
          set(f3,'Position',[850,550,400,1000]);
     set(gca,'YtickLabel', str2mat(names{:}),'YTick',1:6,'fontsize',15,...
        'YMinorTick','on'); title([reg{r}, ' early cc14']); 
        set(gcf,'color','w');


%~~~~ Plot Fraction of missing nuclei distributions  ~~~~
plot_miss13 = cell(1,G); 
for k=1:G;     plot_miss13{k} = foff{k,r}(cc13{k}); end
figure(30+r); subplot(3,1,2);
colordef white; set(gcf,'color','w');
x = linspace(0,1,20);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .05;  % smoothing factor for interpolation
CompDist(plot_miss13,x,xx,method,sigma,names,xlab,F)
ylim([0,ymax]); title([reg{r}, ' cc13']);

figure(20+r); subplot(3,1,2); 
BoxDist(plot_miss13,names,[reg{r}, ' cc13'],{'30C','22C'});
xlim([0,1]);


box_data13{r} = plot_miss13;

  P = zeros(G);
    for j = 1:G
        for k=1:G
        P(j,k) = ranksum(plot_miss13{j},plot_miss13{k});
        end
    end
    figure(10+r); subplot(3,1,2); imagesc(P); colormap(jet); colorbar; colordef white; set(gca,'color','w');
     set(gca,'YtickLabel', str2mat(names{:}),'YTick',1:6,'fontsize',15,...
        'YMinorTick','on'); title([reg{r}, ' cc13']); 
        set(gcf,'color','w');


%~~~~ Plot Fraction of missing nuclei distributions ~~~~~
plot_miss12 = cell(1,G); 
for k=1:G;     plot_miss12{k} = foff{k,r}(cc12{k}|cc11{k} ); end
 figure(30+r); subplot(3,1,3);
colordef white; set(gcf,'color','w');
x = linspace(0,1,15);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .05;  % smoothing factor for interpolation
CompDist(plot_miss12,x,xx,method,sigma,names,xlab,F)
ylim([0,.5*ymax]);  title([reg{r}, ' cc11 & 12']);

figure(20+r); subplot(3,1,3); 
BoxDist(plot_miss12,names,[reg{r}, ' cc11 & 12'],{'30C','22C'});
xlim([0,1]);

box_data12{r} = plot_miss12;

    P = zeros(G);
    for j = 1:G
        for k=1:G
        P(j,k) = ranksum(plot_miss12{j},plot_miss12{k});
        end
    end
    figure(10+r); subplot(3,1,3); imagesc(P); colormap(jet); colorbar; colordef white; set(gca,'color','w');
     set(gca,'YtickLabel', str2mat(names{:}),'YTick',1:6,'fontsize',15,...
        'YMinorTick','on'); title([reg{r}, ' cc11 & 12']); 
        set(gcf,'color','w');
    
% Export Data        
%      saveas(f1,[fout,'hb_graphs_',reg{r},'.ai'],'ai');
%      saveas(f2,[fout,'hb_bars_',reg{r},'.ai'],'ai');
%      saveas(f1,[fout,'hb_graphs_',reg{r},'.jpg'],'jpg');     
%      saveas(f2,[fout,'hb_bars_',reg{r},'.jpg'],'jpg');
%      saveas(f3,[fout,'hb_pvals_',reg{r},'.jpg'],'jpg');

end




%% 

box_cc13 = [box_data13{2}(1),box_data13{3}(1),    box_data13{2}(3),box_data13{3}(3),     box_data13{2}(5),box_data13{3}(5)];
figure(1); clf; BoxDist(box_cc13,names([1,1,3,3,5,5]),'missing',{'Boundary','Central'}); title('cycle 13'); xlim([0,.5]);
figure(1);  set(gcf,'color','w');
  P = zeros(G);
    for j = 1:G
        for k=1:G
        P(j,k) = ranksum(box_cc13{j},box_cc13{k});
        end
    end
    figure(3); clf;  imagesc(P); colormap(jet); colorbar; colordef white; set(gca,'color','w');

box_cc12 = [box_data12{2}(1),box_data12{3}(1),    box_data12{2}(3),box_data12{3}(3),     box_data12{2}(5),box_data12{3}(5)];
figure(2); clf; BoxDist(box_cc12,names([1,1,3,3,5,5]),'missing',{'Boundary','Central'}); title('cycle 12');   xlim([0,.5]);

  P = zeros(G);
    for j = 1:G
        for k=1:G
        P(j,k) = ranksum(box_cc12{j},box_cc12{k});
        end
    end
     figure(4); clf; imagesc(P); colormap(jet); colorbar; colordef white; set(gca,'color','w');


% 
% 
% 
% 
% 
% 
% %%
% 
% xlab = 'fraction of missed nuclei';
% F = 12; % FontSize; 
% ymax = .02; pts = 1000;
% plot_miss = cell(1,G); 
%  for k=1:G;     plot_miss{k} = miss_rate3{k}(cc14{k}); end
% plot_miss{1} = [plot_miss{1}; plot_miss{6}];  
% plot_miss{2} = [plot_miss{2}; plot_miss{9}]; 
% plot_miss{3} = [plot_miss{3}; plot_miss{4}]; 
% plot_miss{6} = [];plot_miss{9} = [];plot_miss{4} = [];
% 
% figure(35); clf;  subplot(3,1,1);
% colordef white; set(gcf,'color','w');
% 
% x = linspace(0,1,20);  % range and number of bins for histogram
% xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
%  method = 'pcubic'; % method for interpolation
% sigma = .05;  % smoothing factor for interpolation
% CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
% ylim([0,.5*ymax]); title('Boundary Early cc14');
% 
% bdata14 = plot_miss;
% 
% %~~~~ Plot Fraction of missing nuclei distributions  ~~~~
% 
% plot_miss = cell(1,G); 
%  for k=1:G;     plot_miss{k} = miss_rate3{k}(cc13{k}); end
% plot_miss{1} = [plot_miss{1}; plot_miss{6}];  
% plot_miss{2} = [plot_miss{2}; plot_miss{9}]; 
% plot_miss{3} = [plot_miss{3}; plot_miss{4}]; 
% plot_miss{6} = [];plot_miss{9} = [];plot_miss{4} = [];
% figure(35); subplot(3,1,2);
% colordef white; set(gcf,'color','w');
% x = linspace(0,1,20);  % range and number of bins for histogram
% xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
%  method = 'pcubic'; % method for interpolation
% sigma = .1;  % smoothing factor for interpolation
% CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
% ylim([0,ymax]); title('Boundary cc13');
% 
% p2enh_v_noS = ranksum(plot_miss{1},plot_miss{2})
% 
% p_noshadow_22v30 = ranksum(plot_miss{3},plot_miss{2})
% p_noprimary_22v30 = ranksum(plot_miss{5},plot_miss{8})
% p_2enhancer_22v30 = ranksum(plot_miss{1},plot_miss{7})
% 
% p_22_shadow_v_primary = ranksum(plot_miss{2},plot_miss{5})
% p_30_shadow_v_primary = ranksum(plot_miss{3},plot_miss{8})
% 
% bdata13 = plot_miss;
% 
% %~~~~ Plot Fraction of missing nuclei distributions ~~~~~
% 
% plot_miss = cell(1,G); 
%  for k=1:G;     plot_miss{k} = miss_rate3{k}(cc12{k}|cc11{k} ); end
% plot_miss{1} = [plot_miss{1}; plot_miss{6}];  
% plot_miss{2} = [plot_miss{2}; plot_miss{9}]; 
% plot_miss{3} = [plot_miss{3}; plot_miss{4}]; 
% plot_miss{6} = [];plot_miss{9} = [];plot_miss{4} = [];
%  figure(35); subplot(3,1,3);
% colordef white; set(gcf,'color','w');
% x = linspace(0,1,20);  % range and number of bins for histogram
% xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
%  method = 'pcubic'; % method for interpolation
% sigma = .1;  % smoothing factor for interpolation
% CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
% ylim([0,.5*ymax]);  title('Boundary cc11 & 12');
% 
% bdata12 = plot_miss;
% 
% %%
% 
% xlab = 'fraction of missed nuclei';
% F = 12; % FontSize; 
% ymax = .02; pts = 1000;
% plot_miss = cell(1,G); 
%  for k=1:G;     plot_miss{k} = miss_rate2{k}(cc14{k}); end
% plot_miss{1} = [plot_miss{1}; plot_miss{6}];  
% plot_miss{2} = [plot_miss{2}; plot_miss{9}]; 
% plot_miss{3} = [plot_miss{3}; plot_miss{4}]; 
% plot_miss{6} = [];plot_miss{9} = [];plot_miss{4} = [];
% 
% figure(36); clf;  subplot(3,1,1);
% colordef white; set(gcf,'color','w');
% x = linspace(0,1,20);  % range and number of bins for histogram
% xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
%  method = 'pcubic'; % method for interpolation
% sigma = .1;  % smoothing factor for interpolation
% CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
% ylim([0,.5*ymax]); title('central Early cc14');
% 
% 
% p_noshadow22_cvb = ranksum(bdata14{2},plot_miss{2})
% 
% p_noshadow30_cvb = ranksum(bdata14{3},plot_miss{3})
% 
% p_noprimary22_cvb = ranksum(bdata14{5},plot_miss{5})
% 
% p_noprimary30_cvb = ranksum(bdata14{8},plot_miss{8})
% 
% 
% %~~~~ Plot Fraction of missing nuclei distributions  ~~~~
% 
% plot_miss = cell(1,G); 
%  for k=1:G;     plot_miss{k} = miss_rate2{k}(cc13{k}); end
% plot_miss{1} = [plot_miss{1}; plot_miss{6}];  
% plot_miss{2} = [plot_miss{2}; plot_miss{9}]; 
% plot_miss{3} = [plot_miss{3}; plot_miss{4}]; 
% plot_miss{6} = [];plot_miss{9} = [];plot_miss{4} = [];
% figure(36); subplot(3,1,2);
% colordef white; set(gcf,'color','w');
% x = linspace(0,1,20);  % range and number of bins for histogram
% xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
%  method = 'pcubic'; % method for interpolation
% sigma = .1;  % smoothing factor for interpolation
% CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
% ylim([0,ymax]); title('central cc13');
% 
% p_noshadow_22v30 = ranksum(plot_miss{3},plot_miss{2})
% p_noprimary_22v30 = ranksum(plot_miss{5},plot_miss{8})
% p_2enhancer_22v30 = ranksum(plot_miss{1},plot_miss{7})
% 
% p_22_shadow_v_primary = ranksum(plot_miss{2},plot_miss{5})
% p_30_shadow_v_primary = ranksum(plot_miss{3},plot_miss{8})
% 
% 
% p_noshadow22_cvb = ranksum(bdata13{2},plot_miss{2})
% 
% p_noshadow30_cvb = ranksum(bdata13{3},plot_miss{3})
% 
% p_noprimary22_cvb = ranksum(bdata13{5},plot_miss{5})
% 
% p_noprimary30_cvb = ranksum(bdata13{8},plot_miss{8})
% 
% %~~~~ Plot Fraction of missing nuclei distributions ~~~~~
% 
% plot_miss = cell(1,G); 
%  for k=1:G;     plot_miss{k} = miss_rate2{k}(cc12{k}|cc11{k} ); end
% plot_miss{1} = [plot_miss{1}; plot_miss{6}];  
% plot_miss{2} = [plot_miss{2}; plot_miss{9}]; 
% plot_miss{3} = [plot_miss{3}; plot_miss{4}]; 
% plot_miss{6} = [];plot_miss{9} = [];plot_miss{4} = [];
%  figure(36); subplot(3,1,3);
% colordef white; set(gcf,'color','w');
% x = linspace(0,1,20);  % range and number of bins for histogram
% xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
%  method = 'pcubic'; % method for interpolation
% sigma = .1;  % smoothing factor for interpolation
% CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
% ylim([0,.5*ymax]);  title('central cc11 & 12');
% 
% 
% p_noshadow22_cvb = ranksum(bdata12{2},plot_miss{2})
% 
% p_noshadow30_cvb = ranksum(bdata12{3},plot_miss{3})
% 
% p_noprimary22_cvb = ranksum(bdata12{5},plot_miss{5})
% 
% p_noprimary30_cvb = ranksum(bdata12{8},plot_miss{8})

