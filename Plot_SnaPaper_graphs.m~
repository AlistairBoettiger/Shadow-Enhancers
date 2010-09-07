clear all; load snail_shadow_data;
[miss_cnt,miss_rate,nd,lowon] = merge_data(11,12,N,miss_cnt,miss_rate,nd,lowon);
[miss_cnt,miss_rate,nd,lowon] = merge_data(9,10,N,miss_cnt,miss_rate,nd,lowon);
[miss_cnt,miss_rate,nd,lowon] = merge_data(5,8,N,miss_cnt,miss_rate,nd,lowon);
[miss_cnt,miss_rate,nd,lowon] = merge_data(1,7,N,miss_cnt,miss_rate,nd,lowon);

ND = cell2mat(nd); 
figure(1); clf;
plot( log2(sort(ND(:))) ,'r.');
%%
cc14 =cell(1,G); cc13 = cell(1,G); cc12 = cell(1,G); 
for z=1:G
    cc14{z} = log2(ND(:,z))>8;
    cc13{z} = log2(ND(:,z))<8 & log2(ND(:,z))>7;
    cc12{z} = log2(ND(:,z))<7;
end


%% Plot Fraction of missing nuclei distributions

xlab = 'fraction of missed nuclei';


plot_miss = cell(1,G); 
for k=1:G;     plot_miss{k} = miss_rate{k}(cc14{k}); end
%for k=1:G;     plot_miss{k} = miss_rate{k}; end


 figure(1); clf;
% colordef black; set(gcf,'color','k');
colordef white; set(gcf,'color','w');

x = linspace(0,1,10);  % range and number of bins for histogram
xx = linspace(0,1,100); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation


%% CompDist fxn code 
G = 6;
data = plot_miss(1:G);
X = cell(1,G); XX = cell(1,G); M = cell(1,G); S = cell(1,G);


for k=1:G
    X{k} = x;
    XX{k} =xx;
    M{k} = method;
    S{k} = sigma;
end


% names = {'no primary, {\color{red}30C}';
%     'no primary, {\color{blue}22C}';
%     'no shadow, {\color{red}30C}';
%     'no shadow, {\color{blue}22C}'; 
%     '2 enhancers, {\color{red}30C}';
%     '2 enhancers, {\color{blue}22C}';
%     'no primary, dl^6';
%     '2 enhancers, dl^6'};


var_hist = cellfun(@nonzeros,data,'UniformOutput',false);
var_hist = cellfun(@hist,var_hist,X,'UniformOutput',false);
N_var = cellfun(@sum,var_hist,'UniformOutput',false);
dist_var = cellfun(@hist2dist,var_hist,X,XX,M,S,'UniformOutput',false);
% C = [1 0 0 ;
%     1 .5 0;
%     1 0 1;
%     .5 0 1;
%     0 1 0;
%    .5 1 0];
C = [1 0 0 ;
    1 0 0;
    .7 0 1;
    .7 0 1;
    0 1 0;
    0 1 0];

S = {'-','--','-','--','-','--'};

    figure(1), clf; set(gcf,'color','w'); 
legend_labels = cell(1,G);
for k=1:G
    legend_labels{k} = [names{k}, ' N=', num2str(N_var{k} )]; 
    plot(xx,dist_var{k},'color',C(k,:),'LineWidth',3,'LineStyle',S{k}); hold on;
end
legend(legend_labels);
xlabel(xlab,'FontSize',16); ylabel('normalized frequency','FontSize',16);
set(gca,'FontSize',15);

%%
%% CompDist fxn code 
G = 4;
trks = [2,6,7,8];
data = plot_miss(trks);
X = cell(1,G); XX = cell(1,G); M = cell(1,G); S = cell(1,G);


for k=1:G
    X{k} = x;
    XX{k} =xx;
    M{k} = method;
    S{k} = sigma;
end





var_hist = cellfun(@nonzeros,data,'UniformOutput',false);
var_hist = cellfun(@hist,var_hist,X,'UniformOutput',false);
N_var = cellfun(@sum,var_hist,'UniformOutput',false);
dist_var = cellfun(@hist2dist,var_hist,X,XX,M,S,'UniformOutput',false);
C = [1 0 0 ;
    0 1 0;
    0 0 1;
    0 0 1];
S = {'--','--','-','--'};

    figure(2), clf; set(gcf,'color','w'); 
legend_labels = cell(1,G);
for k=1:G
    g = trks(k); 
    legend_labels{k} = [names{g}, ' N=', num2str(N_var{k} )]; 
    plot(xx,dist_var{k},'color',C(k,:),'LineWidth',3,'LineStyle',S{k}); hold on;
end
legend(legend_labels);
xlabel(xlab,'FontSize',16); ylabel('normalized frequency','FontSize',16);
set(gca,'FontSize',15);

