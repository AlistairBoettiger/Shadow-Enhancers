%%                          Plot_DotComps2_ConcDept_sna.m                %%
%
% Analyzing New Shadow data
%
%
% Alistair Boettiger                                   Date Begun: 03/05/10
% Levine Lab                                     Functional Since: 03/06/10
%                                                   Last Modified: 01/12/11

%% Description
% comparison
%
%
%% Updates
% Revised 01/12/11 to use most recent  formulation of age structure and
% plotting tools.  
% Really should make this into a function
% Revised 03/17/11  


%% Source Code
clear all;

folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed';

emb_roots = {'MP05_29C_y_sna'; % 1
    'MP05_22C_y_sna';          % 2
    'MP06xYW_30C_y_sna';       % 3
    'MP06xYW_22C_y_sna';       % 4
    'MP10_29C_y_sna';          % 5
    'MP10_22C_y_sna';          % 6
    'MP05xYW_30C_sna_y-full';  % 7
    'MP10xYW_30C_sna_y-full';  % 8
    'MP05xdl6_25C_pt1';        % 9
    'MP05xdl6_25C_pt2';        % 10
    'MP10xdl6_25C_pt1';        % 11 
    'MP10xdl6_25C_pt2'};       % 12
          

 names = {'2 enh 22C';
          '2 enh 30C';
          'no proximal 22C';
          'no proximal 30C';
          'no distal 22C';
          'no distal 30C';
          '2 enh dl6';
          'no proximal dl6'};
 
N = 70;
K = length(emb_roots); 
G= length(names);

layer_off_rate = cell(K,N); 
layer_off_cnt = cell(K,N);
nd = cell(1,K); 

for z=1:K
    nd{z} = zeros(N,1);
end


xmin = .2; xmax = .8; ymin = .25; ymax = .75;
% as fractions of the original image dimensions.  

for z=1:K % k=2;
    for n=1:N
        if n<10
            emb = ['0',num2str(n)];
        else
            emb = num2str(n);
        end

        try
        load([folder,'/',emb_roots{z},emb,'_data.mat']);   
        

        
        if length(pts2) > 700; % only analyze substantial snail regions  
            
            % Setup first border
            Regd = imdilate(Reg1,strel('disk',5))-Reg1; 
            borderNuc = nonzeros(unique(Regd.*H));
            Border{1} = ismember(H,borderNuc);
           %  figure(1); clf; imagesc(Border{1}); 
            [h,w] = size(Reg1); 
            border_ind{1} = borderNuc; 

            % loop through successive internal borders
            j = 1; 
            while length(borderNuc)>20

                Bplot = 0;
                for k=1:j
                    Bplot = Bplot+ k*Border{k} ;
                end
              % figure(1); clf; imagesc(Bplot);  colormap jet; colorbar; caxis([0,10]); pause(.01); 

                j = j+1;
                B = imdilate(Border{j-1},strel('disk',5)) - Border{j-1};
                prevBorder = Reg1;
                for k=1:j
                    blnk = false(h,w); 
                    blnk(Bplot>1) = 1;
                    prevBorder = Reg1 - blnk; % figure(2); clf; imagesc(prevBorder); 
                end

                borderNuc = nonzeros(unique(B.*H.*prevBorder ));
                Border{j} = ismember(H,borderNuc);
                border_ind{j} = borderNuc; 

            end

            J = j-1; 
            layer_off_cnt{z,n} = zeros(1,J); 
            layer_size = zeros(1,J); 
            for j=1:J
                layer_off_cnt{z,n}(j) = length(intersect(border_ind{j},setdiff(pts2,pts1)));  % total number of missing nuclei in each layer
                layer_size(j) = length(intersect(border_ind{j},pts2));
            end
                layer_off_rate{z,n} = layer_off_cnt{z,n}./layer_size;
        end 
        
       

          
           if length(H) > 2000
               im_dim = 2048;
           else
               im_dim = 1024;
           end
               
           lims =  round([xmin,xmax,ymin,ymax]*im_dim); 
           nd{z}(n) = NucDensity(cent,lims,0);

        catch ME
            disp(ME.message); 
            %disp(['can not find file' folder,'/',emb_roots{z},emb,'_data.mat']);
        end
    end
end
data_folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Code_Data/';

 %  save([data_folder,'snail_ConcDept']);


%%
k = 0;

profile = zeros(N*K,11); 

for z=1:K
    for n=1:N
        if isempty(layer_off_rate(z,n)) == 0
            k = k+1;
            profile(k,1:length(layer_off_rate{z,n})) = layer_off_rate{z,n}/mean(layer_off_rate{z,n});
        end
    end
end

% find in which row of cells the maximum miss fraction occured.
profile(isnan(profile)) = 0; 
Hp = profile; Hp(Hp<1) = 0;  Hp(isnan(Hp)) = 0; 
Hp = logical(Hp);
[x,y] = find(Hp);

max_misses = hist(y,1:16); % histogram of number of times each bin was recorded as the max
norm_bins = sum(logical(10000*profile)); % normalize against the number of embryos that contain that bin (not all embryos have 10 bins).  


figure(1); clf; colordef black;

bar(max_misses./norm_bins*max(norm_bins)); xlim([0,10]);
ylabel('Frequency: max misses at distance x','FontSize',14);
xlabel('x, distance (in cells) from edge of expression','FontSize',14); 
xlim([0,10.5]); set(gcf,'color','k'); set(gca,'FontSize',14);

figure(2); clf; colordef black; 
C = jet(length(profile)); 
for c=1:length(profile)
  plot(linspace(1,16,16),profile(c,:),'.','color',C(c,:)); hold on;
end
ylabel('layer miss rate / mean miss rate','FontSize',14);
xlabel('distance from edge of expression','FontSize',14); 
xlim([0,11]); set(gcf,'color','k'); set(gca,'FontSize',14);


% %% 
% % clear all;  load snail_SD_012111
% 
% foff{1} = miss_rate{6};                 Nnuc{1} = nd{6}; % 2 enh 22 C, MP10
% foff{2} = [miss_rate{5},miss_rate{8}];  Nnuc{2} = [nd{5},nd{8}]; % 2 enh 30C MP10
% foff{3} = miss_rate{2};                 Nnuc{3} = nd{2}; % MP05 22C 
% foff{4} = [miss_rate{1}; miss_rate{7}]; Nnuc{4} = [nd{1}, nd{7}]; % MP05 30C
% foff{5} = miss_rate{4};                 Nnuc{5} = nd{4};  % MP06 22C
% foff{6} = miss_rate{3};                 Nnuc{6} = nd{3};  % MP06 30C
% foff{7} = [miss_rate{11},miss_rate{12}]; Nnuc{7} = [nd{11}, nd{12}];% MP10 dl6 25C
% foff{8} = [miss_rate{9},miss_rate{10}]; Nnuc{8} = [nd{9},nd{10}]; % MP05 dl6 25C
% 


%%
clear all;

data_folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Code_Data/';
load([data_folder,'snail_ConcDept']);
%%



foff = cell(8,1); 

 foff{1} = layer_off_rate(6,:) ;   % 2 enh 22 C, MP10
 foff{2} = [layer_off_rate(5,:),layer_off_rate(8,:) ]; %  2 enh 30C MP10
 foff{3} = layer_off_rate(2,:);   % % MP05 22C 
 foff{4} = [layer_off_rate(1,:),layer_off_rate(7,:)  ];  % MP05 30C
 foff{5} = layer_off_rate(4,:) ; % MP06 22C
 foff{6} = layer_off_rate(3,:); %  MP06 30C
 foff{7} = [layer_off_rate(11,:),layer_off_rate(12,:)];% MP10 dl6 25C
 foff{8} = [layer_off_rate(9,:),layer_off_rate(10,:)]; % MP05 dl6 25C

 
Nnucs = cell(8,1); 
 Ls = size(layer_off_rate,2)
 Nnuc = cell(1,6); 
 Nnuc{1} = [nd{6}, zeros(1,Ls-length(nd{6}))] ;  % 
 Nnuc{2} = [nd{5}, zeros(1,Ls-length(nd{5}));  nd{8}, zeros(1,Ls-length(nd{8})) ]; % 
 Nnuc{3} = [nd{2}, zeros(1,Ls-length(nd{2}))]; % 
 Nnuc{4} = [nd{1}, zeros(1,Ls-length(nd{1})); nd{7}, zeros(1,Ls-length(nd{7}))]; %
 Nnuc{5} = [nd{4}, zeros(1,Ls-length(nd{4}))]; % 
 Nnuc{6} = [nd{3}, zeros(1,Ls-length(nd{3}))]; % 
 Nnuc{7} = [nd{11}, zeros(1,Ls-length(nd{11})); nd{12}, zeros(1,Ls-length(nd{12}))]; % 
 Nnuc{8} = [nd{9}, zeros(1,Ls-length(nd{9}));  nd{10}, zeros(1,Ls-length(nd{10}))]; % 

% square off entries
for k=1:G
        foff{k} = [foff{k},cell(1,300-length(foff{k}))];
        data = Nnuc{k};
        Nnuc{k} = [data; zeros(300-length(data),1)];
 end
     



ND = cell2mat(Nnuc); 
age_offset = 5.2;

emb_cycle = age_offset + log2( nonzeros( sort(ND(:)) ) );
figure(2); clf; colordef white;  plot( emb_cycle ,'r.');


title(['hb embryos, N = ',num2str(length(nonzeros(ND(:))) )  ],'FontSize',15);
set(gca,'FontSize',15); grid on;
set(gcf,'color','w'); ylabel('log_2(nuc density)');  xlabel('embryo number'); 
ylim([10,14.99]);


cc14 =cell(1,G); cc13 = cell(1,G); cc12 = cell(1,G); cc11 = cell(1,G); cc10 = cell(1,G); cc9 = cell(1,G); 
for z=1:G
    logage =   age_offset + log2( ND(:,z) );    
    cc14{z} = logage >14;
    cc13{z} = logage <14  & logage> 13;
    cc12{z}  = logage <13 & logage > 12;
    cc11{z} = logage <12 & logage > 0 ;
end

%% Histogram


z = 8
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
xlabel('distance from sna-boundary (cells)','FontSize',F); 
title(['snail ', names{z},'  cc14   N = ',num2str(N)],'FontSize',F);
set(gca,'FontSize',F);
 xlim([1.5,10.5]); ylim([0,.6]); 
 


%% Compare hists
F = 16;

    figure(1); clf; 
    colordef black;
    set(gcf,'color','k');

      zs = [8,3]; 
  % zs = [8,3,7,1]; % no primary +/- dorsal
 %   zs = [4,3,2,1];
%    zs = [6,5,2,1];
    C = flipud(hsv(8));
    A = [1,.2]; 
  leg_lab = cell(length(zs),1);  
for k = 1:length(zs)
    z = zs(k); 
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
    h(k) = bar(miss_dist_hist(1:end),'FaceColor',C(z,:)); hold on; 
 
  %  plot(miss_dist_hist(2:end),'.','color',C(z,:),'linewidth',3); hold on;
    leg_lab{k} = ['snail ', names{z},'  cc14   N = ',num2str(N)];
end

%legend(leg_lab,'Location','SouthOutside');
legend(leg_lab);
ylabel('fraction of missing nulcei','FontSize',F); 
xlabel('distance from sna-boundary (cells)','FontSize',F); 
set(gca,'FontSize',F-4);
xlim([1.5,8.5]); 
ylim([0,.55]); 


%% Compare graphs

F = 16;

    figure(1); clf; 
    colordef black;
    set(gcf,'color','k');

      
    data = zeros(8,10); 
    C = flipud(hsv(8));
    leg_lab = cell(length(zs),1);  
  
for z = 1:8
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
    data(z,1:10) = miss_dist_hist(1:10); 
    leg_lab{z} = ['snail ', names{z},'  cc14   N = ',num2str(N)];
end

noP_temp = data(4,:) - data(3,:)
noP_dl = data(8,:) - data(3,:)
cntl_dl = data(7,:) - data(1,:)
cntl_temp = data(2,:) - data(1,:)


figure(3); clf; set(gcf,'color','k');
plot(2:8,noP_dl(2:8),'m.-'); hold on; xlim([1.5,8.5]);
% plot(2:8,cntl_dl(2:8),'b.-'); 

plot(2:8,data(3,2:8) - data(1,2:8),'g.-'); 

plot(2:8,data(4,2:8) - data(2,2:8),'y.-'); 



[b,bint,r,rint,stats] = regress((2:8)',noP_temp(2:8)')

corrcoef(noP_temp(2:8))

 
 
 
 