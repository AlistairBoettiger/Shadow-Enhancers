%%                          Plot_DotComps.m                             %%
%
% Analyzing New Shadow data
%
%
% Alistair Boettiger                                   Date Begun: 03/05/10
% Levine Lab                                     Functional Since: 03/06/10
%                                                   Last Modified: 03/17/10

%% Description
% Nuclear density for all nuclei
%
%
%% Updates


%% Source Code
clear all;
folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed';

emb_roots = {'MP09_22C_y_hb';
            'MP02_22C_y_hb';
            'MP02_30C_y_hb';
            'MP02_30C_LacZ_hb';   % yes it's actually yellow
            'MP01_22C_y_hb';
            'MP05_29C_y_sna'; % 1
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
            'MP10xdl6_25C_pt2';
            'C33_29C_LacZ_hb';
            'hb2enh_29C_LacZ_hb';
            'hb2enh_22C_LacZ_hb';
            'C33_22C_LacZ_hb';
            'hbP_22C_LacZ_hb';
            'C55_29C_LacZ_hb';
            'kr2enh_22C_LacZ_kr';  
            'krCD1_22C_LacZ_kr';
            'krCD2_22C_LacZ_kr'
            };  
          


     

N = 90;
K = length(emb_roots); 

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

          [hW,hL] = size(H);
           if hW*hL > 1E6;
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


%  save hb_SD-9-13-10
 %  save hb_shadow_yellow_data;   % load hb_shadow_yellow_data;

save all_nuc_data;

%%
 % clear all; load snail_shadow_data;

ND = cell2mat(nd); 

emb_cycle = 4.25 + log2( nonzeros( sort(ND(:)) ) );
figure(4); clf; plot( emb_cycle ,'r.');

title(['all embryos, N = ',num2str(length(emb_cycle) )  ],'FontSize',15);
set(gca,'FontSize',15); grid on;

figure(2); clf;
plot(nonzeros(sort(ND(:))),'r.');


NDi = ND(:);

ND_lab = [NDi,(1:length(NDi))'];

ND_sort = sortrows(ND_lab);

%%
fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Results/';

nuc_folder = {'cc10','cc11','cc12','cc13','cc14'}; 

 for j = find(ND_sort(:,1)>1,1):K*N;  % j=N*18+7;
%for j = K*N-50:K*N; 
ez = 1+floor(ND_sort(j,2)/N);  
n = ND_sort(j,2) - (ez-1)*N;

% ez = 18; n=4;
        if n<10
            emb = ['0',num2str(n)];
        else
            emb = num2str(n);
        end
        clear handles;
        load([folder,'/',emb_roots{ez},emb,'_data.mat']);   
        
% figure(5); clf; imshow(handles.It(:,:,3)); 
% title(['cycle = ',num2str(4.2 + log2( ND_sort(j,1) ) )  ],'FontSize',15);
% pause(.1);

cc = round( 4.2 + log2( ND_sort(j,1) ))

Io = imresize(handles.It(:,:,3),.5);
imwrite(Io,[fout,nuc_folder{cc-9},'/','im_',emb_roots{ez},emb,'.jpg']);


         [hW,hL] = size(H);
           if hW*hL > 1E6;
               im_dim = 2048;
           else
               im_dim = 1024;
           end
        
       lims =  round([xmin,xmax,ymin,ymax]*im_dim); 
           ndi = NucDensity(cent,lims,1)


              FigOut = figure(1); 
              FigOut = imresize(FigOut,.5); 
        saveas(FigOut,[fout,nuc_folder{cc-9},'/','nd_',emb_roots{ez},...
          emb,'.jpg']);   
           
end


%%
