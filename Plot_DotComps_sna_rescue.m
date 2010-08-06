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

folder =  '/Volumes/Data/Lab Data/Shadow_data';

emb_roots = {'MP12-21p_22C_LacZ_vn_sna';
    'MP7-5Ap_22C_LacZ_vn_sna';
    'MP7-5Ap_29C_LacZ_vn_sna';
    'MP7-5IIp_22C_LacZ_vn_sna';
    'MP7-5IIp_29C_LacZ_vn_sna'};


% 12-21 01~, 06, 08~ 10~ 12 homozygous snail
%  hets 09, 02~ 03~, 05 07

% 7-5IIp 03 04 06 08  homo sna
% hets 05 09 07

N = 15;
K = length(emb_roots); 

ectopic_cnt = cell(1,K); 
ectopic_rate = cell(1,K); 
nd = cell(1,K); 

for z=1:K
    ectopic_cnt{z} = zeros(N,1);
    ectopic_rate{z} = zeros(N,1); 
    nd{z} = zeros(N,1);
end


xmin =200; xmax =800; ymin=150; ymax =300;
for z=1:K % k=2;
    for n=1:N
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
        ectopic_cnt{z}(n) =  length(intersect(unique(pts1), ptr_nucin2)); 
       ectopic_rate{z}(n) =ectopic_cnt{z}(n)/length(ptr_nucin2); 
        nd{z}(n) = NucDensity(cent,[xmin,xmax,ymin,ymax]);

        catch ME
            disp(['can not find file' folder,'/',emb_roots{z},emb,'_data.mat']);
        end
    end
end

%%
ND = cell2mat(nd); 
figure(1); clf;
plot( log2(sort(ND(:))) ,'r.');
%%
cc14 =cell(1,K); cc13 = cell(1,K); cc12 = cell(1,K); 
for z=1:K
    cc14{z} = log2(ND(:,z))>8;
    cc13{z} = log2(ND(:,z))<8 & log2(ND(:,z))>7;
    cc12{z} = log2(ND(:,z))<7;
end

%%

 x = linspace(0,1,28);
  xx = linspace(0,1,100); 
 method = 'pcubic';
sigma = .1; 

X = cell(1,K); XX = cell(1,K); M = cell(1,K); S = cell(1,K);


for k=1:K
    X{k} = x;
    XX{k} =xx;
    M{k} = method;
    S{k} = sigma;
end

var_hist = cellfun(@nonzeros,ectopic_rate,'UniformOutput',false);
var_hist = cellfun(@hist,var_hist,X,'UniformOutput',false);
N_var = cellfun(@sum,var_hist,'UniformOutput',false);
dist_var = cellfun(@hist2dist,var_hist,X,XX,M,S,'UniformOutput',false);
figure(1); clf; 

C = hsv(K);

for k=1:K
    plot(dist_var{k},'color',C(k,:),'LineWidth',3); hold on;
end
legend(emb_roots);

