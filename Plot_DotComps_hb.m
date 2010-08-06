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
s29_root = 'C33_29C_LacZ_hb';
sp29_root = 'hb2enh_29C_LacZ_hb';
sp22_root = 'hb2enh_22C_LacZ_hb';
s22_root = 'C33_22C_LacZ_hb';
p22_root = 'hbP_22C_LacZ_hb';
p29_root = 'C55_29C_LacZ_hb';

N = 40;

p29_miss_cnt = zeros(N,1);
p29_miss_rate = zeros(N,1); 
p29_lowon = zeros(N,1);
p29_nd = zeros(N,1);

p22_miss_cnt = zeros(N,1);
p22_miss_rate = zeros(N,1); 
p22_lowon = zeros(N,1);
p22_nd = zeros(N,1);

s29_miss_cnt = zeros(N,1);
s29_miss_rate = zeros(N,1); 
s29_lowon = zeros(N,1);
s29_nd = zeros(N,1); 

sp29_miss_cnt = zeros(N,1);
sp29_miss_rate = zeros(N,1); 
sp29_lowon = zeros(N,1);
sp29_nd = zeros(N,1);

sp22_miss_cnt = zeros(N,1);
sp22_miss_rate = zeros(N,1); 
sp22_lowon = zeros(N,1);
sp22_nd = zeros(N,1);

s22_miss_cnt = zeros(N,1);
s22_miss_rate = zeros(N,1); 
s22_lowon = zeros(N,1);
s22_nd = zeros(N,1);

xmin =200; xmax =800; ymin=150; ymax =300;

for n=1:N

    if n<10
        emb = ['0',num2str(n)];
    else
        emb = num2str(n);
    end
   
    try
    load([folder,'/',s29_root,emb,'_data.mat']); 
    
    % get the indices of all nuclei in green that are not also red.  
    % require these nuclei also fall in the 'region' for red nuclei.  
   % s29_miss_cnt(n) =  length(intersect(setdiff(pts2,pts1), ptr_nucin2));
   s29_miss_cnt(n) =  length(intersect(intersect(pts2,ptr_nucin1), ptr_nucin2)); 
    s29_miss_rate(n) = s29_miss_cnt(n)/length(pts2); 
    s29_lowon(n) = lowon_fxn(H,handles,nin2);
    s29_nd(n) = NucDensity(cent,[xmin,xmax,ymin,ymax]);

    catch ME
        disp(['can not find file' folder,'/',s29_root,emb,'_data.mat']);
    end
    
    try   
    load([folder,'/',sp29_root,emb,'_data.mat']); 
    % sp29_miss_cnt(n) =  length(intersect(setdiff(pts2,pts1), ptr_nucin2)); 
    sp29_miss_cnt(n) = length(intersect(intersect(pts2,ptr_nucin1), ptr_nucin2)); 
    sp29_miss_rate(n) = sp29_miss_cnt(n)/length(pts2); 
    sp29_lowon(n) = lowon_fxn(H,handles,nin2); 
    sp29_nd(n) = NucDensity(cent,[xmin,xmax,ymin,ymax]);
    catch ME
          disp(['can not find file' folder,'/',sp29_root,emb,'_data.mat']);
    end
    
        
    try   
    load([folder,'/',sp22_root,emb,'_data.mat']); 
    % sp22_miss_cnt(n) =  length(intersect(setdiff(pts2,pts1), ptr_nucin2)); 
    sp22_miss_cnt(n) = length(intersect(intersect(pts2,ptr_nucin1), ptr_nucin2)); 
    sp22_miss_rate(n) = sp22_miss_cnt(n)/length(pts2); 
    sp22_lowon(n) = lowon_fxn(H,handles,nin2); 
    sp22_nd(n) = NucDensity(cent,[xmin,xmax,ymin,ymax]);
    
     figure(2); clf;
    imshow(handles.It); hold on;
    
      plot(bndrys1{1}(:,1),bndrys1{1}(:,2),'r','linewidth',3 ); 
    for k=1:length(bndrys2)
     plot(bndrys2{k}(:,1),bndrys2{k}(:,2),'g','linewidth',3 );
    end
    pause(1); 
    
    catch ME
          disp(['can not find file' folder,'/',sp22_root,emb,'_data.mat']);
    end
    
    try   
    load([folder,'/',s22_root,emb,'_data.mat']); 
    s22_miss_cnt(n)= length(intersect(intersect(pts2,ptr_nucin1), ptr_nucin2)); %  =  length(intersect(setdiff(pts2,pts1), ptr_nucin2)); 
    s22_miss_rate(n) = sp22_miss_cnt(n)/length(pts2); 
    s22_lowon(n) = lowon_fxn(H,handles,nin2); 
    s22_nd(n) = NucDensity(cent,[xmin,xmax,ymin,ymax]);
    catch ME
          disp(['can not find file' folder,'/',s22_root,emb,'_data.mat']);
    end
        
    try   
    load([folder,'/',p22_root,emb,'_data.mat']); 
    p22_miss_cnt(n)= length(intersect(intersect(pts2,ptr_nucin1), ptr_nucin2)); %  =  length(intersect(setdiff(pts2,pts1), ptr_nucin2)); 
    p22_miss_rate(n) = p22_miss_cnt(n)/length(pts2); 
    p22_lowon(n) = lowon_fxn(H,handles,nin2);
    p22_nd(n) = NucDensity(cent,[xmin,xmax,ymin,ymax]);
    catch ME
          disp(['can not find file' folder,'/',p22_root,emb,'_data.mat']);
    end
    
    try   
    load([folder,'/',p29_root,emb,'_data.mat']); 
    p29_miss_cnt(n)= length(intersect(intersect(pts2,ptr_nucin1), ptr_nucin2));  % =  length(intersect(setdiff(pts2,pts1), ptr_nucin2)); 
    p29_miss_rate(n) = p29_miss_cnt(n)/length(pts2); 
    p29_lowon(n) = lowon_fxn(H,handles,nin2); 
    p29_nd(n) = NucDensity(cent,[xmin,xmax,ymin,ymax]);
    catch ME
          disp(['can not find file' folder,'/',p29_root,emb,'_data.mat']);
    end
end

%%
figure(1); clf;
plot( log2(sort([p22_nd;p29_nd;s29_nd;sp22_nd;sp29_nd])) ,'r.');
figure(1); clf;
plot( log2(sort([p22_nd, p29_nd, s29_nd, sp22_nd, sp29_nd])),'.');

p22_cc14 = log2(p22_nd)>8; p22_cc13 = log2(p22_nd)<8 & log2(p22_nd)>7; p22_cc12 = log2(p22_nd) < 7;
p29_cc14 = log2(p29_nd)>8; p29_cc13 = log2(p29_nd)<8 & log2(p29_nd)>7; p29_cc12 = log2(p29_nd) < 7; 

sp22_cc14 = log2(sp22_nd)>8; sp22_cc13 = log2(sp22_nd)<8 & log2(sp22_nd)>7; sp22_cc12 = log2(sp22_nd) < 7;
sp29_cc14 = log2(sp29_nd)>8; sp29_cc13 = log2(sp29_nd)<8 & log2(sp29_nd)>7; sp29_cc12 = log2(sp29_nd) < 7; 

s22_cc14 = log2(s22_nd)>8; s22_cc13 = log2(s22_nd)<8 & log2(s22_nd)>7; s22_cc12 = log2(s22_nd) < 7;
s29_cc14 = log2(s29_nd)>8; s29_cc13 = log2(s29_nd)<8 & log2(s29_nd)>7; s29_cc12 = log2(s29_nd) < 7;
%%

Tvar_data = {p29_lowon; p22_lowon; s29_lowon; s22_lowon; sp29_lowon; sp22_lowon};
G = length(Tvar_data);
 x = linspace(0,1,25);
  xx = linspace(0,1,100); 
 method = 'pcubic';
sigma = .1; 

X = cell(G,1); XX = cell(G,1); M = cell(G,1); S = cell(G,1); 


for k=1:G
    X{k} = x;
    XX{k} =xx;
    M{k} = method;
    S{k} = sigma;
end

var_hist = cellfun(@nonzeros,Tvar_data,'UniformOutput',false);
var_hist = cellfun(@hist,var_hist,'UniformOutput',false);
N_var = cellfun(@sum,var_hist,'UniformOutput',false);
dist_var = cellfun(@hist2dist,var_hist,X,XX,M,S,'UniformOutput',false);
figure(1); clf; 



for k=1:G
    plot(dist_var{k},'color',C(1,:),'LineWidth',3); hold on;
end


% figure(3); clf; subplot(6,1,1); set(gcf,'color','w');
%  hist(nonzeros(sp29_lowon),x); legend('control 29C');
%  set(gca,'FontSize',15); xlim([0,max(x)]); ylabel('embryos'); xlabel('variation in total mRNA'); 
%  subplot(6,1,2); hist(nonzeros(s29_lowon),x); legend('no primary 29C');
%  set(gca,'FontSize',15); xlim([0,max(x)]); ylabel('embryos'); xlabel('variation in total mRNA'); 
%  subplot(6,1,3); hist(nonzeros(sp22_lowon),x); legend('control 22C');
%  set(gca,'FontSize',15); xlim([0,max(x)]); ylabel('embryos'); xlabel('variation in total mRNA'); 
%  subplot(6,1,4); hist(nonzeros(s22_lowon),x); legend('no primary 22C');
%  set(gca,'FontSize',15); xlim([0,max(x)]); ylabel('embryos'); xlabel('variation in total mRNA'); 
%   subplot(6,1,5); hist(nonzeros(p22_lowon),x); legend('no shadow 22C');
%  set(gca,'FontSize',15); xlim([0,max(x)]); ylabel('embryos'); xlabel('variation in total mRNA'); 
%    subplot(6,1,6); hist(nonzeros(p29_lowon),x); legend('no shadow 29C');
%  set(gca,'FontSize',15); xlim([0,max(x)]); ylabel('embryos'); xlabel('variation in total mRNA'); 
 
 sp29_var =  hist(nonzeros(sp29_lowon(sp29_cc14)),x);
 s29_var = hist(nonzeros(s29_lowon(s29_cc14)),x);
 sp22_var = hist(nonzeros(sp22_lowon(sp22_cc14)),x); 
 s22_var =  hist(nonzeros(s22_lowon(s22_cc14)),x);
 p29_var =  hist(nonzeros(p29_lowon(p29_cc14)),x);
  p22_var =  hist(nonzeros(p22_lowon(p22_cc14)),x);
%   
%    sp29_var =  hist(nonzeros(sp29_lowon(sp29_cc14 | sp29_cc13 | sp29_cc12)),x);
%  s29_var = hist(nonzeros(s29_lowon(s29_cc14 | s29_cc13 | s29_cc12)),x);
%  sp22_var = hist(nonzeros(sp22_lowon(sp22_cc14 | sp22_cc13 | sp22_cc12)),x); 
%  s22_var =  hist(nonzeros(s22_lowon(s22_cc14 | s22_cc13 | s22_cc12)),x);
%  p29_var =  hist(nonzeros(p29_lowon(p29_cc14 | p29_cc13 | p29_cc12)),x);
%   p22_var =  hist(nonzeros(p22_lowon(p22_cc14 | p22_cc13 | p22_cc12)),x);
% 
%     
%    sp29_var =  hist(nonzeros(sp29_lowon( sp29_cc13 | sp29_cc12)),x);
%  s29_var = hist(nonzeros(s29_lowon( s29_cc13 | s29_cc12)),x);
%  sp22_var = hist(nonzeros(sp22_lowon( sp22_cc13 | sp22_cc12)),x); 
%  s22_var =  hist(nonzeros(s22_lowon( s22_cc13 | s22_cc12)),x);
%  p29_var =  hist(nonzeros(p29_lowon( p29_cc13 | p29_cc12)),x);
%   p22_var =  hist(nonzeros(p22_lowon( p22_cc13 | p22_cc12)),x);
  
  
N_sp29 = sum(sp29_var);
N_s29 = sum(s29_var);
N_sp22 = sum(sp22_var); 
N_s22 = sum(s22_var);
N_p22 = sum(p22_var);
N_p29 = sum(p29_var); 



 
  gn = {['2-enh 29C N=', num2str(N_sp29)],['no primary 29C N=',num2str(N_s29)]...
      ,['2-enh 22C N=',num2str(N_sp22)],['no primary 22C N=',num2str(N_s22)],...
        ['no shadow 29C N=',num2str(N_p29)],['no shadow 22C N=',num2str(N_p22)]};

 xx = linspace(0,1,100); 
 method = 'pcubic';
sigma = .1; 
   
 sp29_var_n =  hist2dist(sp29_var,x,xx,method,sigma);
 s29_var_n = hist2dist(s29_var,x,xx,method,sigma);
 sp22_var_n = hist2dist(sp22_var,x,xx,method,sigma);
 s22_var_n =  hist2dist(s22_var,x,xx,method,sigma);
 p29_var_n =  hist2dist(p29_var,x,xx,method,sigma);
 p22_var_n =  hist2dist(p22_var,x,xx,method,sigma);
 
  figure(2); clf; 
    C = hsv(6); 
    plot(xx,sp29_var_n,'color',C(1,:),'LineWidth',3); hold on;
    plot(xx,s29_var_n,'color',C(2,:),'LineWidth',3); 
    plot(xx,sp22_var_n,'color',C(3,:),'LineWidth',3); 
    plot(xx,s22_var_n,'color',C(4,:),'LineWidth',3); 
     plot(xx,p29_var_n,'color',C(5,:),'LineWidth',3);
      plot(xx,p22_var_n,'color',C(6,:),'LineWidth',3);
    legend('2-enh 29C','no primary 29C','2-enh 22C','no primary 22C',...
        'no shadow 29C', 'no shadow 22C');
    set(gcf,'color','w'); set(gca,'FontSize',15); 
    xlabel('variance of transcript intensity'); 
    xlim([0,.7]); 
    
    V{1} = nonzeros(sp29_lowon);
    V{2} = nonzeros(s29_lowon);
    V{3} = nonzeros(sp22_lowon);
    V{4} = nonzeros(s22_lowon); 
    
    P_var = zeros(4,4); 
    for i=1:4
        for j=1:4   
   P_var(i,j) =  log10(ranksum(V{i},V{j}));
        end
    end
    gn = {'2-enh 29C','no primary 29C','2-enh 22C','no primary 22C'};
    
    figure(6); clf; imagesc(P_var); colorbar;  colormap('cool');
    set(gca,'YtickLabel', str2mat(gn{:}),'YTick',1:4,'fontsize',15,...
    'YMinorTick','on'); title('Variance in Total mRNA');
    

%%

miss_rate_data = {p29_miss_rate; p22_miss_rate;  
    s29_miss_rate; s22_miss_rate;
    sp29_miss_rate; sp22_miss_rate};

x = linspace(0,1,7);
figure(1); clf; subplot(6,1,1); set(gcf,'color','w');
 hist(nonzeros(sp29_miss_rate),x);title('control 29C','FontSize',15);% legend('control 29C');
  set(gca,'FontSize',15); xlim([0,1]); ylabel('embryos'); xlabel('fraction missing nuclei'); 
 subplot(6,1,2); hist(nonzeros(s29_miss_rate),x); title('no primary 29C','FontSize',15); % legend('no primary 29C');
  set(gca,'FontSize',15); xlim([0,1]); ylabel('embryos'); xlabel('fraction missing nuclei'); 
 subplot(6,1,3); hist(nonzeros(sp22_miss_rate),x); title('control 22C','FontSize',15); % legend('control 22C');
 set(gca,'FontSize',15); xlim([0,1]); ylabel('embryos'); xlabel('fraction missing nuclei'); 
 subplot(6,1,4); hist(nonzeros(s22_miss_rate),x); title('no primary 22C','FontSize',15);% legend('no primary 22C');
 set(gca,'FontSize',15); xlim([0,1]); ylabel('embryos'); xlabel('fraction missing nuclei'); 
  subplot(6,1,5); hist(nonzeros(p22_miss_rate),x); title('no shadow 22C','FontSize',15);% legend('no primary 22C');
 set(gca,'FontSize',15); xlim([0,1]); ylabel('embryos'); xlabel('fraction missing nuclei'); 
  subplot(6,1,6); hist(nonzeros(p29_miss_rate),x); title('no shadow 29C','FontSize',15);% legend('no primary 22C');
 set(gca,'FontSize',15); xlim([0,1]); ylabel('embryos'); xlabel('fraction missing nuclei'); 
 
 
%  sp29_miss =   hist(nonzeros(sp29_miss_rate(sp29_cc13 | sp29_cc12 )),x); 
%  s29_miss = hist(nonzeros(s29_miss_rate(s29_cc13 | s29_cc12 )),x);
%  sp22_miss =  hist(nonzeros(sp22_miss_rate(sp22_cc13 | sp22_cc12 )),x);
%  s22_miss =  hist(nonzeros(s22_miss_rate(s22_cc13 | s22_cc12 )),x);
%  p29_miss = hist(nonzeros(p29_miss_rate(p29_cc13 | p29_cc12  )),x);
%  p22_miss =  hist(nonzeros(p22_miss_rate(p22_cc13 | p22_cc12 )),x);



 

sp29_miss =   hist(nonzeros(sp29_miss_rate(sp29_cc14)),x); 
 s29_miss = hist(nonzeros(s29_miss_rate(s29_cc14)),x);
 sp22_miss =  hist(nonzeros(sp22_miss_rate(sp22_cc14)),x);
 s22_miss =  hist(nonzeros(s22_miss_rate(s22_cc14)),x);
 p29_miss = hist(nonzeros(p29_miss_rate(p29_cc14)),x);
 p22_miss =  hist(nonzeros(p22_miss_rate(p22_cc14)),x);


 N_sp29 = sum(sp29_miss);
N_s29 = sum(s29_miss);
N_sp22 = sum(sp22_miss); 
N_s22 = sum(s22_miss);
N_p22 = sum(p22_miss);
N_p29 = sum(p29_miss); 
 
 xx = linspace(0,1,100); 
 method = 'pcubic';
sigma = .01 % .15; 
   
 sp29_miss =  hist2dist(sp29_miss,x,xx,method,sigma);
 s29_miss = hist2dist(s29_miss,x,xx,method,sigma);
 sp22_miss = hist2dist(sp22_miss,x,xx,method,sigma);
 s22_miss =  hist2dist(s22_miss,x,xx,method,sigma);
 p29_miss = hist2dist(p29_miss,x,xx,method,sigma);
 p22_miss =  hist2dist(p22_miss,x,xx,method,sigma);
 

 
  gn = {['2-enh 29C N=', num2str(N_sp29)],['no primary 29C N=',num2str(N_s29)]...
      ,['2-enh 22C N=',num2str(N_sp22)],['no primary 22C N=',num2str(N_s22)],...
        ['no shadow 29C N=',num2str(N_p29)],['no shadow 22C N=',num2str(N_p22)]};
 
  figure(4); clf; set(gcf,'color','w');
    C = hsv(6); 
    plot(xx,sp29_miss,'color',C(1,:),'LineWidth',3); hold on;
    plot(xx,s29_miss,'color',C(2,:),'LineWidth',3); 
    plot(xx,sp22_miss,'color',C(3,:),'LineWidth',3); 
    plot(xx,s22_miss,'color',C(4,:),'LineWidth',3); 
    plot(xx,p29_miss,'color',C(5,:),'LineWidth',3); 
    plot(xx,p22_miss,'color',C(6,:),'LineWidth',3); 
  %  legend('2-enh 29C','no primary 29C','2-enh 22C','no primary 22C',...
  %      'no shadow 29C','no shadow 22C');
  legend(gn{:});
  
    set(gca,'FontSize',15); xlabel('fraction inactive nuclei'); 

    
    % Stats
    
     
    V{1} = nonzeros(sp29_miss_rate);
    V{2} = nonzeros(s29_miss_rate);
    V{3} = nonzeros(sp22_miss_rate);
    V{4} = nonzeros(s22_miss_rate); 
     V{5} = nonzeros(p29_miss_rate);
     V{6} = nonzeros(p22_miss_rate);
    
    
    P_dot = zeros(6,6); 
    for i=1:6
        for j=1:6   
   P_dot(i,j) =  log10(ranksum(V{i},V{j}));
        end
    end
  
    
    figure(5); clf; imagesc(P_dot); colorbar; colormap('cool');
    set(gca,'YtickLabel', str2mat(gn{:}),'YTick',1:6,'fontsize',15,...
    'YMinorTick','on'); title('dots'); 

%%


  load([folder,'/',s29_root,'03','_data.mat']); %'06'
%load([folder,'/',sp29_root,'23','_data.mat']); 

figure(1); clf; imshow(handles.It);

% I = handles.It;
% I(:,:,3) = 1.4*I(:,:,3);
% figure(1); clf; imshow(1.5*I);

ints = sort(double(handles.Im1(:)));
mV = ints(round(.996*length(handles.Im1(:))));
Im1 = imadjust(handles.Im1,[12/255,mV/255]); 

Cell_bnd = false(h,w); Cell_bnd(handles.cell_bords) = 1;
Cell_bnd = imdilate(Cell_bnd,strel('disk',1)); 

 figure(2); clf;
    Io = uint8(zeros(h,w,3));
    Io(:,:,1) = Im1; 
    Io(:,:,2) =  imadd(uint8(0*L2.*Cell_bnd),0*handles.Im2); % imadd(uint8(255*Cell_bnd),1*handles.Im2);  %
    Io(:,:,3) =   imadd(uint8(255*Cell_bnd),0*handles.In); % uint8(255*Cell_bnd); %
   % Io = uint8(bsxfun(@times,double(Io)/155,double(handles.In)));
  clf; imshow(Io); hold on;
  title(['embryo ',emb],'FontSize',15); pause(2); 
  
  figure(5); 
 Io = uint8(zeros(h,w,3));
    Io(:,:,2) =  imadd(uint8(0*L2.*Cell_bnd),1*handles.Im2); % imadd(uint8(255*Cell_bnd),1*handles.Im2);  %
    % Io = uint8(bsxfun(@times,double(Io)/155,double(handles.In)));
  clf; imshow(Io); hold on;
  title(['embryo ',emb],'FontSize',15); pause(2); 

  %%

  for n=1:N

    if n<10
        emb = ['0',num2str(n)];
    else
        emb = num2str(n);
    end
   
  
% emb =  '02', 11,12,27;  % sp22 
% emb =  '10';  % s22 
% emb =  '26';  % sp29
% emb = '06' % s29

try
load([folder,'/',s29_root,emb,'_data.mat']); 

%figure(3); clf; hist(double(handles.Im1(:)),200);     

% 

cell_int = H; 
       Nnucs = max(H(:));
        mu_I = zeros(Nnucs,1); 
%         sig_I = zeros(Nnucs,1);
         mu_N = zeros(Nnucs,1);
        reg_data = regionprops(H,'PixelIdxList');
            for i=1:Nnucs; 
                pixes = reg_data(i).PixelIdxList;      
                mu_I(i) = median(handles.Im1(pixes));
                 mu_N(i) = median(handles.In(pixes));
%                 sig_I(i) = std(double(handles.Im1(pixes)));
                cell_int(pixes) = mu_I(i); 
            end

            cell_int = uint8( 255*cell_int./max(cell_int(:))  );
           % figure(4); clf; imshow(cell_int);

%
Cell_bnd = false(h,w);
Cell_bnd(handles.cell_bords) = 1;
 Cell_bnd = imdilate(Cell_bnd,strel('disk',1)); 

Im1 = flipud(fliplr(Im1));
L1 = flipud(fliplr(L1));
L2 = flipud(fliplr(L2));
Cell_bnd = flipud(fliplr(Cell_bnd));
cell_int = flipud(fliplr(cell_int));
handles.In = flipud(fliplr(handles.In));
handles.Im2 = flipud(fliplr(handles.Im2));

 figure(3); clf;
    Io = uint8(zeros(h,w,3));
    Io(:,:,1) = imadd(cell_int,Im1); 
    Io(:,:,2) =  imadd(uint8(155*L2.*Cell_bnd),0*handles.Im2); % imadd(uint8(255*Cell_bnd),1*handles.Im2);  %
    Io(:,:,3) =   imadd(uint8(0*L2.*Cell_bnd),handles.In); % uint8(255*Cell_bnd); %
   % Io = uint8(bsxfun(@times,double(Io)/155,double(handles.In)));
  clf; imshow(Io); hold on;
  title(['embryo ',emb],'FontSize',15); pause(2); 
  
  
catch
end

  end
  
  %%
  
  load([folder,'/',sp29_root,'26','_data.mat']);  % s29_12  sp29 26


  figure(1); clf; imshow(handles.It);
  
cell_int = H; 
       Nnucs = max(H(:));
        mu_I = zeros(Nnucs,1); 
        reg_data = regionprops(H,'PixelIdxList');
            for i=1:Nnucs; 
                pixes = reg_data(i).PixelIdxList;      
                mu_I(i) = median(handles.Im1(pixes));
                cell_int(pixes) = mu_I(i); 
            end
            cell_int = uint8( 255*cell_int./max(cell_int(:)));
            

  
  s29_DI = uint8(bsxfun(@times,double(cell_int)/255,double(handles.Im1)));
   C = jet(70);  C = [zeros(3,3);C;ones(300,3)];
  figure(4); clf; imagesc(s29_DI); colormap(C);

  
    s29_DI = uint8(bsxfun(@times,double(cell_int)/255,double(handles.In)));
   C = jet(200);  C = [zeros(3,3);C;ones(50,3)];
  figure(5); clf; imagesc(s29_DI); colormap(C);
  
%%
x = linspace(0,1,10);
figure(1); clf; 
 hist(nonzeros(sp29_miss_rate),x); 
hblock =  findobj(gca,'type','patch');
set(hblock,'FaceColor','r','EdgeColor','k');
hold on; 
hist(nonzeros(s29_miss_rate),x); alpha(.6);
xlim([min(x)-.05,max(x)]); 
set(gca,'FontSize',17);
xlabel('Fraction missing reporter'); 
set(gcf,'Color','w');
legend('sog2enh, 29C','sogS, 29C');

x = linspace(0,1,10);
figure(2); clf; 
 hist(nonzeros(sp29_miss_rate),x); 
hblock =  findobj(gca,'type','patch');
set(hblock,'FaceColor','r','EdgeColor','k');
hold on; 
hist(nonzeros(sp22_miss_rate),x); alpha(.6);
xlim([min(x)-.05,max(x)]); 
set(gca,'FontSize',17);
xlabel('Fraction missing reporter'); 
set(gcf,'Color','w');
legend('sog2enh, 29C','sog2enh, 22C');

%%
 


 figure(2); clf;
     Io = uint8(zeros(h,w,3));
    Io(:,:,1) = uint8(255*L1);
    Io(:,:,2) = uint8(255*L2);
    Io(:,:,3) = 30*handles.In;
     DI = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
  clf; imshow(DI);
    for k=1:handles.Rs1
 hold on; plot(bndrys1{k}(:,1),bndrys1{k}(:,2),'r','LineWidth',3);
    end
    for k=1:handles.Rs2
 hold on; plot(bndrys2{k}(:,1),bndrys2{k}(:,2),'g','LineWidth',3);
    end
 hold on;
 plot(rna_x1,rna_y1,'r.');
 plot(rna_x2,rna_y2,'g.');

         Nnucs_on1 = length(pts1)
         Nnucs_on2 = length(pts2)
        on1_not_2 = setdiff(pts1,pts2);
        on2_not_1 = setdiff(pts2,pts1);
        cnt_1not2 = length(on1_not_2);
        cnt_2not1 = length(on2_not_1);