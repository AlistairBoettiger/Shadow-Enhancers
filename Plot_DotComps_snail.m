%%                          Plot_DotComps_sna.m                          %%
%
% Analyzing New Shadow data
%
%
% Alistair Boettiger                                   Date Begun: 03/05/10
% Levine Lab                                     Functional Since: 03/05/10
%                                                   Last Modified: 06/25/10 

%% Description
% comparison
%
%
%% Updates


%% Source Code
clear all;

load sna_emb_annot;  % 
% load notes called eg. s22_age, recording the cell cycle of each embryo 
% and s22_emb, which has a 1 for all embryos we want to analyze and a 0 for
% numbers that should be skipped.  



folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed';
s29_rootB = 'MP05xYW_30C_sna_y-full';
s29_root = 'MP05_29C_y_sna';
s22_root = 'MP05_22C_y_sna';
sp29_root = 'MP10_29C_y_sna';
sp22_root = 'MP10_22C_y_sna';


N = 60;
s29_age = 14*ones(N,1); s29_emb = ones(N,1); 

s29_miss_cnt = zeros(N,1);
s29_miss_rate = zeros(N,1); 
s29_lowon = zeros(N,1);

s29_miss_cntB = zeros(N,1);
s29_miss_rateB = zeros(N,1); 
s29_lowonB = zeros(N,1);


sp29_miss_cnt = zeros(N,1);
sp29_miss_rate = zeros(N,1); 
sp29_lowon = zeros(N,1);

sp22_miss_cnt = zeros(N,1);
sp22_miss_rate = zeros(N,1);
sp22_lowon = zeros(N,1);

s22_miss_cnt = zeros(N,1);
s22_miss_rate = zeros(N,1); 
s22_lowon = zeros(N,1);

for n=1:N

    if n<10
        emb = ['0',num2str(n)];
    else
        emb = num2str(n);
    end
   
    try
    load([folder,'/',s29_root,emb,'_data.mat']); 
    
        if s29_age(n)>12 && s29_emb(n)==1
            % get the indices of all nuclei in green that are not also red.  
            % require these nuclei also fall in the 'region' for red nuclei.  
            s29_miss_cnt(n) =  length(intersect(setdiff(pts2,pts1), ptr_nucin2));  
            s29_miss_rate(n) = s29_miss_cnt(n)/length(pts2); 
            s29_lowon(n) = lowon_fxn(H,handles,nin2,ptr_nucin2,emb);
        end
    catch ME
               disp(ME.message);
    end
    
    
    try
    load([folder,'/',s29_rootB,emb,'_data.mat']); 
    
        if s29_age(n)>12 && s29_emb(n)==1
            % get the indices of all nuclei in green that are not also red.  
            % require these nuclei also fall in the 'region' for red nuclei.  
            s29_miss_cntB(n) =  length(intersect(setdiff(pts2,pts1), ptr_nucin2));  
            s29_miss_rateB(n) = s29_miss_cnt(n)/length(pts2); 
            s29_lowonB(n) = lowon_fxn(H,handles,nin2,ptr_nucin2,emb);
        end
    catch ME
               disp(ME.message);
    end
    
%     
    try   
    load([folder,'/',sp29_root,emb,'_data.mat']); 
   
        if sp29_age(n)>12 && sp29_emb(n)==1
            sp29_miss_cnt(n) =  length(intersect(setdiff(pts2,pts1), ptr_nucin2)); 
            sp29_miss_rate(n) = sp29_miss_cnt(n)/length(pts2); 
            sp29_lowon(n) = lowon_fxn(H,handles,nin2,ptr_nucin2,emb); 
        end
    catch ME
               disp(ME.message);
     end
    
    
        
    try   
    load([folder,'/',sp22_root,emb,'_data.mat']);    
        if sp22_age(n)>12 && sp22_emb(n)==1
            sp22_miss_cnt(n) =  length(intersect(setdiff(pts2,pts1), ptr_nucin2)); 
            sp22_miss_rate(n) = sp22_miss_cnt(n)/length(pts2); 
            sp22_lowon(n) = lowon_fxn(H,handles,nin2,ptr_nucin2,emb);  
        end
    catch ME
               disp(ME.message);
    end
      
    
    try   
    load([folder,'/',s22_root,emb,'_data.mat']); 
        if s22_age(n)>12 && s22_emb(n)==1
            s22_miss_cnt(n) =  length(intersect(setdiff(pts2,pts1), ptr_nucin2)); 
            s22_miss_rate(n) = s22_miss_cnt(n)/length(pts2); 
            s22_lowon(n) = lowon_fxn(H,handles,nin2,ptr_nucin2,emb);  
        end
    catch ME
        disp(ME.message);  
    end
    
        
end


s29_miss_cnt = [s29_miss_cnt; s29_miss_cntB];
s29_miss_rate = [s29_miss_rate; s29_miss_rateB];
s29_lowon = [s29_lowon ; s29_lowonB ];


%%

x = linspace(0,1,14);
figure(3); clf; subplot(4,1,1); set(gcf,'color','w');
 hist(nonzeros(sp29_lowon),x); legend('control 30C');
  set(gca,'FontSize',15); xlim([0,1]); ylabel('embryos'); xlabel('fraction missing nuclei'); 
 subplot(4,1,2); hist(nonzeros(s29_lowon),x); legend('no primary 30C');
  set(gca,'FontSize',15); xlim([0,1]); ylabel('embryos'); xlabel('fraction missing nuclei'); 
 subplot(4,1,3); hist(nonzeros(sp22_lowon),x); legend('control 22C');
 set(gca,'FontSize',15); xlim([0,1]); ylabel('embryos'); xlabel('fraction missing nuclei'); 
 subplot(4,1,4); hist(nonzeros(s22_lowon),x); legend('no primary 22C');
 set(gca,'FontSize',15); xlim([0,1]); ylabel('embryos'); xlabel('fraction missing nuclei'); 
%%

  sp29_miss =   hist(nonzeros(sp29_miss_rate),x); 
 s29_miss = hist(nonzeros(s29_miss_rate),x);
 sp22_miss =  hist(nonzeros(sp22_miss_rate),x);
 s22_miss =  hist(nonzeros(s22_miss_rate),x);

 xx = linspace(0,1,100); 
 method = 'pcubic';
sigma = .3;%.15; 
   
 sp29_miss =  hist2dist(sp29_miss,x,xx,method,sigma);
 s29_miss = hist2dist(s29_miss,x,xx,method,sigma);
 sp22_miss = hist2dist(sp22_miss,x,xx,method,sigma);
 s22_miss =  hist2dist(s22_miss,x,xx,method,sigma);
 
  figure(4); clf; set(gcf,'color','w');
    C = hsv(4); 
    plot(xx,sp29_miss,'color',C(1,:),'LineWidth',3); hold on;
    plot(xx,s29_miss,'color',C(2,:),'LineWidth',3); 
    plot(xx,sp22_miss,'color',C(3,:),'LineWidth',3); 
    plot(xx,s22_miss,'color',C(4,:),'LineWidth',3); 
    legend('2-enh 29C','no primary 29C','2-enh 22C','no primary 22C');
    set(gca,'FontSize',15); xlabel('fraction inactive nuclei'); 

    
    % Stats
    
     
    V{1} = nonzeros(sp29_miss_rate);
    V{2} = nonzeros(s29_miss_rate);
    V{3} = nonzeros(sp22_miss_rate);
    V{4} = nonzeros(s22_miss_rate); 
    
    P_dot = zeros(4,4); 
    for i=1:4
        for j=1:4   
   P_dot(i,j) =  log10(ranksum(V{i},V{j}));
        end
    end
    gn = {'2-enh 30C','no primary 30C','2-enh 22C','no primary 22C'};
    
    figure(5); clf; imagesc(P_dot); colorbar; colormap('cool');
    set(gca,'YtickLabel', str2mat(gn{:}),'YTick',1:4,'fontsize',15,...
    'YMinorTick','on'); title('dots'); 

%%
x = linspace(0,1,15);
figure(3); clf; subplot(4,1,1); set(gcf,'color','w');
 hist(nonzeros(sp29_miss_rate),x); legend('control 29C');
  set(gca,'FontSize',15); xlim([0,1]); ylabel('embryos'); xlabel('fraction missing nuclei'); 
 subplot(4,1,2); hist(nonzeros(s29_miss_rate),x); legend('no primary 29C');
  set(gca,'FontSize',15); xlim([0,1]); ylabel('embryos'); xlabel('fraction missing nuclei'); 
 subplot(4,1,3); hist(nonzeros(sp22_miss_rate),x); legend('control 22C');
 set(gca,'FontSize',15); xlim([0,1]); ylabel('embryos'); xlabel('fraction missing nuclei'); 
 subplot(4,1,4); hist(nonzeros(s22_miss_rate),x); legend('no primary 22C');
 set(gca,'FontSize',15); xlim([0,1]); ylabel('embryos'); xlabel('fraction missing nuclei'); 
 
  mp10_29 = hist(nonzeros(sp29_miss_rate),x); 
 mp05_29 =  hist(nonzeros(s29_miss_rate),x); 
 mp10_22 =  hist(nonzeros(sp22_miss_rate),x); 
 mp05_22 = hist(nonzeros(s22_miss_rate),x); 
 
 col = hsv(4); 
 
  
  % normalize
      n_mp10_29 = mp10_29./sum(mp10_29); 
      n_mp05_29 = mp05_29./sum(mp05_29); 
      n_mp10_22 = mp10_22./sum(mp10_22); 
      n_mp05_22 = mp05_22./sum(mp05_22); 
      
 figure(2); clf;
 plot(x,n_mp10_29,'.','color',col(1,:)); hold on;
  plot(x,n_mp05_29,'.','color',col(2,:));
  plot(x,n_mp10_22,'.','color',col(3,:));
  plot(x,n_mp05_22,'.','color',col(4,:));
%  
    [estimates, fit_mp10_29] = fxn_fit_mixedexp(x,n_mp10_29,[1,-10,-1,-20]);
    [estimates, fit_mp05_29] = fxn_fit_mixedexp(x,n_mp05_29,[1,-10,-1,-20]);
    [estimates, fit_mp10_22] = fxn_fit_mixedexp(x,n_mp10_22,[1,-10,-1,-20]);
    [estimates, fit_mp05_22] = fxn_fit_mixedexp(x,n_mp05_22,[1,-10,-1,-20]);
  figure(2); 
%    plot(x,1*exp(-10*x)-1*exp(-20*x),'r');
 
 hold on;
 plot(x,fit_mp10_29,'color',col(1,:),'linewidth',3); 
  plot(x,fit_mp05_29,'color',col(2,:),'linewidth',3); 
  plot(x,fit_mp10_22,'color',col(3,:),'linewidth',3); 
   plot(x,fit_mp05_22,'color',col(4,:),'linewidth',3); 
 
   
   legend('MP10 30C','MP05 30C', 'MP10 22C', 'MP05 22C')
   
   sum(fit_mp10_29)
   sum(fit_mp05_29)
   sum(fit_mp10_22)
   sum(fit_mp05_22)


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
legend('MP10 30C','MP05, 30C');



x = linspace(0,1,10);
figure(2); clf; 
 hist(nonzeros(s22_miss_rate),x); 
hblock =  findobj(gca,'type','patch');
set(hblock,'FaceColor','r','EdgeColor','k');
hold on; 
hist(nonzeros(s29_miss_rate),x); alpha(.6);
xlim([min(x)-.05,max(x)]); 
set(gca,'FontSize',17);
xlabel('Fraction missing reporter'); 
set(gcf,'Color','w');
legend('MP05 22C','MP05, 30C');


x = linspace(0,1,10);
figure(3); clf; 
 hist(nonzeros(s22_miss_rate),x); 
hblock =  findobj(gca,'type','patch');
set(hblock,'FaceColor','r','EdgeColor','k');
hold on; 
hist(nonzeros(sp22_miss_rate),x); alpha(.6);
xlim([min(x)-.05,max(x)]); 
set(gca,'FontSize',17);
xlabel('Fraction missing reporter'); 
set(gcf,'Color','w');
legend('MP05 22C','MP10, 22C');

% %%  Detailed Views of choice embryos
% 
% age = zeros(1,N);
% good_emb = zeros(1,N); 
% 
% for n=1:N % n=23
% 
%     if n<10
%         emb = ['0',num2str(n)];
%     else
%         emb = num2str(n);
%     end
%    
%     try
% load([folder,'/',s22_root,emb,'_data.mat']); % 02, 08, 14, 15 18 19 22 23 26   
%  
% Cell_bnd = false(h,w);
% Cell_bnd(cellbords) = 1;
% Cell_bnd = imdilate(Cell_bnd,strel('disk',3)); 
% 
% 
% 
% Nuc_line = H; Nuc_line(Nuc_line>0)=1; Nuc_line(Cell_bnd == 1)=0;
% Nuc_cntr = H; Nuc_cntr(Nuc_cntr>0)=1; Nuc_cntr(imdilate(Cell_bnd,strel('disk',1))==1)=0;
% Cell_bnd = logical(Nuc_line-Nuc_cntr); 
% 
%  figure(2); clf;
%     Io = uint8(zeros(h,w,3));
%     Io(:,:,1) = imadd(uint8(0*L2n1a.*Cell_bnd),3*handles.Im1); 
%     Io(:,:,2) =  imadd(uint8(255*L2n1.*Cell_bnd),.8*handles.Im2); % imadd(uint8(255*Cell_bnd),1*handles.Im2);  %
%     Io(:,:,3) =   imadd(uint8(255*L2n1.*Cell_bnd),.7*handles.In); % uint8(255*Cell_bnd); %
%     % DI = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
%   clf; imshow(Io); hold on;
% 
%    % plot(rna_x1,rna_y1,'r.','MarkerSize',5);
%  % plot(rna_x2,rna_y2,'g.','MarkerSize',5);
% 
%  title(['embryo', emb],'FontSize',15); 
%  
%  age(n) = input('embryo age');
%  good_emb(n) = input('good image?'); 
% 
%  
%  
% %      
% %       figure(3); clf;
% %     Io = uint8(zeros(h,w,3));
% %     Io(:,:,1) = imadd(uint8(0*L2n1a.*Cell_bnd),3*handles.Im1); 
% %     Io(:,:,2) =  imadd(uint8(0*L2n1.*Cell_bnd),.8*handles.Im2); % imadd(uint8(255*Cell_bnd),1*handles.Im2);  %
% %     Io(:,:,3) =   imadd(uint8(0*L2n1.*Cell_bnd),.7*handles.In); % uint8(255*Cell_bnd); %
% %     % DI = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
% %   clf; imshow(Io); hold on;
% % 
% %  % plot(rna_x1,rna_y1,'r.','MarkerSize',5);
% %  % plot(rna_x2,rna_y2,'g.','MarkerSize',5);
% %      title('disambiguation of mRNA 1'); 
%     catch
%     
%     end
% end
% 
% % save sna_s22 age good_emb
% 
%  %%    Detailed views of choice embryos
%      %  % s29 02, 08, 14, 15 18 19 22 23 26   
%      load([folder,'/',s29_root,'15','_data.mat']); %sp22 14 08 28 07  12 15 17 18 23 25 26 27 29 31 34
%  
% Cell_bnd = false(h,w);
% Cell_bnd(cellbords) = 1;
% Cell_bnd = imdilate(Cell_bnd,strel('disk',3)); 
%  
% 
% 
% Nuc_line = H; Nuc_line(Nuc_line>0)=1; Nuc_line(Cell_bnd == 1)=0;
% Nuc_cntr = H; Nuc_cntr(Nuc_cntr>0)=1; Nuc_cntr(imdilate(Cell_bnd,strel('disk',1))==1)=0;
% Cell_bnd = logical(Nuc_line-Nuc_cntr); 
% 
%  figure(2); clf;
%     Io = uint8(zeros(h,w,3));
%     Io(:,:,1) = imadd(uint8(0*L2n1a.*Cell_bnd),2.2*handles.Im1); 
%     Io(:,:,2) =  imadd(uint8(255*L2n1.*Cell_bnd),.9*handles.Im2); % imadd(uint8(255*Cell_bnd),1*handles.Im2);  %
%     Io(:,:,3) =   imadd(uint8(255*L2n1.*Cell_bnd),.9*handles.In); % uint8(255*Cell_bnd); %
%     % DI = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
%   clf; imshow(Io); hold on;
% 
%   % plot(rna_x1,rna_y1,'r.','MarkerSize',5);
%  % plot(rna_x2,rna_y2,'g.','MarkerSize',5);
%      title('disambiguation of mRNA 1'); 
%      
%      
%       figure(3); clf;
%     Io = uint8(zeros(h,w,3));
%     Io(:,:,1) = imadd(uint8(0*L2n1a.*Cell_bnd),2.6*handles.Im1); 
%     Io(:,:,2) =  imadd(uint8(0*L2n1.*Cell_bnd),.8*handles.Im2); % imadd(uint8(255*Cell_bnd),1*handles.Im2);  %
%     Io(:,:,3) =   imadd(uint8(0*L2n1.*Cell_bnd),.9*handles.In); % uint8(255*Cell_bnd); %
%     % DI = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
%   clf; imshow(Io); hold on;
% % 
% %  % plot(rna_x1,rna_y1,'r.','MarkerSize',5);
% %  % plot(rna_x2,rna_y2,'g.','MarkerSize',5);
% %      title('disambiguation of mRNA 1'); 
% 
% %%
% 
% x = linspace(0,1,10);
% figure(2); clf; 
%  hist(nonzeros(sp29_miss_rate),x); 
% hblock =  findobj(gca,'type','patch');
% set(hblock,'FaceColor','r','EdgeColor','k');
% hold on; 
% hist(nonzeros(sp22_miss_rate),x); alpha(.6);
% xlim([min(x)-.05,max(x)]); 
% set(gca,'FontSize',17);
% xlabel('Fraction missing reporter'); 
% set(gcf,'Color','w');
% legend('sog2enh, 29C','sog2enh, 22C');
% 
% %%
%  
% 
% 
%  figure(2); clf;
%      Io = uint8(zeros(h,w,3));
%     Io(:,:,1) = uint8(255*L1);
%     Io(:,:,2) = uint8(255*L2);
%     Io(:,:,3) = 30*handles.In;
%      DI = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
%   clf; imshow(DI);
%     for k=1:handles.Rs1
%  hold on; plot(bndrys1{k}(:,1),bndrys1{k}(:,2),'r','LineWidth',3);
%     end
%     for k=1:handles.Rs2
%  hold on; plot(bndrys2{k}(:,1),bndrys2{k}(:,2),'g','LineWidth',3);
%     end
%  hold on;
%  plot(rna_x1,rna_y1,'r.');
%  plot(rna_x2,rna_y2,'g.');
% 
%          Nnucs_on1 = length(pts1)
%          Nnucs_on2 = length(pts2)
%         on1_not_2 = setdiff(pts1,pts2);
%         on2_not_1 = setdiff(pts2,pts1);
%         cnt_1not2 = length(on1_not_2);
%         cnt_2not1 = length(on2_not_1);