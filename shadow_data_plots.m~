
% Alistair Boettiger                            Date Begun: 02/06/10
%
%
% Image shadow data exported from im_nucdots_v8


% shadow_data_plots.m 


% data order 
%  [RNA_cnt1,RNA_cnt2,Nucs_cnt1,Nucs_cnt2,FAC1,FAC2,handles.nucT,cnt_1not2,
%  cnt_2not1,str2double(handles.emb)];


% Sog Primary Data at 22C, analyzed 02/05-06/10

% trimmed to include only like stages 
sogP_22C = [474,581,590,652,0.80339,0.8911,977,98,197,1
753,897,868,975,0.86751,0.92,977,92,229,2
666,780,885,877,0.75254,0.8894,977,133,245,5
813,798,992,891,0.81956,0.89562,977,163,140,6
711,737,956,824,0.74372,0.89442,977,182,205,7
944,920,1115,1007,0.84664,0.9136,977,162,143,8
603,541,762,602,0.79134,0.89867,977,179,112,10
902,761,1014,880,0.88955,0.86477,977,221,70,11
621,481,896,588,0.69308,0.81803,977,229,80,12
504,548,760,728,0.66316,0.75275,977,185,227,13
761,825,1028,999,0.74027,0.82583,977,148,208,14
668,660,883,775,0.75651,0.85161,977,151,148,17
642,904,879,964,0.73038,0.93776,977,58,324,18
615,646,888,847,0.69257,0.76269,977,196,217,21
704,697,895,798,0.78659,0.87343,977,178,168,22];

% Sog Primary Data at 30C, analyzed 02/05/10
% trimmed to include only like stages 
sogP_30C = [464,659,668,780,0.69461,0.84487,977,76,263,3
397,791,799,930,0.49687,0.85054,977,48,444,6
515,767,767,837,0.67145,0.91637,977,40,296,7
540,725,801,846,0.67416,0.85697,977,73,262,8
442,546,794,664,0.55668,0.82229,977,120,227,10
678,887,893,1094,0.75924,0.81079,977,105,317,11
585,412,1037,545,0.56413,0.75596,977,321,145,12
543,510,758,634,0.71636,0.80442,977,166,133,13
740,955,1020,1165,0.72549,0.81974,977,128,341,15
349,499,678,631,0.51475,0.79081,977,91,239,17
431,601,686,742,0.62828,0.80997,977,72,240,19
899,855,1072,1073,0.83862,0.79683,977,240,194,22
744,807,911,927,0.81668,0.87055,977,115,177,23
813,880,1038,1001,0.78324,0.87912,977,144,211,24
670,838,914,943,0.73304,0.88865,977,86,255,25
518,756,744,847,0.69624,0.89256,977,49,286,26
592,740,795,836,0.74465,0.88517,977,74,219,28
541,684,940,878,0.57553,0.77904,977,158,296,30
754,919,976,1056,0.77254,0.87027,977,96,260,31
562,743,751,842,0.74834,0.88242,977,75,251,33
395,449,663,614,0.59578,0.73127,977,134,184,35];


% sogP staged
sogP22C_10cell = [621,481,896,588,0.69308,0.81803,977,229,80,12
                504,548,760,728,0.66316,0.75275,977,185,227,13
                458,471,696,591,0.65805,0.79695,977,159,176,16
                615,646,888,847,0.69257,0.76269,977,196,217,21];

sogP22C_15cell = [753,897,868,975,0.86751,0.92,977,92,229,2
                666,780,885,877,0.75254,0.8894,977,133,245,5
                813,798,992,891,0.81956,0.89562,977,163,140,6
                711,737,956,824,0.74372,0.89442,977,182,205,7
                603,541,762,602,0.79134,0.89867,977,179,112,10
                902,761,1014,880,0.88955,0.86477,977,221,70,11
                761,825,1028,999,0.74027,0.82583,977,148,208,14
                668,660,883,775,0.75651,0.85161,977,151,148,17
                642,904,879,964,0.73038,0.93776,977,58,324,18
                615,646,888,847,0.69257,0.76269,977,196,217,21
                704,697,895,798,0.78659,0.87343,977,178,168,22];

sogP22C_20cell = [944,920,1115,1007,0.84664,0.9136,977,162,143,8];


% sogS staged
sogP30C_10cell = [
471,337,858,509,0.54895,0.66208,977,259,121,21
442,546,794,664,0.55668,0.82229,977,120,227,10
349,499,678,631,0.51475,0.79081,977,91,239,17
390,382,733,543,0.53206,0.7035,977,182,170,18
431,601,686,742,0.62828,0.80997,977,72,240,19
283,466,629,619,0.44992,0.75283,977,107,298,32
];

sogP30C_15cell = [464,659,668,780,0.69461,0.84487,977,76,263,3
397,791,799,930,0.49687,0.85054,977,48,444,6
540,725,801,846,0.67416,0.85697,977,73,262,8
678,887,893,1094,0.75924,0.81079,977,105,317,11
744,807,911,927,0.81668,0.87055,977,115,177,23
813,880,1038,1001,0.78324,0.87912,977,144,211,24
670,838,914,943,0.73304,0.88865,977,86,255,25
518,756,744,847,0.69624,0.89256,977,49,286,26
592,740,795,836,0.74465,0.88517,977,74,219,28
1029,689,1244,1023,0.82717,0.67351,977,426,86,29
541,684,940,878,0.57553,0.77904,977,158,296,30
754,919,976,1056,0.77254,0.87027,977,96,260,31
562,743,751,842,0.74834,0.88242,977,75,251,33
306,240,345,268,0.88696,0.89552,977,92,31,34
395,449,663,614,0.59578,0.73127,977,134,184,35];

sogP30C_20cell = [515,767,767,837,0.67145,0.91637,977,40,296,7
740,955,1020,1165,0.72549,0.81974,977,128,341,15
899,855,1072,1073,0.83862,0.79683,977,240,194,22
670,838,914,943,0.73304,0.88865,977,86,255,25
754,919,976,1056,0.77254,0.87027,977,96,260,31];

x = linspace(0,1,10); 
figure(2); clf; subplot(4,1,1); hist(sogP_22C(:,6),x); title('22C sog'); xlim([0,1]); xlabel('FAC');
subplot(4,1,2); hist(sogP_22C(:,5),x); title('22C LacZ'); xlim([0,1]); xlabel('FAC');
subplot(4,1,3); hist(sogP_30C(:,6),x); title('30C sog'); xlim([0,1]); xlabel('FAC');
subplot(4,1,4); hist(sogP_30C(:,5),x); title('30C LacZ'); xlim([0,1]); xlabel('FAC');



ten30 = sogP30C_10cell(:,9)./sogP30C_10cell(:,4); 
ten22 = sogP22C_10cell(:,9)./sogP22C_10cell(:,4);
fifteen30 = sogP30C_15cell(:,9)./sogP30C_15cell(:,4); 
fifteen22 = sogP22C_15cell(:,9)./sogP22C_15cell(:,4); 
twenty30 = sogP30C_20cell(:,9)./sogP30C_20cell(:,4); 
twenty22 = sogP22C_20cell(:,9)./sogP22C_20cell(:,4); 

figure(2); clf; colordef white; set(gcf,'color','w'); 
subplot(3,1,1); hist(ten30,x); hold on;
h =  findobj(gca,'type','patch'); set(h,'FaceColor','r','EdgeColor','w');
subplot(3,1,1); hist(ten22,x); title('[#sogP off]/[# sog] 10cell width');
alpha(.6); legend('30C','22C'); xlim([0,1]);
subplot(3,1,2); hist(fifteen30,x);  hold on;
h =  findobj(gca,'type','patch'); set(h,'FaceColor','r','EdgeColor','w');
subplot(3,1,2); hist(fifteen22,x); title('[#sogP off]/[# sog] 15cell width');
alpha(.6);legend('30C','22C'); xlim([0,1]);
subplot(3,1,3); hist(twenty30,x);  hold on;
h =  findobj(gca,'type','patch'); set(h,'FaceColor','r','EdgeColor','w');
subplot(3,1,3); hist(twenty22,x); title('[#sogP off]/[# sog] 20cell width');
alpha(.6); legend('30C','22C'); xlim([0,1]);


%%


ten30 = sogP30C_10cell(:,5)./sogP30C_10cell(:,6); 
ten22 = sogP22C_10cell(:,5)./sogP22C_10cell(:,6);
fifteen30 = sogP30C_15cell(:,5)./sogP30C_15cell(:,6); 
fifteen22 = sogP22C_15cell(:,5)./sogP22C_15cell(:,6); 
twenty30 = sogP30C_20cell(:,5)./sogP30C_20cell(:,6); 
twenty22 = sogP22C_20cell(:,5)./sogP22C_20cell(:,6); 

x = linspace(0,1,14); 

figure(2); clf; colordef white; set(gcf,'color','w'); 
subplot(3,1,1); hist(ten30,x); hold on;
h =  findobj(gca,'type','patch'); set(h,'FaceColor','r','EdgeColor','w');
subplot(3,1,1); hist(ten22,x); title('[FAC sogP]/[FAC sog] 10cell width');
alpha(.6); legend('30C','22C','Location','NorthWest'); xlim([0,1]);
subplot(3,1,2); hist(fifteen30,x);  hold on;
h =  findobj(gca,'type','patch'); set(h,'FaceColor','r','EdgeColor','w');
subplot(3,1,2); hist(fifteen22,x); title('[FAC sogP]/[FAC sog] 15cell width');
alpha(.6);legend('30C','22C','Location','NorthWest'); xlim([0,1]);
subplot(3,1,3); hist(twenty30,x);  hold on;
h =  findobj(gca,'type','patch'); set(h,'FaceColor','r','EdgeColor','w');
subplot(3,1,3); hist(twenty22,x); title('[FAC sogP]/[FAC sog] 20cell width');
alpha(.6); legend('30C','22C','Location','NorthWest'); xlim([0,1]);


%% 

  A =  sogP_22C(:,5)./sogP_22C(:,6); 
  B =  sogP_30C(:,5)./sogP_30C(:,6);
 x = linspace(0,2,25); 
% 
% 
% A =  sogP_30C(:,5); 
%  B = sogP_22C(:,6);

%A = sogP_22C(:,9)./sogP_22C(:,4);  % # sog no LacZ / # sog
% B = sogP_30C(:,9)./sogP_30C(:,4); 
%x = linspace(0,500,10);  % number of endogenous on with reporter off
x = linspace(0,1,15); 
% 

figure(3); clf; colordef black;  set(gcf,'color','k')
hist(A,x); hold on;
h =  findobj(gca,'type','patch'); set(h,'FaceColor','r','EdgeColor','w'); 

hist(B,x) 
alpha(.6);

xlim([0,1.2]);
xlabel('[FAC sogP]/[FAC sog]','FontSize',15);
ylabel('Frequency','FontSize',15); 
set(gca,'FontSize',15);
legend('SogP 22C','SogP 30C','Location','NorthWest');

%%  Normalize Histogram

An = hist(A,x); An = An./sum(An);
Bn = hist(B,x); Bn = Bn./sum(Bn);

figure(5); clf; bar(x,An);  
h =  findobj(gca,'type','patch');
set(h,'FaceColor','r','EdgeColor','w'); 
colordef black; set(gcf,'color','k');  hold on; 

bar(x,Bn); 
alpha(.6);

xlim([0,.6]);
xlabel('Fraction of Unactivated Nuclei','FontSize',15);
ylabel('Frequency','FontSize',15); 
set(gca,'FontSize',15);
legend('SogP 22C','SogP 30C','Location','NorthEast');


%% binomial fits 

% doesn't fit binomial.  Variance is too large for the number of total
% cells participating.  
Na = 4; Nb = 2;

pA = mean(Na*A)/Na
var(Na*A)
vA = Na*pA*(1-pA)
mA = Na*pA

pB = mean(Nb*B)/Nb
var(Nb*B)
vB = Nb*pB*(1-pB)
mB = Nb*pB

%%



% gene_names = {'P-YW:End ratio'; 'P-dl:End ratio'; 'S-YW: End ratio';
%    'S-dl:End ratio';'P YW 22 domain';'P dl 25 domain'; 'S YW 22 domain';'S dl 25 domain'};

%---


%     %------------
%        Data = {
%         sort(SogPxYW_22C(:,6));
%         sort(SogPxdl_25C(:,6));
%         sort( SogPxYW_22C(:,5)./SogPxYW_22C(:,6) ); 
%         sort( SogPxdl_25C(:,5)./SogPxdl_25C(:,6) );
%         sort( SogSxYW_22C(:,6) );
%         sort( SogSxYW_22C(:,5)./SogSxYW_22C(:,6) );        
%         };
%     
%     
% gene_names = {'Endgonous YW P';  'Endogenous dl';
%     'P-YW:End ratio'; 'P-dl:End ratio';'End YW S';'S-YW: End ratio' };
%     %--------
    
        G = length(Data); 
  
% rescale 
embs = sum(cellfun('size', Data,1)) ;
L0 = 100;
DA = zeros(G,L0-1);
for g=1:G
    for i=1:L0-1
        k= 1+floor(length(Data{g})*i/L0);
        DA(g,i) = Data{g}(k); 
   end
end
 sat = 1;% 1.4;
 DA(DA>sat) = sat; % saturate to account for different strength signals that merge due to lack of resolution 

 
 
% sort by asynchrony score
asyncT = 1;% 1; % .8
Dt = DA; Dt(DA>asyncT)=NaN;
sv = nanmean(Dt,2);
DAs = sortrows([sv,DA]); 
gs = sortrows([sv,(1:G)'])  ;
gn = gene_names(gs(:,2)); 
gene_names2 = [num2cell(sv),gene_names];

 D = DAs(:,2:end);

 C = colormap(1-hot(256));
 C = [ones(round(.3*256),3); C; zeros(round(.1*256),3)];
 
 F=18;
% 
 figure(2); clf; imagesc(D); figure(2);
 colorbar;  colordef white;
  set(gcf,'color','w');
set(gca,'YtickLabel',...
 str2mat(gn{:}),...
     'YTick',[1:G],'fontsize',3/4*F,...
    'YMinorTick','on');
xlabel('percent of embryos','Fontsize',F);
 colormap(C);
colordef black; set(gcf,'color','k');
