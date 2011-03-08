
%% Figs_Revised_GapShadows.m

clear all
fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Results/';
folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed/';


%% ------------ Analyze Missed Nuclei, hunchback ---------------------- %%
%edat =  'MP09_22C_y_hb12_data.mat'; % cc13  11
emb = '12';  eroot = 'MP09_22C_y_hb';% 11*,12,22, 25,  26 27, 28, 31, 36
MP09 = [eroot,emb,'_data.mat']; % cc13   

% eroot = 'BAC01b_22C_y_hb'; % 5,  21   22,30-cc12,  27
% eroot = 'BAC01b_30C_y_hb';  % 09  10
emb = '09';   eroot = 'MP01_22C_y_hb';% 04  09*  17    cc12: 10-13
MP01 = [eroot,emb,'_data.mat']; % cc13   

% eroot = 'BAC02_30C_y_hb'; %eroot = 'MP02_22C_y_hb'; % % eroot =
% 'MP02_30C_LacZ_hb';%
emb = '05';  
eroot = 'BAC02_22C_y_hb';% 10% = 02   20% =  03*  06   10
MP02 = [eroot,emb,'_data.mat']; % cc13   

figure(1); set(gcf,'color','k');  title('Kr'); 
subplot(3,1,1); f = [0,2]; I_MP01 = endVrept(folder,MP01,1.3,.3,f);
subplot(3,1,2); f = [1,0]; I_MP02 = endVrept(folder,MP02,3,.3,f);
subplot(3,1,3); f = [1,2]; I_MP09 = endVrept(folder,MP09,1.5,.3,f);

     imwrite(I_MP01,[fout,'MP01_comp.tif'],'tif');
     imwrite(I_MP02,[fout,'MP02_comp.tif'],'tif');
     imwrite(I_MP09,[fout,'MP09_comp.tif'],'tif');
%-------------------------------------------------------------------------%
%



%% ------------ Analyze Ectopic Nuclei, hunchback ---------------------- %%


%MP09 = 'MP09_22C_y_hb41_data.mat'; % cc14 41 44


%MP01 = 'BAC01b_22C_y_hb45_data.mat'; %   37
MP01 = 'BAC01b_30C_y_hb15_data.mat'; % 18  19**  02  ~13   15**
%MP01 = 'MP01_22C_y_hb09_data.mat'; % cc14  06 07**  08 27
MP02 = 'MP02_22C_y_hb51_data.mat'; % cc14  44*  51*  42*    52  83  89**
MP09 = 'BAC09_30C_y_hb16_data.mat'; % cc14  16* ~31 35 43
%MP09 = 'BAC09_22C_y_hb43_data.mat'; % cc14   ~31 ~35 43**



figure(1); set(gcf,'color','k');   
 subplot(3,1,1); f = [0,0]; I_MP01 = ShowEctopNuc(folder,MP01,2.4,.3,f,1);
 subplot(3,1,2); f = [1,0]; I_MP02 = ShowEctopNuc(folder,MP02,1.8,.3,f,1);
subplot(3,1,3); f = [1,0]; I_MP09 = ShowEctopNuc(folder,MP09,2.4,.3,f,1);
%%
     imwrite(I_MP01,[fout,'MP01_ectop.tif'],'tif');
     imwrite(I_MP02,[fout,'MP02_ectop.tif'],'tif');
     imwrite(I_MP09,[fout,'MP09_ectop.tif'],'tif');
%-------------------------------------------------------------------------%




%% ------------ Analyze Missed Nuclei, Kr and kni ---------------------- %%

kr1 = 'krCD1_22C_LacZ_kr08_data.mat'; % 06 08
kr2 = 'krCD2_22C_LacZ_kr21_data.mat'; %   12* , 19, 21, 22
kr2x = 'kr2enh_22C_LacZ_kr09_data.mat';% 09
kni5 = 'kni_5p_22C_LacZ_kni18_data.mat'; % 18  05 12 21  23
kniI = 'kni_int_22C_LacZ_kni15_data.mat'; % 06 22*, 23, 11, 13, 15, 18, 33
kni2x = 'kni_2enh_22C_LacZ_kni19_data.mat';% 16 19 04 06 10 22, 24


figure(1); set(gcf,'color','k');  title('Kr'); 
subplot(3,1,1); f = [1,2]; Ikr1 = endVrept(folder,kr1,1.9,.3,f);
subplot(3,1,2); f = [0,0]; Ikr2 = endVrept(folder,kr2,1.9,.3,f);
subplot(3,1,3); f = [1,2]; Ikr2x = endVrept(folder,kr2x,1.3,.3,f);

figure(2); set(gcf,'color','k');  title('kni'); 
subplot(3,1,1); f = [0,2]; Ikni5 = endVrept(folder,kni5,2.5,.3,f);
subplot(3,1,2); f = [0,2]; IkniI = endVrept(folder,kniI,2,.3,f);
subplot(3,1,3); f = [0,2]; Ikni2x = endVrept(folder,kni2x,2.5,.3,f);  
           
     
 %% write images to folder
 
 
     imwrite(Ikr1,[fout,'krCD1_comp.tif'],'tif');
     imwrite(Ikr2,[fout,'krCD2_comp.tif'],'tif');
     imwrite(Ikr2x,[fout,'kr2enh_comp.tif'],'tif');
     imwrite(Ikni5,[fout,'kni_5p_comp.tif'],'tif');
     imwrite(IkniI,[fout,'kni_int_comp.tif'],'tif');
     imwrite(Ikni2x,[fout,'kni_2enh_comp.tif'],'tif');
%-------------------------------------------------------------------------%     
     




%% ----------------- Analyze Ectopic Nuclei kni ------------------------ %%
kni5 = 'kni_5p_22C_LacZ_kni18_data.mat'; % 18* 12   %   12 21  23   % cc13   14 16  19 20   cc11=17
kniI = 'kni_int_22C_LacZ_kni07_data.mat'; % 06~  22  24  02e  04e  07e ~14 16 %   22*, 23, 11, 13, 15, 18, 33  % cc13   01  10 11  12
kni2x = 'kni_2enh_22C_LacZ_kni16_data.mat';% 16 06 19* ~22  %    % cc13 08   17
  
figure(2); set(gcf,'color','k');
subplot(3,1,1); f = [0,2]; Ikni5 = ShowEctopNuc(folder,kni5,2.9,.3,f,0);
subplot(3,1,2); f = [0,0]; IkniI = ShowEctopNuc(folder,kniI,1.8,.3,f,0);
subplot(3,1,3); f = [0,0]; Ikni2x = ShowEctopNuc(folder,kni2x,3.5,.3,f,0);

%%
     imwrite(Ikni5,[fout,'kni_5p_ectop.tif'],'tif');
     imwrite(IkniI,[fout,'kni_int_ectop.tif'],'tif');
     imwrite(Ikni2x,[fout,'kni_2enh_ectop.tif'],'tif');
%-------------------------------------------------------------------------%






%% ----------------- Analyze Ectopic Nuclei Kr ------------------------ %%
kr1 = 'krCD1_22C_LacZ_kr06_data.mat'; % 06* 08~   04  ~17  18
kr2 = 'krCD2_22C_LacZ_kr01_data.mat'; %  01* 12* , 19*, 21~, 22  13
kr2x = 'kr2enh_22C_LacZ_kr09_data.mat';% 09  08 10 14


figure(1); set(gcf,'color','k');  title('Kr'); 
subplot(3,1,1); f = [0,2]; Ikr1 = ShowEctopNuc(folder,kr1,1.9,.3,f,0);
subplot(3,1,2); f = [1,2]; Ikr2 = ShowEctopNuc(folder,kr2,1.9,.3,f,0);
subplot(3,1,3); f = [1,2]; Ikr2x = ShowEctopNuc(folder,kr2x,1.5,.3,f,0);

%%
imwrite(Ikr1,[fout,'KrCD1_ectop.tif'],'tif');
imwrite(Ikr2,[fout,'KrCD2_ectop.tif'],'tif');
imwrite(Ikr2x,[fout,'Kr_2enh_ectop.tif'],'tif');





