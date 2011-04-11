
%% Figs_Revised_GapShadows.m

clear all
fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Results/';
folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed/';


%% ------------ Analyze Missed Nuclei, sog 30C ---------------------- %%


sP_30 =  'SogP30C_LacZ_sog15_data.mat';  %  07 11 'sogP_C81_22C_LacZ_sog';
sS_30 =  'sogSxYW_29C_LacZ_sog12_data.mat'; % 02 07 12*
s2_30 =  'sog2enh_29C_LacZ_sog13_data.mat'; % 02 06 13


figure(1); set(gcf,'color','k');  title('sog 30C'); 
subplot(3,1,1); f = [0,0]; Is1 = endVrept(folder,sP_30,1.4,.3,f);
subplot(3,1,2); f = [0,2]; Is2 = endVrept(folder,sS_30,1.3,.3,f);
subplot(3,1,3); f = [1,2]; Is2x = endVrept(folder,s2_30,1.9,.3,f);


     imwrite(Is1,[fout,'sogP_30C_comp.tif'],'tif');
     imwrite(Is2,[fout,'sogS_30C_comp.tif'],'tif');
     imwrite(Is2x,[fout,'sog2enh_30C_comp.tif'],'tif');
     
%%
sP_22 =  'sogPxYW_22C_LacZ_sog18_data.mat'; %02  08
sS_22 =  'sogSxYW_22C_LacZ_sog17_data.mat' ; % 01
s2_22 =  'sog2enh_22C_LacZ_sog02_data.mat';%  02

figure(1); set(gcf,'color','k'); % title('sog 22C'); colordef black;
subplot(3,1,1); f = [1,2]; Is1_22 = endVrept(folder,sP_22,1.8,.3,f);
subplot(3,1,2); f = [1,0]; Is2_22 = endVrept(folder,sS_22,1.3,.3,f);
subplot(3,1,3); f = [0,0]; Is2x_22 = endVrept(folder,s2_22,3.5,.3,f);


     imwrite(Is1_22,[fout,'sogP_22C_comp.tif'],'tif');
     imwrite(Is2_22,[fout,'sogS_22C_comp.tif'],'tif');
     imwrite(Is2x_22,[fout,'sog2enh_22C_comp.tif'],'tif');
     

     
%% ------------ Analyze Missed Nuclei, hb-LacZ 30C ---------------------- %%


hb2_22 =  'hb2enh_22C_LacZ_hb01_data.mat';
hbP_22 =  'hbP_22C_LacZ_hb01_data.mat';
% 'C33_22C_LacZ_hb';

hbP_30 =  'C55_29C_LacZ_hb13_data.mat'; % 03  09  10
hbS_30 =  'C33_29C_LacZ_hb03_data.mat'; % 03
hb2_30 =  'hb2enh_29C_LacZ_hb02_data.mat'; % 01 02

figure(1); set(gcf,'color','k');  title('sog 30C'); 
subplot(3,1,1); f = [1,2]; Ihb1 = endVrept(folder,hbP_30,1.4,.3,f);
subplot(3,1,2); f = [0,0]; Ihb2 = endVrept(folder,hbS_30,1.3,.3,f);
subplot(3,1,3); f = [1,2]; Ihb2x = endVrept(folder,hb2_30,1.9,.3,f);

  imwrite(Ihb1,[fout,'hbP_LacZ_30C_comp.tif'],'tif');
     imwrite(Ihb2,[fout,'hbS_LacZ_30C_comp.tif'],'tif');
     imwrite(Ihb2x,[fout,'hb2enh_30C_comp.tif'],'tif');
     
     
     
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
     



   
   
   
   

