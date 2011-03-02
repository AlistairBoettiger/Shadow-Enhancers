
%% Figs_hb-embs

% Figure plotting of hb embryos for presentation
 
clear all



fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Results/';

folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed/';
kr1 = 'krCD1_22C_LacZ_kr08_data.mat'; % 06 08
kr2 = 'krCD2_22C_LacZ_kr21_data.mat'; %   12* , 19, 21, 22
kr2x = 'kr2enh_22C_LacZ_kr09_data.mat';% 09
kni5 = 'kni_5p_22C_LacZ_kni18_data.mat'; % 18  05 12 21  23
kniI = 'kni_int_22C_LacZ_kni15_data.mat'; % 06 22*, 23, 11, 13, 15, 18, 33
kni2x = 'kni_2enh_22C_LacZ_kni19_data.mat';% 16 19 04 06 10 22, 24


    

figure(1); set(gcf,'color','k');
subplot(3,1,1); f = [1,2]; Ikr1 = endVrept(folder,kr1,1.9,.3,f);
subplot(3,1,2); f = [0,0]; Ikr2 = endVrept(folder,kr2,1.9,.3,f);
subplot(3,1,3); f = [1,2]; Ikr2x = endVrept(folder,kr2x,1.3,.3,f);
% 
figure(2); set(gcf,'color','k');
subplot(3,1,1); f = [0,2]; Ikni5 = endVrept(folder,kni5,2.5,.3,f);
subplot(3,1,2); f = [0,2]; IkniI = endVrept(folder,kniI,2,.3,f);
subplot(3,1,3); f = [0,2]; Ikni2x = endVrept(folder,kni2x,2.5,.3,f);
     
     

     
 %%     
     
 % write images to folder
        imwrite(Ikr1,[fout,'krCD1_comp.tif'],'tif');
     imwrite(Ikr2,[fout,'krCD2_comp.tif'],'tif');
     imwrite(Ikr2x,[fout,'kr2enh_comp.tif'],'tif');
     imwrite(Ikni5,[fout,'kni_5p_comp.tif'],'tif');
        imwrite(IkniI,[fout,'kni_int_comp.tif'],'tif');
     imwrite(Ikni2x,[fout,'kni_2enh_comp.tif'],'tif');
     
  %%
 
  kni5 = 'kni_5p_22C_LacZ_kni18_data.mat'; % 18* 12   %   12 21  23   % cc13   14 16  19 20   cc11=17
kniI = 'kni_int_22C_LacZ_kni07_data.mat'; % 06~  22  24  02e  04e  07e ~14 16 %   22*, 23, 11, 13, 15, 18, 33  % cc13   01  10 11  12
kni2x = 'kni_2enh_22C_LacZ_kni16_data.mat';% 16 06 19* ~22  %    % cc13 08   17
  
  
  figure(2); set(gcf,'color','k');
subplot(3,1,1); f = [0,2]; Ikni5 = ShowEctopNuc(folder,kni5,2.9,.3,f);
subplot(3,1,2); f = [0,0]; IkniI = ShowEctopNuc(folder,kniI,1.8,.3,f);
subplot(3,1,3); f = [0,0]; Ikni2x = ShowEctopNuc(folder,kni2x,3.5,.3,f);

     imwrite(Ikni5,[fout,'kni_5p_ectop.tif'],'tif');
        imwrite(IkniI,[fout,'kni_int_ectop.tif'],'tif');
     imwrite(Ikni2x,[fout,'kni_2enh_ectop.tif'],'tif');

     
     
%%


load([folder,kr1]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) = .8*handles.In- Io(:,:,1);
     Io(:,:,3) = .8*handles.In - Io(:,:,1);
     Ikr1 = uint8(bsxfun(@times,double(Io)/205,double(handles.In)));
     Ikr1 = imflip(Ikr1,2);
     age = getage(H,cent)
           
     load([folder,kr2]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(160*L1);
     Io(:,:,2) =  .8*handles.In- Io(:,:,1);
     Io(:,:,3) = .8*handles.In - Io(:,:,1);
     Ikr2 = uint8(bsxfun(@times,double(Io)/95,double(handles.In)));
     Ikr2 = imflip(imflip(Ikr2,2),1);
     age = getage(H,cent)
     
     load([folder,kr2x]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) =  .8*handles.In- Io(:,:,1);
     Io(:,:,3) = .8*handles.In - Io(:,:,1);
     Ikr2x= uint8(bsxfun(@times,double(Io)/185,double(handles.In)));
       Ikr2x = imflip(imflip(Ikr2x,2),1);
     age = getage(H,cent)
     
     
     load([folder,kni5]); 
     Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) = .8*handles.In- Io(:,:,1);
     Io(:,:,3) = .8*handles.In - Io(:,:,1);
     Ikni5 = uint8(bsxfun(@times,double(Io)/75,double(handles.In)));
     Ikni5 = imflip(Ikni5,2); 
     age = getage(H,cent)
           
     load([folder,kniI]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) =  .8*handles.In- Io(:,:,1);
     Io(:,:,3) = .8*handles.In - Io(:,:,1);
     IkniI = uint8(bsxfun(@times,double(Io)/165,double(handles.In)));
     IkniI = imflip(IkniI,2); 
     age = getage(H,cent)
     
     load([folder,kni2x]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) =  .8*handles.In- Io(:,:,1);
     Io(:,:,3) = .8*handles.In - Io(:,:,1);
     Ikni2x= uint8(bsxfun(@times,double(Io)/55,double(handles.In)));
     age = getage(H,cent)
     
     figure(1); clf;
     subplot(2,3,1); imshow(Ikr1);
     subplot(2,3,2); imshow(Ikr2);
     subplot(2,3,3); imshow(Ikr2x);
     subplot(2,3,4); imshow(Ikni5);
     subplot(2,3,5); imshow(IkniI);
     subplot(2,3,6); imshow(Ikni2x);
     set(gcf,'color','k');
     
     imwrite(Ikr1,[fout,'krCD1.tif'],'tif');
     imwrite(Ikr2,[fout,'krCD2.tif'],'tif');
     imwrite(Ikr2x,[fout,'kr2enh.tif'],'tif');
     imwrite(Ikni5,[fout,'kni_5p.tif'],'tif');
        imwrite(IkniI,[fout,'kni_int.tif'],'tif');
     imwrite(Ikni2x,[fout,'kni_2enh.tif'],'tif');
     