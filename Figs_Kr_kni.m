
%% Figs_hb-embs

% Figure plotting of hb embryos for presentation
 
clear all



fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Results/';

folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed/';
kr1 = 'krCD1_22C_LacZ_kr06_data.mat'; % 
kr2 = 'krCD2_22C_LacZ_kr08_data.mat'; %  
kr2x = 'kr2enh_22C_LacZ_kr09_data.mat';% 
kni5 = 'kni_5p_22C_LacZ_kni18_data.mat'; %
kniI = 'kni_int_22C_LacZ_kni22_data.mat'; % 
kni2x = 'kni_2enh_22C_LacZ_kni16_data.mat';%

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
     