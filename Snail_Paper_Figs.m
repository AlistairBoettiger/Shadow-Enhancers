%%
clear all; 


fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Results/';

folder =  '/Volumes/Data/Lab Data/Shadow_data';
%imn = '12.2_snaD_CyO-HbZ_29C_Bgal-sim_sna_01.tif'; % ,01,14,13
imn = 'BAC7_snaD_HbZ_29C_sim-exon_sna_07.tif'; 
I = imread([folder,'/',imn]);

figure(1); clf; imshow(I);

[h,w] = size(I(:,:,1));

C = [1,1,0;
    0,0,0;
    0,.5,.5];

T = [.3,1;
    .12,1;
    .26,.65];

f = [1,2];

I = im_recolor(I,C,T,f) ;
figure(1); clf; imshow(I);

Io = imresize(I,.5);


figure(1); clf;  imshow(Io); 

imwrite(Io,[fout,'/',imn]); 

%%

folder =  '/Volumes/Data/Lab Data/Shadow_data';
%imn = '12.2_snaD_CyO-HbZ_29C_Bgal-sim_sna_01.tif'; % ,01,14,13
imn = 'MP12-21p_22C_LacZ_vn_sna_06.tif'; 
I = imread([folder,'/',imn]);

[h,w] = size(I(:,:,1));

C = [0,0,0 ;
    2.45,0,0;
    0,1,1]; 

I2 = I; 
I2(:,:,1) = imadjust(I(:,:,1),[.2,1],[0,1]);
I2(:,:,2) = imadjust(I(:,:,2),[.12,1],[0,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.14,.51],[0,1]);

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

 Ic = imflip(imflip(Ic,1),2); 

figure(1); clf;  imshow(Ic); 

imwrite(Ic,imn); 
%% Fig 2D  v2

 save_folder = '/Users/alistair/Documents/Berkeley/Levine Lab/Projects/Shadow Enhancers/Snail_Paper';
 folder =  '/Volumes/Data/Lab Data/Shadow_data';
% imn = '7.5_snaD_CyO-HbZ_22C_Bgal-sim_sna_01.tif';% 02, 10, 04
imn = 'MP07Hz_snaD_HbZ_22C_sim_sna-LacZ_01.tif';
I = imread([folder,'/',imn]);

[h,w] = size(I(:,:,1));

C = [1,1,0 ;
    1.8,0,0;
    0,1,1]; 

I2 = I; 
I2(:,:,1) = imadjust(I(:,:,1),[.25,.9],[0,1]);
I2(:,:,2) = imadjust(I(:,:,2),[.25,.6],[0,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.12,1],[0,1]);

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);
 Ic = imflip(imflip(Ic,2),1); 

figure(1); clf;  imshow(Ic); 
 imwrite(Ic,[save_folder,'/','rc_',imn]); 
 
 %% Fig 2C v2 

 save_folder = '/Users/alistair/Documents/Berkeley/Levine Lab/Projects/Shadow Enhancers/Snail_Paper';
 folder =  '/Volumes/Data/Lab Data/Shadow_data';
% imn = '7.5_snaD_CyO-HbZ_22C_Bgal-sim_sna_01.tif';% 02, 10, 04
imn = 'MP12Hz_snaD_HbZ_30C_sim_sna-LacZ_01.tif';
I = imread([folder,'/',imn]);

[h,w] = size(I(:,:,1));

C = [.81,.81,0 ;
    2.8,0,0;
    0,2,2]; 

I2 = I; 
I2(:,:,1) = imadjust(I(:,:,1),[.28,.8],[0,1]);
I2(:,:,2) = imadjust(I(:,:,2),[.051,.6],[0,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.1,1],[0,1]);

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);
 Ic = imflip(imflip(Ic,2),1); 

figure(1); clf;  imshow(Ic); 
 imwrite(Ic,[save_folder,'/','rc_',imn]); 
 
  %% Fig 2C v2 

 save_folder = '/Users/alistair/Documents/Berkeley/Levine Lab/Projects/Shadow Enhancers/Snail_Paper';
 folder =  '/Volumes/Data/Lab Data/Shadow_data';
% imn = '7.5_snaD_CyO-HbZ_22C_Bgal-sim_sna_01.tif';% 02, 10, 04
imn = 'MP12Hz_snaD_HbZ_22C_sim_sna-LacZ_02.tif';
I = imread([folder,'/',imn]);

[h,w] = size(I(:,:,1));

C = [1.4,1.4,0 ;
    1.4,0,0;
    0,1,1]; 

I2 = I; 
I2(:,:,1) = imadjust(I(:,:,1),[.2,.9],[0,1]);
I2(:,:,2) = imadjust(I(:,:,2),[.2,.6],[0,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.1,1],[0,1]);

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);
% Ic = imflip(imflip(Ic,2),1); 

figure(1); clf;  imshow(Ic); 
 imwrite(Ic,[save_folder,'/','rc_',imn]); 
  %% Fig 2A v2 

 save_folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Snail_Paper';
 folder =  '/Volumes/Data/Lab Data/Shadow_data';
% imn = '7.5_snaD_CyO-HbZ_22C_Bgal-sim_sna_01.tif';% 02, 10, 04
imn = 'MP12Hz_snaD_HbZ_22C_sim_sna-LacZ_08.tif';
I = imread([folder,'/',imn]);

[h,w] = size(I(:,:,1));

C = [1.4,1.4,0 ;
    1.7,0,0;
    0,1,1]; 

I2 = I; 
I2(:,:,1) = imadjust(I(:,:,1),[.22,.75],[0,1]);
I2(:,:,2) = imadjust(I(:,:,2),[.2,.55],[0,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.1,1],[0,1]);

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);
 Ic = imflip(imflip(Ic,2),1); 

figure(1); clf;  imshow(Ic); 
 imwrite(Ic,[save_folder,'/','rc_',imn]); 
   %% Fig 2E v2 

 save_folder = '/Users/alistair/Documents/Berkeley/Levine Lab/Projects/Shadow Enhancers/Snail_Paper';
 folder =  '/Volumes/Data/Lab Data/Shadow_data';
imn = 'BAC7_snaD_HbZ_29C_sim-exon_sna_07.tif';
I = imread([folder,'/',imn]);

[h,w] = size(I(:,:,1));

C = [1.,1.,0 ;
    1,0,0;
    0,1.25,1.25]; 

I2 = I; 
I2(:,:,1) = imadjust(I(:,:,1),[.3,1],[0,1]);
I2(:,:,2) = imadjust(I(:,:,2),[.3,.6],[0,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.23,1],[0,1]);

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);
 Ic = imflip(imflip(Ic,1),0); 

figure(1); clf;  imshow(Ic); 
 imwrite(Ic,[save_folder,'/','rc_',imn]); 
 
    %% Fig 2B v2 

 save_folder = '/Users/alistair/Documents/Berkeley/Levine Lab/Projects/Shadow Enhancers/Snail_Paper';
 folder =  '/Volumes/Data/Lab Data/Shadow_data';
imn = 'BAC7_snaD_HbZ_22C_sim-exon_sna_06.tif';
I = imread([folder,'/',imn]);

[h,w] = size(I(:,:,1));

C = [1.,1.,0 ;
    1,0,0;
    0,.71,.71]; 

I2 = I; 
I2(:,:,1) = imadjust(I(:,:,1),[.34,.62],[0,1]);
I2(:,:,2) = imadjust(I(:,:,2),[.18,.55],[0,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.23,.75],[0,1]);

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);
 Ic = imflip(imflip(Ic,1),0); 

figure(1); clf;  imshow(Ic); 
 imwrite(Ic,[save_folder,'/','rc_',imn]); 
 
 
 %% Fig 2C v2
 % exon sim sna[D] BAC7_snaD_HbZ_29C_sim-exon_sna_07.tif;
 
 % 22C sim, sna, gastrulate BAC7 % ;
 % 22C sim sna 
 
 
 
 %% Fig 4D

 save_folder = '/Users/alistair/Documents/Berkeley/Levine Lab/Projects/Shadow Enhancers/Snail_Paper';
 folder =  '/Volumes/Data/Lab Data/Shadow_data';

imn = 'MP07Hz_snaD_HbZ_30C_sim_sna-LacZ_07.tif';
I = imread([folder,'/',imn]);

[h,w] = size(I(:,:,1));

C = [1,1,0 ;
    1.6,0,0;
    0,1.2,1.2]; 

I2 = I; 
I2(:,:,1) = imadjust(I(:,:,1),[.35,1],[0,1]);
I2(:,:,2) = imadjust(I(:,:,2),[.25,.8],[0,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.05,1],[0,1]);

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);
 %Ic = imflip(imflip(Ic,2),1); 

figure(1); clf;  imshow(Ic); 
 imwrite(Ic,[save_folder,'/','rc_',imn]); 
 
%% snail Paper Figs
clear all;

folder =  '/Volumes/Data/Lab Data/Shadow_data';
save_folder = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Snail_Paper';
imn = 'snaSh1_sna_LacZ_04.tif';

I = imread([folder,'/',imn]);

[h,w] = size(I(:,:,1));

C = [0,.9,0 ;
    0,0,0;
    0,.8,.8]; 
T = [.2,1;
    .12 1;
    .1 1];

f = [0,2];

Ic  = im_recolor(I,C,T,f) ;
figure(1); clf; imshow(Ic); 
imwrite(Ic,[save_folder,'/','g-shadow_',imn],'tif');

C = [0,0,0 ;
    1.7,0,0;
    0,.8,.8]; 
Ic  = im_recolor(I,C,T,f) ;

figure(2);  imshow(Ic); 
figure(2); set(gcf,'color','k'); 
imwrite(Ic,[save_folder,'/','snail_',imn],'tif');

%%
folder =  '/Volumes/Data/Lab Data/Shadow_data';
imn = 'snaSh1_sna_LacZ_05.tif';

I = imread([folder,'/',imn]);

[h,w] = size(I(:,:,1));

C = [4,0,0 ;
    0,0,0;
    0,1.2,1.2]; 

I2 = I; 
I2(:,:,1) = imadjust(I(:,:,1),[.1,1],[0,1]);


% C = [1.1,0,0 ;
%     0,.9,.9;
%     .7,.7,.7];

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

Ic = imflip(imflip(Ic,2),1); 

figure(1); clf; % subplot(1,2,1);
imshow(Ic); 
 imwrite(Ic,'snaS_04.tif','tif'); 

C = [0,0,0 ;
    4,0,0;
    0,1.2,1.2]; 

I2 = I; 
I2(:,:,2) = imadjust(I(:,:,2),[.12,1],[0,1]);

% C = [1.1,0,0 ;
%     0,.9,.9;
%     .7,.7,.7];

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);


Ic = imflip(imflip(Ic,2),1); 

figure(2);%  subplot(1,2,2); 
imshow(Ic); 
figure(1); set(gcf,'color','k'); 
 imwrite(Ic,'sna_04.tif','tif'); 

%%  Figure 3
%% Fig 3A
clear all;
C = [2.5,2.5,0 ;
    0,0,0;
    0,.65,.65]; 

  rt = 'MP05_22C_y'; 
   num = '24'; % 47
   
folder =  '/Volumes/Data/Lab Data/Shadow_data';
save_folder = '/Users/alistair/Documents/Berkeley/Levine Lab/Snail_Paper/';
imn = [rt,'_sna_',num,'.tif']; % 

I = imread([folder,'/',imn]);

[h,w] = size(I(:,:,1));


I2 = I; 
I2(:,:,1) = imadjust(I(:,:,1),[.08,1],[0,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.14,1],[0,1]);

% C = [1.1,0,0 ;
%     0,.9,.9;
%     .7,.7,.7];

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

 Ic = imflip(Ic,1); 

figure(1); clf; % subplot(1,2,1);
imshow(Ic); 
 imwrite(Ic,[save_folder,'/',rt,'.tif'],'tif'); 
 
 load([folder,'/',rt,'_sna',num,'_data.mat']);  
 
      Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(255*L1);
     Io(:,:,2) = 1*uint8(255*L1) + .8*handles.In;
     Io(:,:,3) = .8*handles.In - Io(:,:,1);
     Ired = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
     
 Ired = imflip(Ired,1); 
 figure(2); clf; 
 imshow(Ired); hold on;
 
 imwrite(Ired,[save_folder,'/',rt,'_vis.tif'],'tif'); 
 
 Cell_bnd = imdilate(Cell_bnd,strel('disk',1)); 
 
   figure(3); clf;
    Io = uint8(zeros(h,w,3));
    Io(:,:,1) = imadd(uint8(255*L2n1a.*Cell_bnd),C(1,1)*I2(:,:,1)-uint8(255*L2n1a.*Cell_bnd) ); 
    Io(:,:,2) =  imadd(uint8(10*L2n1a.*Cell_bnd),C(1,2)*I2(:,:,1)-uint8(255*L2n1a.*Cell_bnd)  +C(3,2)*I2(:,:,3)); % imadd(uint8(255*Cell_bnd),1*handles.Im2);  %
    Io(:,:,3) =   imadd(uint8(10*L2n1a.*Cell_bnd),C(3,3)*I2(:,:,3)); % uint8(255*Cell_bnd); %
    % DI = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
    Izoom = Io(185:285,295:445,:);
  clf; imshow(Izoom);
  

 imwrite(Izoom,[save_folder,'/',rt,'_lab_zoom.tif'],'tif'); 
%%

% updated 06/30/10
 clear all;
C = [1.5,1.5,0 ;
    0,0,0;
    0,.65,.65]; 

  rt = 'MP06xYW_30C_y'; 
   num = '12'; % 47
i_folder =     '/Volumes/Data/Lab Data/Shadow_data';
d_folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed';
save_folder = '/Users/alistair/Documents/Berkeley/Levine Lab/Projects/Shadow Enhancers/Snail_Paper';
imn = [rt,'_sna_',num,'.tif']; % 

I = imread([i_folder,'/',imn]);

[h,w] = size(I(:,:,1));


I2 = I; 
I2(:,:,1) = imadjust(I(:,:,1),[.2,.8],[0,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.14,1],[0,1]);

% C = [1.1,0,0 ;
%     0,.9,.9;
%     .7,.7,.7];

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

 Ic = imflip(Ic,0); 

figure(1); clf; % subplot(1,2,1);
imshow(Ic); 
 imwrite(Ic,[save_folder,'/',rt,'.tif'],'tif'); 
 
 load([d_folder,'/',rt,'_sna',num,'_data.mat']);  
 
      Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(255*L1);
     Io(:,:,2) = 1*uint8(255*L1) + .8*handles.In;
     Io(:,:,3) = .8*handles.In - Io(:,:,1);
     Ired = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
     
 Ired = imflip(Ired,0); 
 figure(2); clf; 
 imshow(Ired); hold on;
 
 imwrite(Ired,[save_folder,'/',rt,'_vis.tif'],'tif'); 
 
 Cell_bnd = imdilate(Cell_bnd,strel('disk',1)); 
 
   figure(3); clf;
    Io = uint8(zeros(h,w,3));
    Io(:,:,1) = imadd(uint8(255*L2n1a.*Cell_bnd),C(1,1)*I2(:,:,1)-uint8(255*L2n1a.*Cell_bnd) ); 
    Io(:,:,2) =  imadd(uint8(10*L2n1a.*Cell_bnd),C(1,2)*I2(:,:,1)-uint8(255*L2n1a.*Cell_bnd)  +C(3,2)*I2(:,:,3)); % imadd(uint8(255*Cell_bnd),1*handles.Im2);  %
    Io(:,:,3) =   imadd(uint8(10*L2n1a.*Cell_bnd),C(3,3)*I2(:,:,3)); % uint8(255*Cell_bnd); %
    % DI = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
    Izoom = Io(2*185:2*285,2*250:2*400,:);
  clf; imshow(Izoom);
  

 imwrite(Izoom,[save_folder,'/',rt,'_lab_zoom.tif'],'tif'); 
 
 
 %%
clear all;

C = [2.5,2.5,0 ;
    0,0,0;
    0,.7,.7]; 

  rt = 'MP10_22C_y'; 
   num = '26'; % 07
   
i_folder =     '/Volumes/Data/Lab Data/Shadow_data';
d_folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed';
save_folder = '/Users/alistair/Documents/Berkeley/Levine Lab/Projects/Shadow Enhancers/Snail_Paper';

imn = [rt,'_sna_',num,'.tif']; % 

I = imread([i_folder,'/',imn]);

[h,w] = size(I(:,:,1));


I2 = I; 
I2(:,:,1) = imadjust(I(:,:,1),[.02,1],[0,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.06,1],[0,1]);

% C = [1.1,0,0 ;
%     0,.9,.9;
%     .7,.7,.7];

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

 Ic = imflip(Ic,2) ;

figure(1); clf; % subplot(1,2,1);
imshow(Ic); 
 
 load([d_folder,'/',rt,'_sna',num,'_data.mat']);  
 
      Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(255*L1);
     Io(:,:,2) = 1*uint8(255*L1) + .8*handles.In;
     Io(:,:,3) = .8*handles.In - Io(:,:,1);
     Ired = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
     
  Ired = imflip(Ired,2); 
 figure(2); clf; 
 imshow(Ired); hold on;
 
 
  figure(1); clf;
    Io = uint8(zeros(h,w,3));
    Io(:,:,1) = imadd(uint8(0*L2n1a.*Cell_bnd),C(1,1)*I2(:,:,1)); 
    Io(:,:,2) =  imadd(uint8(255*L2n1a.*Cell_bnd),C(1,2)*I2(:,:,1)+C(3,2)*I2(:,:,3)); % imadd(uint8(255*Cell_bnd),1*handles.Im2);  %
    Io(:,:,3) =   imadd(uint8(255*L2n1a.*Cell_bnd),C(3,3)*I2(:,:,3)); % uint8(255*Cell_bnd); %
    % DI = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
    Izoom = Io(85:185,375:525,:);
  clf; imshow(Izoom);

  
   imwrite(Ic,[save_folder,'/',rt,'.tif'],'tif'); 
 imwrite(Ired,[save_folder,'/',rt,'_vis.tif'],'tif'); 
 imwrite(Izoom,[save_folder,'/',rt,'_lab_zoom.tif'],'tif'); 
 
   %%
  clear all;

   rt = 'MP05xYW_30C_sna_y-full'; 
   num = '21'; % 18
   
i_folder =     '/Volumes/Data/Lab Data/Shadow_data';
d_folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed';
save_folder = '/Users/alistair/Documents/Berkeley/Levine Lab/Projects/Shadow Enhancers/Snail_Paper';

imn = [rt,'_',num,'.tif']; % 

 I = imread([i_folder,'/',imn]);

 load([d_folder,'/',rt,num,'_data.mat']);  

C = [0,0,0 ;
    1.5,1,0;
    0,1,1]; 

I = handles.It; 
[h,w] = size(I(:,:,1));

I2 = I; 
I2(:,:,2) = imadjust(I(:,:,2),[.18,1],[0,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.3,1],[0,1]);

% C = [1.1,0,0 ;
%     0,.9,.9;
%     .7,.7,.7];

figure(1); clf; imshow(I2);

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

Ic = imflip(imflip(Ic,1),2); 

figure(1); clf; % subplot(1,2,1);
imshow(Ic); 
 

 
      Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1.6*uint8(255*L1);
     Io(:,:,2) = 1.6*uint8(255*L1) + 1.25*handles.In;
     Io(:,:,3) = 1.25*handles.In - Io(:,:,1);
     Ired = uint8(bsxfun(@times,double(Io)/205,double(handles.In)));
   
     Ic = imflip(imflip(Ic,1),2); 
 Ired = imflip(imflip(Ired,1),2); 
 figure(2); clf; 
 imshow(Ired); hold on;  % (100:200,400:600,:)
 
 Cell_bnd = imdilate(Cell_bnd,strel('disk',2)); 

 figure(1); clf;
    Io = uint8(zeros(h,w,3));
    Io(:,:,1) = imadd(uint8(255*L2n1a.*Cell_bnd),C(2,1)*Ic(:,:,1)-uint8(255*L2n1a.*Cell_bnd) ); 
    Io(:,:,2) =  imadd(uint8(10*L2n1a.*Cell_bnd),C(2,2)*Ic(:,:,1)-uint8(255*L2n1a.*Cell_bnd)  +C(3,2)*Ic(:,:,3)); % imadd(uint8(255*Cell_bnd),1*handles.Im2);  %
    Io(:,:,3) =   imadd(uint8(10*L2n1a.*Cell_bnd),C(3,3)*Ic(:,:,3));
    Izoom =  Io(2*60:2*160,2*500:2*650,:);
    figure(1); % clf; imshow(Io);
  clf; imshow(Izoom);
  
   imwrite(Ic,[save_folder,'/',rt,'.tif'],'tif'); 
 imwrite(Ired,[save_folder,'/',rt,'_vis.tif'],'tif'); 
 imwrite(Izoom,[save_folder,'/',rt,'_lab_zoom.tif'],'tif'); 
  
 
  %%  old version 3 D
  clear all;

   rt = 'MP05_29C_y'; 
   num = '14'; % 26, 22 18, 14, 08
   
folder =  '/Volumes/Data/Lab Data/Shadow_data';
save_folder = '/Users/alistair/Documents/Berkeley/Levine Lab/Snail_Paper/';

imn = [rt,'_sna_',num,'.tif']; % 

% I = imread([folder,'/',imn]);

 load([folder,'/',rt,'_sna',num,'_data.mat']);  

C = [3,3,0 ;
    0,0,0;
    0,1,1]; 

I = handles.It; 
[h,w] = size(I(:,:,1));

I2 = I; 
I2(:,:,1) = imadjust(I(:,:,1),[.1,1],[0,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.1,1],[0,1]);

% C = [1.1,0,0 ;
%     0,.9,.9;
%     .7,.7,.7];

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

Ic = imflip(Ic,1); 

figure(1); clf; % subplot(1,2,1);
imshow(Ic); 
 

 
      Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1.6*uint8(255*L1);
     Io(:,:,2) = 1.6*uint8(255*L1) + 1.25*handles.In;
     Io(:,:,3) = 1.25*handles.In - Io(:,:,1);
     Ired = uint8(bsxfun(@times,double(Io)/145,double(handles.In)));
     
 Ired = imflip(Ired,1); 
 figure(2); clf; 
 imshow(Ired); hold on;  % (100:200,400:600,:)
 

 figure(1); clf;
    Io = uint8(zeros(h,w,3));
    Io(:,:,1) = imadd(uint8(255*L2n1a.*Cell_bnd),C(1,1)*I2(:,:,1)-uint8(255*L2n1a.*Cell_bnd) ); 
    Io(:,:,2) =  imadd(uint8(10*L2n1a.*Cell_bnd),C(1,2)*I2(:,:,1)-uint8(255*L2n1a.*Cell_bnd)  +C(3,2)*I2(:,:,3)); % imadd(uint8(255*Cell_bnd),1*handles.Im2);  %
    Io(:,:,3) =   imadd(uint8(10*L2n1a.*Cell_bnd),C(3,3)*I2(:,:,3));
    Izoom = Io(120:220,325:475,:);
  clf; imshow(Izoom);
  
   imwrite(Ic,[save_folder,'/',rt,'.tif'],'tif'); 
 imwrite(Ired,[save_folder,'/',rt,'_vis.tif'],'tif'); 
 imwrite(Izoom,[save_folder,'/',rt,'_lab_zoom.tif'],'tif'); 
  
 %%
   clear all;
   rt = 'MP10_29C_y'; 
   num = '18';% 08, 15, 18
   
folder =  '/Volumes/Data/Lab Data/Shadow_data';
imn = [rt,'_sna_',num,'.tif']; % 

I = imread([folder,'/',imn]);

[h,w] = size(I(:,:,1));

C = [2.5,2.5,0 ;
    0,0,0;
    0,.65,.65]; 

I2 = I; 
I2(:,:,1) = imadjust(I(:,:,1),[.02,1],[0,1]);


% C = [1.1,0,0 ;
%     0,.9,.9;
%     .7,.7,.7];

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

%Ic = imflip(Ic,1); 

figure(1); clf; % subplot(1,2,1);
imshow(Ic); 
 imwrite(Ic,[rt,'.tif'],'tif'); 
 
 load([folder,'/',rt,'_sna',num,'_data.mat']);  
 
      Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 2*uint8(255*L1);
     Io(:,:,2) = 2*uint8(255*L1) + 1*handles.In;
     Io(:,:,3) = 1*handles.In - Io(:,:,1);
     Ired = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
     
 % Ired = imflip(Ired,1); 
 figure(2); clf; 
 imshow(Ired); hold on;
 
 imwrite(Ired,[rt,'_vis.tif'],'tif'); 
 
  
   figure(3); clf;
    Io = uint8(zeros(h,w,3));
    Io(:,:,1) = imadd(uint8(0*L2n1a.*Cell_bnd),C(1,1)*I2(:,:,1)); 
    Io(:,:,2) =  imadd(uint8(255*L2n1a.*Cell_bnd),C(1,2)*I2(:,:,1)+C(3,2)*I2(:,:,3)); % imadd(uint8(255*Cell_bnd),1*handles.Im2);  %
    Io(:,:,3) =   imadd(uint8(255*L2n1a.*Cell_bnd),C(3,3)*I2(:,:,3)); % uint8(255*Cell_bnd); %
    % DI = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
    Izoom = Io(170:270,325:475,:);
  clf; imshow(Izoom);

 imwrite(Izoom,[rt,'_lab_zoom.tif'],'tif'); 
 
 
 
 
%% Figure 4 
 %% Snail rescue figs
folder = '/Volumes/Data/Lab Data/Shadow_data';
%root = '7.5_snaD_CyO-HbZ_29C_Bgal-sim_sna_14.tif'; red = [.14,.65]; 
root = '7.5_snaD_CyO-HbZ_29C_Bgal-sim_sna_01.tif'; red = [.17,.52]; 

I = imread([folder,'/',root]);
figure(1); clf; imshow(I(:,:,1));
[h,w] = size(I(:,:,1)); 

I2 = I;
I2(:,:,1) = imadjust(I(:,:,1),red);
I2(:,:,2) = imadjust(I(:,:,2),[0,1],[0,1]);
figure(1); clf; imshow(I2(:,:,1));
figure(1); clf; imshow(I2);

C = [1.4,1.4,0 ;
    1.8,0,0;
    0,1.2,1.2]; 

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

Ic = imflip(Ic,2);

figure(2); clf; imshow(Ic);
 imwrite(Ic,['snaBAC',root],'tif'); 

Izoom = Ic(110:260,700:950,:);
figure(3); clf; imshow(Izoom);

imwrite(Izoom,['snaBAC_zoom_',root],'tif'); 

%% alt fig 4F

 save_folder = '/Users/alistair/Desktop/Projects/Shadow Enhancers/Snail_Paper';
folder = '/Volumes/Data/Lab Data/Shadow_data';
root = '7.5_snaD_CyO-HbZ_29C_Bgal-sim_sna_14.tif'; red = [.14,.65]; 


I = imread([folder,'/',root]);
figure(1); clf; imshow(I(:,:,1));
[h,w] = size(I(:,:,1)); 

I2 = I;
I2(:,:,1) = imadjust(I(:,:,1),red);
I2(:,:,2) = imadjust(I(:,:,2),[0,1],[0,1]);
figure(1); clf; imshow(I2(:,:,1));
figure(1); clf; imshow(I2);

C = [1.4,1.4,0 ;
    1.8,0,0;
    0,1.2,1.2]; 

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

Ic = imflip(Ic,2);

figure(2); clf; imshow(Ic);
 imwrite(Ic,[save_folder,'/','snaBAC',root],'tif'); 

Izoom = Ic(110:260,700:950,:);
figure(3); clf; imshow(Izoom);

imwrite(Izoom,[save_folder,'/','snaBAC_zoom_',root],'tif'); 


%%  Rescue 29C
root2 = '12.2_snaD_CyO-HbZ_29C_Bgal-sim_sna_06.tif'; red = [.11,.45];
green = [0,1]; blue = [0,1]; 
C = [1,1,0 ;
    2.5,0,0;
    0,.54,.54]; 

I = imread([folder,'/',root2]);

[h,w] = size(I(:,:,1)); 

I2 = I;
I2(:,:,1) = imadjust(I(:,:,1),[.22,.8]);
I2(:,:,2) = imadjust(I(:,:,2),[.1,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.2,1]);


Ic2 = uint8(zeros(h,w,3));
Ic2(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic2(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic2(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

% Ic2 = imflip(Ic2,2); 

figure(4); clf; imshow(Ic2);
 imwrite(Ic2,[save_folder,'/','snaBAC',root2],'tif'); 
Izoom = Ic2(280:440,700:950,:);
figure(3); clf; imshow(Izoom);
imwrite(Izoom,[save_folder,'/','snaBAC_zoom_',root2],'tif'); 
%% rescue 29C
root2 = '12.2_snaD_CyO-HbZ_29C_Bgal-sim_sna_13.tif'; red = [.11,.75];
green = [0,1]; blue = [0,1]; 
C = [1.4,1.4,0 ;
    1.4,0,0;
    0,1,1]; 

I = imread([folder,'/',root2]);

[h,w] = size(I(:,:,1)); 

I2 = I;
I2(:,:,1) = imadjust(I(:,:,1),[.22,.8]);
I2(:,:,2) = imadjust(I(:,:,2),[.1,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.2,1]);


Ic2 = uint8(zeros(h,w,3));
Ic2(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic2(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic2(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

 Ic2 = imflip(Ic2,1); 

figure(4); clf; imshow(Ic2);
 imwrite(Ic2,[save_folder,'/','snaBAC_gastrulate',root2],'tif'); 
Izoom = Ic2(300:440,700:950,:);
figure(3); clf; imshow(Izoom);
imwrite(Izoom,[save_folder,'/','snaBAC_g_zoom_',root2],'tif'); 

%% rescue 29C  Fig 4E
%root2 = 'BAC12_snaD_HbZ_29C_sim-exon_sna_02.tif'; red = [.11,.75];

%root2 = 'BAC7_snaD_HbZ_29C_sim-exon_sna_08.tif'
root2 = 'BAC12_snaD_HbZ_29C_sim-exon_sna_01.tif'

green = [0,1]; blue = [0,1]; 
C = [0,1,0 ;
    1.8,0,0;
    0,1.3,1.3]; 

I = imread([folder,'/',root2]);

[h,w] = size(I(:,:,1)); 

I2 = I;
I2(:,:,1) = imadjust(I(:,:,1),[.22,.8]);
I2(:,:,2) = imadjust(I(:,:,2),[.1,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.2,1]);


Ic2 = uint8(zeros(h,w,3));
Ic2(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic2(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic2(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

  Ic2 = imflip(Ic2,1); 

figure(4); clf; imshow(Ic2);
 imwrite(Ic2,[save_folder,'/','snaBAC_gastrulate',root2],'tif'); 
Izoom = Ic2(300:440,700:950,:);
figure(3); clf; imshow(Izoom);
imwrite(Izoom,[save_folder,'/','snaBAC_g_zoom_',root2],'tif'); 

 %% Snail rescue late 
folder = '/Volumes/Data/Lab Data/Shadow_data';
root = '7.5_snaD_CyO-HbZ_29C_Bgal-sim_sna_10.tif'; red = [.14,.65]; 
%root = '7.5_snaD_CyO-HbZ_29C_Bgal-sim_sna_01.tif'; red = [.17,.52]; 

I = imread([folder,'/',root]);
figure(1); clf; imshow(I(:,:,1));
[h,w] = size(I(:,:,1)); 

I2 = I;
I2(:,:,1) = imadjust(I(:,:,1),[.25,.5]);
I2(:,:,2) = imadjust(I(:,:,2),[.2,1],[0,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.03,1],[0,1]);


C = [1.25,1.25,0 ;
    1.8,0,0;
    0,1.9,1.9]; 

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

% Ic = imflip(Ic,2);

figure(2); clf; imshow(Ic);
imwrite(Ic,[save_folder,'/','snaBAC',root],'tif'); 




%% dorsal hets
clear all;
rt = 'MP10xdl6_25C_y'; 

   num = '02';% 08, 15, 18
   
folder =  '/Volumes/Data/Lab Data/Shadow_data';
save_folder = '/Users/alistair/Documents/Berkeley/Levine Lab/Snail_Paper';
imn = [rt,'_hb_',num,'.tif']; % 

I = imread([folder,'/',imn]);

[h,w] = size(I(:,:,1));

C = [1.5,1.5,0 ;
    0,0,0;
    0,.8,.8]; 

I2 = I; 
I2(:,:,1) = imadjust(I(:,:,1),[.08,1],[0,1]);


% C = [1.1,0,0 ;
%     0,.9,.9;
%     .7,.7,.7];

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

Ic = imflip(Ic,1); 

figure(1); clf; % subplot(1,2,1);
imshow(Ic); 
 imwrite(Ic,[save_folder,'/',rt,'.tif'],'tif'); 
 
 load([folder,'/',rt,'_hb',num,'_data.mat']);  
 
      Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1.5*C(1,1)*uint8(255*L1);
     Io(:,:,2) = 1.5*C(1,2)*uint8(255*L1) + 1.5*C(3,2)*handles.In;
     Io(:,:,3) = 1.5*C(3,3)*handles.In - Io(:,:,1);
     Ired = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
     
  Ired = imflip(Ired,1); 
 figure(2); clf; 
 imshow(Ired); hold on;
 
 imwrite(Ired,[save_folder,'/',rt,'_vis.tif'],'tif'); 
 
    L2n1a = (Reg1 - L1) > 0; 
 
   figure(3); clf;
   Cell_bnd = imdilate(Cell_bnd,strel('disk',1));
   
    Io = uint8(zeros(h,w,3));
    Io(:,:,1) = imadd(uint8(255*L2n1a.*Cell_bnd),C(1,1)*I2(:,:,1)-uint8(255*L2n1a.*Cell_bnd) ); 
    Io(:,:,2) =  imadd(uint8(10*L2n1a.*Cell_bnd),C(1,2)*I2(:,:,1)-uint8(255*L2n1a.*Cell_bnd)  +C(3,2)*I2(:,:,3)); % imadd(uint8(255*Cell_bnd),1*handles.Im2);  %
    Io(:,:,3) =   imadd(uint8(10*L2n1a.*Cell_bnd),C(3,3)*I2(:,:,3)); % uint8(255*Cell_bnd); %
    % DI = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
    Izoom = Io(160:300,425:625,:); % 170:270,325:475,
  clf; imshow(Izoom);

 imwrite(Izoom,[save_folder,'/',rt,'_lab_zoom.tif'],'tif'); 
 
 
 %% dorsal hets
clear all;
rt = 'MP05xdl6_25C_y'; 

   num = '01';% 
   
folder =  '/Volumes/Data/Lab Data/Shadow_data';
save_folder = '/Users/alistair/Documents/Berkeley/Levine Lab/Snail_Paper';
imn = [rt,'_hb_',num,'.tif']; % 

I = imread([folder,'/',imn]);

[h,w] = size(I(:,:,1));

C = [1.5,1.5,0 ;
    0,0,0;
    0,.8,.8]; 

I2 = I; 
I2(:,:,1) = imadjust(I(:,:,1),[.1,1],[0,1]);


% C = [1.1,0,0 ;
%     0,.9,.9;
%     .7,.7,.7];

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

Ic = imflip(imflip(Ic,1),2); 

figure(1); clf; % subplot(1,2,1);
imshow(Ic); 
 imwrite(Ic,[save_folder,'/',rt,'.tif'],'tif'); 
 
 load([folder,'/',rt,'_hb',num,'_data.mat']);  
 
      Io = uint8(zeros(h,w,3));
     Io(:,:,1) = C(1,1)*uint8(255*L1);
     Io(:,:,2) = C(1,2)*uint8(255*L1) + 1.5*C(3,2)*handles.In;
     Io(:,:,3) = 1.5*C(3,3)*handles.In - Io(:,:,1);
     Ired = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
     
  Ired = imflip(imflip(Ired,1),2); 
 figure(2); clf; 
 imshow(Ired); hold on;
 
 imwrite(Ired,[save_folder,'/',rt,'_vis.tif'],'tif'); 
 
    L2n1a = (Reg1 - L1) > 0; 
 
   figure(3); clf;
   Cell_bnd = imdilate(Cell_bnd,strel('disk',1));
   
    Io = uint8(zeros(h,w,3));
    Io(:,:,1) = imadd(uint8(255*L2n1a.*Cell_bnd),C(1,1)*I2(:,:,1)-uint8(255*L2n1a.*Cell_bnd) ); 
    Io(:,:,2) =  imadd(uint8(10*L2n1a.*Cell_bnd),C(1,2)*I2(:,:,1)-uint8(255*L2n1a.*Cell_bnd)  +C(3,2)*I2(:,:,3)); % imadd(uint8(255*Cell_bnd),1*handles.Im2);  %
    Io(:,:,3) =   imadd(uint8(10*L2n1a.*Cell_bnd),C(3,3)*I2(:,:,3)); % uint8(255*Cell_bnd); %
    % DI = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
    Izoom = Io(60:200,600:800,:); % 170:270,325:475,
  clf; imshow(Izoom);

 imwrite(Izoom,[save_folder,'/',rt,'_lab_zoom.tif'],'tif'); 
 
 
 
  %% Ip Snail rescue late 
folder = '/Volumes/Data/Lab Data/Shadow_data';
root = 'snaE_snaN_sna_sim_03.tif';
save_folder = '/Users/alistair/Documents/Berkeley/Levine Lab/Snail_Paper';
%root = '7.5_snaD_CyO-HbZ_29C_Bgal-sim_sna_01.tif'; red = [.17,.52]; 

I = imread([folder,'/',root]);
figure(1); clf; imshow(I(:,:,1));
[h,w] = size(I(:,:,1)); 

I2 = I;
I2(:,:,1) = imadjust(I(:,:,1),[.2,1]);
I2(:,:,2) = imadjust(I(:,:,2),[.1,.38],[0,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.03,1],[0,1]);


C = [1.2,0,0 ;
    1.3,1.3,0;
    0,1.2,1.2]; 

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

 Ic = imflip(imflip(Ic,0),0);


figure(2); clf; imshow(Ic);
imwrite(Ic,[save_folder,'/','snaBAC',root],'tif'); 


    Izoom = Ic(120:280,700:900,:); % 170:270,325:475,
 figure(3); clf; imshow(Izoom);

 imwrite(Izoom,[save_folder,'/',root,'_lab_zoom.tif'],'tif'); 

%% Ip Snail rescue late 
folder = '/Volumes/Data/Lab Data/Shadow_data';
root = 'snaE_snaN_sna_sim_02.tif';
save_folder = '/Users/alistair/Documents/Berkeley/Levine Lab/Snail_Paper';
%root = '7.5_snaD_CyO-HbZ_29C_Bgal-sim_sna_01.tif'; red = [.17,.52]; 

I = imread([folder,'/',root]);
figure(1); clf; imshow(I(:,:,1));
[h,w] = size(I(:,:,1)); 

I2 = I;
I2(:,:,1) = imadjust(I(:,:,1),[.2,1]);
I2(:,:,2) = imadjust(I(:,:,2),[.3,1],[0,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.3,1],[0,1]);


C = [1,0,0 ;
    1.3,1.3,0;
    0,.7,.7]; 

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

 Ic = imflip(imflip(Ic,1),2);

figure(2); clf; imshow(Ic);
imwrite(Ic,[save_folder,'/','snaBAC',root],'tif'); 


    Izoom = Ic(240:400,700:900,:); % 170:270,325:475,
 figure(3); clf; imshow(Izoom);

 imwrite(Izoom,[save_folder,'/',root,'_lab_zoom.tif'],'tif'); 

 %% Supplement
 save_folder = '/Users/alistair/Desktop/Projects/Shadow Enhancers/Snail_Paper';
 
 folder = '/Volumes/Data/Lab Data/Shadow_data';
 yintR = 'MP10_22C_y_sna_25.tif'; % 'MP05_29C_y_sna_26.tif'; %'
 yexonR =  'MP05_29C_sna_y-full_04.tif'; %'MP05_29C_sna_y-full_02.tif';
 
 yint = imread([folder,'/',yintR]);
 yexon = imread([folder,'/',yexonR]);
 
 
 
 
 yint2 = yint;
 yint2(:,:,1) = imadjust(yint(:,:,1), [.1,.34]);
 yint2(:,:,3) = imadjust(yint(:,:,3), [.25,1]); 
 
 yint2(:,:,1) = yint2(:,:,1);
 yint2(:,:,2) = yint2(:,:,1) + .8*yint2(:,:,3);
 yint2(:,:,3) = .8*yint2(:,:,3);
 yint2 = imflip(yint2,2); 
 
 yizoom = yint2(200:300,100:300,:);
 
 figure(1); clf; imshow(yint2);
 figure(3); clf; imshow(yizoom);
 
 
 yexon2 = yexon;
 yexon2(:,:,2) = imadjust(yexon(:,:,2),[.031,.5]);
 yexon2(:,:,3) = imadjust(yexon(:,:,3),[.3,1]); 
 
 yexon2(:,:,1) = yexon2(:,:,2);
 yexon2(:,:,2) = yexon2(:,:,2) + yexon2(:,:,3);
 yexon2(:,:,3) = yexon2(:,:,3);
 
 yexon2 = imflip(yexon2,2); 
 
 figure(2); clf; imshow(yexon2);
 
 yezoom = yexon2(300:400,100:300,:);
 figure(4); clf; imshow(yezoom);
 
 imwrite(yexon2,[save_folder,'/','ye_',yexonR],'tif');
 imwrite(yezoom,[save_folder,'/','yezoom_',yexonR],'tif');
 
  imwrite(yint2,[save_folder,'/','yi_',yintR],'tif');
  imwrite(yizoom,[save_folder,'/','yizoom_',yintR],'tif');
 
  bk = yizoom; bk(:,:,1) = zeros(101,201);
  bk(:,:,2) = zeros(101,201); 
  bk(:,:,3) = zeros(101,201); 
  
 imwrite(bk,[save_folder,'/','blk.tif']); 
 
 
%%
 save_folder = '/Users/alistair/Desktop/Projects/Shadow Enhancers/Snail_Paper';
 
 folder = '/Volumes/Data/Lab Data/Shadow_data';
B7_name = 'BAC7_snaD_HbZ_29C_sim-exon_sna_05.tif'; % 'MP05_29C_y_sna_26.tif'; %'
B12_name = 'BAC12_snaD_HbZ_29C_sim-exon_sna_04.tif'; %'MP05_29C_sna_y-full_02.tif';
 
 B7 = imread([folder,'/',B7_name]);
 B12 = imread([folder,'/',B12_name]);
 
 
 
I = B7; 

I2 = I; 
[h,w] = size(I(:,:,1));
 
I2(:,:,1) = imadjust(I(:,:,1),[.45,1]);
I2(:,:,2) = imadjust(I(:,:,2),[.4,1],[0,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.2,1],[0,1]);


 figure(1); clf; imshow(I2);

C = [0,1.3,0 ;
    1.6,0,0;
    0,.8,.8]; 

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

% Ic = imflip(Ic,2);
 
B7_out = Ic;
 figure(1); clf; imshow(B7_out);
  B7zoom = B7_out(75:200,850:1050,:);
 figure(3); clf; imshow(B7zoom);



I = imflip(B12,1) ;
I2 = I; 
[h,w] = size(I(:,:,1));
 
I2(:,:,1) = imadjust(I(:,:,1),[.3,1]);
I2(:,:,2) = imadjust(I(:,:,2),[.2,1],[0,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.2,1],[0,1]);


 figure(2); clf; imshow(I2);

C = [0,1.4,0 ;
    1.6,0,0;
    0,1,1]; 

Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

 B12_out = Ic; 
 figure(2); clf; imshow(B12_out);
   B12zoom = B12_out(75:200,750:950,:);
 figure(4); clf; imshow(B12zoom);

 
 
   imwrite(B7_out,[save_folder,'/','B7_',B7_name],'tif');
  imwrite(B12_out,[save_folder,'/','B12_',B12_name],'tif');
  
     imwrite(B7zoom,[save_folder,'/','B7zoom_',B7_name],'tif');
  imwrite(B12zoom,[save_folder,'/','B12zoom_',B12_name],'tif');
  
  %% Supplement
  
  %% rescue 29C
root2 = '7.5_snaD_CyO-HbZ_22C_Bgal-sim_sna_04.tif'; red = [.11,.75];
green = [0,1]; blue = [0,1]; 
C = [2,2,0 ;
    2.5,0,0;
    0,1.3,1.3]; 

I = imread([folder,'/',root2]);

[h,w] = size(I(:,:,1)); 

I2 = I;
I2(:,:,1) = imadjust(I(:,:,1),[.22,1]);
I2(:,:,2) = imadjust(I(:,:,2),[.1,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.2,1]);


Ic2 = uint8(zeros(h,w,3));
Ic2(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic2(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic2(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

 Ic2 = imflip(Ic2,1); 

figure(4); clf; imshow(Ic2);
 imwrite(Ic2,[save_folder,'/','snaBAC_gastrulate',root2],'tif'); 
Izoom = Ic2(160:300,100:350,:);
figure(3); clf; imshow(Izoom);
imwrite(Izoom,[save_folder,'/','snaBAC_g_zoom_',root2],'tif'); 
  

  %% rescue 29C
root2 = 'BAC12_snaD_HbZ_29C_sim-exon_sna_06.tif'; 
C = [0,1,0 ;
    1.3,0,0;
    0,1,1]; 

I = imread([folder,'/',root2]);

[h,w] = size(I(:,:,1)); 

I2 = I;
I2(:,:,1) = imadjust(I(:,:,1),[.53,.89]);
I2(:,:,2) = imadjust(I(:,:,2),[.1,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.2,1]);


Ic2 = uint8(zeros(h,w,3));
Ic2(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic2(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic2(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

 Ic2 = imflip(Ic2,2); 

figure(4); clf; imshow(Ic2);
 imwrite(Ic2,[save_folder,'/','snaBAC_gastrulate',root2],'tif'); 
Izoom = Ic2(60:200,100:350,:);
figure(3); clf; imshow(Izoom);
imwrite(Izoom,[save_folder,'/','snaBAC_g_zoom_',root2],'tif'); 
 

  %% rescue 29C
root2 = 'BAC7_snaD_HbZ_29C_sim-exon_sna_02.tif'; 
C = [0,1,0 ;
    5,0,0;
    0,1,1]; 

I = imread([folder,'/',root2]);

[h,w] = size(I(:,:,1)); 

I2 = I;
I2(:,:,1) = imadjust(I(:,:,1),[.53,.89]);
I2(:,:,2) = imadjust(I(:,:,2),[.1,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.2,1]);


Ic2 = uint8(zeros(h,w,3));
Ic2(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic2(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic2(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

 % Ic2 = imflip(Ic2,2); 

figure(4); clf; imshow(Ic2);
 imwrite(Ic2,[save_folder,'/','snaBAC_gastrulate',root2],'tif'); 




  %% rescue 29C
root2 = 'BAC12_snaD_HbZ_29C_c_vnd_sna-LacZ_06.tif'; 
C = [0,0,0 ;
   3.4,0,0;
    0,1.5,1.5]; 

I = imread([folder,'/',root2]);

[h,w] = size(I(:,:,1)); 

I2 = I;
I2(:,:,1) = imadjust(I(:,:,1),[.53,.89]);
I2(:,:,2) = imadjust(I(:,:,2),[.1,1]);
I2(:,:,3) = imadjust(I(:,:,3),[.35,1]);


Ic2 = uint8(zeros(h,w,3));
Ic2(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic2(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic2(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);

  Ic2 = imflip(Ic2,1); 

figure(4); clf; imshow(Ic2);
 imwrite(Ic2,[save_folder,'/','snaBAC_gastrulate',root2],'tif'); 

  