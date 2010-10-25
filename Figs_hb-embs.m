
%% Figs_hb-embs

% Figure plotting of hb embryos for presentation
 
clear all



fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Results/';


%%

   rt = 'MP02_22C_y'; 
   num = '62';%
   
folder =  '/Volumes/Data/Lab Data/Shadow_data';
imn = [rt,'_hb_',num,'.tif'];  

I = imread([folder,'/',imn]);

[h,w] = size(I(:,:,1));

C = [2.5,2.5,0 ;
    0,0,0;
    0,.65,.65]; 

I2 = I; 
I2(:,:,1) = imadjust(I(:,:,1),[.02,1],[0,1]);


Ic = uint8(zeros(h,w,3));
Ic(:,:,1) = C(1,1)*I2(:,:,1) + C(2,1)*I2(:,:,2) + C(3,1)*I2(:,:,3);
Ic(:,:,2) = C(1,2)*I2(:,:,1) + C(2,2)*I2(:,:,2) + C(3,2)*I2(:,:,3);
Ic(:,:,3) = C(1,3)*I2(:,:,1) + C(2,3)*I2(:,:,2) + C(3,3)*I2(:,:,3);


figure(1); clf; % subplot(1,2,1);
imshow(Ic); 
 imwrite(Ic,[rt,'.tif'],'tif'); 
 
 load([folder,'/',rt,'_hb',num,'_data.mat']);  
 
      Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 2*uint8(255*L1);
     Io(:,:,2) = 2*uint8(255*L1) + 1*handles.In;
     Io(:,:,3) = 1*handles.In - Io(:,:,1);
     Ired = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
     

 figure(2); clf; 
 imshow(Ired); hold on;
 
 imwrite(Ired,[fout,rt,'_vis.tif'],'tif'); 
 
  
   figure(3); clf;
    Io = uint8(zeros(h,w,3));
    Io(:,:,1) = imadd(uint8(0*L2n1a.*Cell_bnd),C(1,1)*I2(:,:,1)); 
    Io(:,:,2) =  imadd(uint8(255*L2n1a.*Cell_bnd),C(1,2)*I2(:,:,1)+C(3,2)*I2(:,:,3)); % imadd(uint8(255*Cell_bnd),1*handles.Im2);  %
    Io(:,:,3) =   imadd(uint8(255*L2n1a.*Cell_bnd),C(3,3)*I2(:,:,3)); % uint8(255*Cell_bnd); %
    % DI = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
    Izoom = Io(170:270,325:475,:);
  clf; imshow(Izoom);

 imwrite(Izoom,[fout,rt,'_lab_zoom.tif'],'tif'); 
 































%% Young stochastic MP09 at 30C


folder = '/Volumes/Data/Lab Data/Shadow_data';
emb = 'BAC09_30C_y_hb_73.tif';
name = 'BAC09_30C_y_hb_73.tif';

I = imread([folder,'/',emb]);  

C = [1,1,0;
    3,0,0;
    0,.5,.5];

T = [.1,1;
    .12,1;
    .1,.45];

f = [1,2];

I = im_recolor(I,C,T,f) ;
figure(1); clf; imshow(I);

Io = imresize(I,1);


imwrite(Io,[fout,name]);


%%  Young stochastic MP09 at 30C b

folder = '/Volumes/Data/Lab Data/Shadow_data';
emb = 'BAC09_30C_y_hb_69.tif';
name = 'BAC09_30C_y_hb_69.tif';

I = imread([folder,'/',emb]);  

C = [1,1,0;
    3,0,0;
    0,.5,.5];

T = [.1,.7;
    .22,.7;
    .1,.45];

f = [1,2];

I = im_recolor(I,C,T,f) ;
figure(1); clf; imshow(I);

Io = imresize(I,1);


imwrite(Io,[fout,name]);


%% MP02 reaches tip

folder = '/Volumes/Data/Lab Data/Shadow_data';
emb = 'MP02_30C_LacZ_hb_40.tif';
name = 'MP02_Aexpr.tif';

I = imread([folder,'/',emb]);  

C = [1,1,0;
    0,0,0;
    0,.5,.5];

T = [.13,1;
    .1,1;
    .1,1];

f = [1,1];

I = im_recolor(I,C,T,f) ;
figure(1); clf; imshow(I);

Io = imresize(I,.5);


imwrite(Io,[fout,name]);

%% MP02 reaches tip B

folder = '/Volumes/Data/Lab Data/Shadow_data';
emb = 'MP02_22C_y_hb_39.tif';
name = 'MP02_Aexpr_C.tif';

I = imread([folder,'/',emb]);  

C = [0,0,0;
    1.3,1.3,0;
    0,.65,.65];

T = [.16,1;
    .1,1;
    .1,1];

f = [1,1];

I = im_recolor(I,C,T,f) ;
figure(1); clf; imshow(I);

Io = imresize(I,.5);


imwrite(Io,[fout,name]);



%% MP01 does not reach tip

folder = '/Volumes/Data/Lab Data/Shadow_data';
emb = 'MP01_22C_y_hb_14.tif';
name = 'MP01_noAexpr.tif';

I = imread([folder,'/',emb]);  

C = [1,1,0;
    0,0,0;
    0,.5,.5];

T = [.13,1;
    .1,1;
    .1,1];

f = [2,0];

I = im_recolor(I,C,T,f) ;
figure(1); clf; imshow(I);

Io = imresize(I,.5);


imwrite(Io,[fout,name]);

%% MP01 does not reach tip

folder = '/Volumes/Data/Lab Data/Shadow_data';
emb = 'MP01_22C_y_hb_24.tif';
name = 'MP01_noAexpr_B.tif';

I = imread([folder,'/',emb]);  

C = [1,1,0;
    0,0,0;
    0,.45,.45];

T = [.1,.57;
    .1,1;
    .1,1];

f = [2,0];

I = im_recolor(I,C,T,f) ;
figure(1); clf; imshow(I);

Io = imresize(I,.5);


imwrite(Io,[fout,name]);


