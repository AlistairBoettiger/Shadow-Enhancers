


clear all; 


%%
clear all;
folder =  '/Volumes/Data/Lab Data/Shadow_data';
prj = '/Users/alistair/Documents/Berkeley/Levine Lab/Projects/Shadow Enhancers/';
fdir = 'Results-July/';
imn = 'Kr2enh_LacZ_kr_11.tif'; 
I = imread([folder,'/',imn]);

C = [3,0,0 ;
    0,2,0;
    0,.75,.75]; 

T = [.18,.7;
     .12, .65;
     .22, 1];  % min max thresholds for channels
 
 f = [0,1]; % imflip horizontal = 1, vertical = 2, no flip = 0

 Io = im_recolor(I,C,T,f); 
Io = imresize(Io,.3); 

figure(1); clf; imshow(Io); 

In = I; 

In(:,:,2) = uint8(zeros(size(I(:,:,1)) ));

figure(1); clf; imshow(In);


% imwrite(Io,[prj,fdir,'rc_',imn],'tif');

%%
clear all;
folder =  '/Volumes/Data/Lab Data/Shadow_data';
prj = '/Users/alistair/Documents/Berkeley/Levine Lab/Projects/Shadow Enhancers/';
fdir = 'Results-July/';
imn = 'KrCD1_LacZ_kr_03.tif'; 
I = imread([folder,'/',imn]);

C = [3,0,0 ;
    0,2,0;
    0,.75,.75]; 

T = [.05,.5;
     .1, .75;
     .24, 1];  % min max thresholds for channels
 
 f = [0,0]; % imflip horizontal = 1, vertical = 2, no flip = 0

 Io = im_recolor(I,C,T,f); 
Io = imresize(Io,.3); 

figure(1); clf; imshow(Io); 

% imwrite(Io,[prj,fdir,'rc_',imn],'tif');

%%
clear all;
folder =  '/Volumes/Data/Lab Data/Shadow_data';
prj = '/Users/alistair/Documents/Berkeley/Levine Lab/Projects/Shadow Enhancers/';
fdir = 'Results-July/';
imn = 'Kr2enh_LacZ_kr_04.tif'; 
I = imread([folder,'/',imn]);

C = [3,0,0 ;
    0,2,0;
    0,.75,.75]; 

T = [.1,.7;
     .12, .45;
     .22, 1];  % min max thresholds for channels
 
 f = [1,0]; % imflip horizontal = 1, vertical = 2, no flip = 0

 Io = im_recolor(I,C,T,f); 
Io = imresize(Io,.3); 

figure(1); clf; imshow(Io); 

imwrite(Io,[prj,fdir,'rc_',imn],'tif');
%%
clear all;
folder =  '/Volumes/Data/Lab Data/Shadow_data';
prj = '/Users/alistair/Documents/Berkeley/Levine Lab/Projects/Shadow Enhancers/';
fdir = 'Results-July/';
imn = 'KrCD1_LacZ_kr_01.tif'; 
I = imread([folder,'/',imn]);

C = [2,0,0 ;
    0,1,0;
    0,.75,.75]; 

T = [.1,.6;
     .12, .7;
     .22, 1];  % min max thresholds for channels
 
 f = [1,0]; % imflip horizontal = 1, vertical = 2, no flip = 0

 Io = im_recolor(I,C,T,f); 
Io = imresize(Io,.3); 

figure(1); clf; imshow(Io); 

imwrite(Io,[prj,fdir,'rc_',imn],'tif');
%%
clear all;
folder =  '/Volumes/Data/Lab Data/Shadow_data';
prj = '/Users/alistair/Documents/Berkeley/Levine Lab/Projects/Shadow Enhancers/';
fdir = 'Results-July/';
imn ='Kr2enh_LacZ_kr_03.tif'; % 'KrCD2_LacZ_kr_01.tif'; 
I = imread([folder,'/',imn]);

C = [0,0,0 ;
    0,2,0;
    .8,0,.0]; 

T = [.1,.8;
     .1, .7;
     .22, 1];  % min max thresholds for channels
 
 f = [0,2]; % imflip horizontal = 1, vertical = 2, no flip = 0

 Io = im_recolor(I,C,T,f); 
Io = imresize(Io,.3); 

figure(1); clf; imshow(Io); 

 imwrite(Io,[prj,fdir,'rc2_',imn],'tif');
%%
clear all;
folder =  '/Volumes/Data/Lab Data/Shadow_data';
prj = '/Users/alistair/Documents/Berkeley/Levine Lab/Projects/Shadow Enhancers/';
fdir = 'Results-July/';
imn = 'KrCD2_LacZ_kr_02.tif'; 
I = imread([folder,'/',imn]);

C = [2,0,0 ;
    0,1,0;
    0,.75,.75]; 

T = [.1,.9;
     .3, .7;
     .22, 1];  % min max thresholds for channels
 
 f = [0,2]; % imflip horizontal = 1, vertical = 2, no flip = 0

 Io = im_recolor(I,C,T,f); 
Io = imresize(Io,.3); 

figure(1); clf; imshow(Io); 

imwrite(Io,[prj,fdir,'rc_',imn],'tif');
