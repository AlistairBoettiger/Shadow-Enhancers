
%% 

% Alistair Boettiger
% 

clear all;

folder = '/Volumes/Data/Lab Data/Raw_Data/02-17-11/MP08-snaD/';
%fname = 'MP08-snaD_sna_LacZ_sim_02';
fname = 'MP08-snaD_sna_LacZ_sim_04';


fin = [folder,fname,'_max.tif'];
I = imread(fin);

[h,w] = size(I(:,:,1));
%%
Iz = uint8(zeros(h,w,3));
Iz(:,:,1) = I(:,:,4); % snail
Iz(:,:,2) = I(:,:,2); % sim
Iz(:,:,3) = I(:,:,1) ; % nuclei


C = [1.5,0,0;
    1,1,0;
    0,.6,.6];

T = [.06,.6;
    .1,.35;
    .2,.6];

f = [0,0];

Iz = im_recolor(Iz,C,T,f) ;

ICyO = imadjust(I(:,:,3),[.12,.25] ,[0,1]);

Iz(:,:,2) = Iz(:,:,2)+ ICyO;


Io = imresize(Iz,.5); 
figure(1); clf; imshow(Io);





