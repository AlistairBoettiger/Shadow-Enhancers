

clear all;

% im = 'APtoll_dlGFP_Series009_t00_z0';
% im = 'APtoll_twi_sna_sim_07_z0';

% chn = 0; chn2 = 1;
%folder = '/07.01.09/APtoll_dlGFP/';
%folder = '/11.13.08/APtoll/';
% folder =  '/04-16-10/12.2_snaD_CyO-HbZ_29C/';
% im = '12.2_snaD_CyO-HbZ_29C_Bgal-sim_sna_14_z0';
% T = 120; 
% 
% folder =  '/04-16-10/7.5_snaD_CyO-HbZ_29C/';
% im = '7.5_snaD_CyO-HbZ_29C_Bgal-sim_sna_14_z0';
%T = 30; 

% % good 
% folder = '/04-23-10/7.5_snaD_CyO-HbZ_22C/';
% im = '7.5_snaD_CyO-HbZ_22C_Bgal-sim_sna_3D_02_z';  chn = 0; chn2 = 1;% 
% % im = '7.5_snaD_CyO-HbZ_22C_sna_nuc_3D_01_z';

% folder = '/04-23-10/rhoN-vnN_TM3-HbZ_22C/';
% im = 'rhoN-vnN_TM3-HbZ_22C_rho_LacZ-tup_01_z';

% % nice
% folder = '/04-30-10/snaE_snaN/';
% im = 'snaE_snaN_sna_sim_stack_02_z';  chn2 = 2; chn = 0; 
% % im = '7.5_snaD_CyO-HbZ_22C_sna_nuc_3D_01_z';

% nice
folder = '/04-30-10/snaE_snaN/';
im = 'snaE_snaN_sna_sim_stack_02_z';  chn2 = 2; chn = 0; 

% folder = '/05-07-10/BAC7_snaD_HbZ_29C/';
% im = 'BAC7_snaD_HbZ_29C_sim-exon_sna_st_08_z';
% chn2 = 2; chn = 1; L = 60;
% 
% 
% folder = '/05-07-10/BAC12_snaD_HbZ_29C/';
% im = 'BAC12_snaD_HbZ_29C_sim-exon_sna_st_01_z';

% folder = '/05-07-10/BAC7_snaD_HbZ_22C/';
% im = 'BAC7_snaD_HbZ_22C_sim-exon_sna_st_06_z';
% chn2 = 2; chn = 1; L = 60;
% 
% folder = '/05-07-10/BAC7_snaD_HbZ_29C_c/';
% im = 'BAC7_snaD_HbZ_29C_c_vnd_sna-LacZ_st_01_z';
% chn2 = 2; chn = 1; L = 115;


% 150
% 400 approx nuclear
res = 150; % 320;



I = cell(1,20); I2 = cell(1,20); 
watcher = 0;  N = 0;

for z=1:300
  
    if z<11
        zs = ['00', num2str(z-1)];
    elseif z<101
        zs = ['0',num2str(z-1)]; 
    elseif z>=101
        zs = num2str(z-1);
    end
      fn = [im,zs,'_ch0',num2str(chn),'.tif'];
      fn2 = [im,zs,'_ch0',num2str(chn2),'.tif'];
    
    try
    I{z}=imread(['/Volumes/Data/Lab Data/Raw_Data',folder,fn]); 
    I2{z} = imread(['/Volumes/Data/Lab Data/Raw_Data',folder,fn2]); 
    
    catch E
        watcher = 1; 
        disp(E.message); 
    end
    if watcher == 1;
        break;
    end
      N = N+1;  
end

V = uint8(zeros(res,res,N)); 
V2 = V; 
    
for z=1:N
    V(:,:,z) = imresize(I{z},[res,res]); 
    V(:,:,z) = imdilate(V(:,:,z),strel('disk',2)); 
    V2(:,:,z) = imresize(I2{z},[res,res]);
    V2(:,:,z) = imdilate(V2(:,:,z),strel('disk',2)); 
  
    figure(1); clf; imshow(V(:,:,z)); pause(.02);
    
end


for z=1:N
    V3 = uint8(zeros(res,res,3));
    V3(:,:,1) = V(:,:,z);
    V3(:,:,3) = V2(:,:,z);
    figure(1); clf; imshow(V3); pause(.02);
end

[xs,ys,zs] = size(V);

[X,Y,Z] = meshgrid(1:xs,1:ys,1:zs);

%%
L = 55;

figure(3); clf; 
%      p = patch(isosurface(X,Y,Z,V2,100)); hold on; % 130
%      set(p,'FaceColor','blue','EdgeColor','none','FaceAlpha',.2);
   p2 = patch(isosurface(X,Y,Z,V,L)); %55 40
   set(p2,'FaceColor','red','EdgeColor','none');
   set(p2,'FaceAlpha',.95);

lighting none; camlight(180,260);

view(-145,-56);  %camlight(140,-68);

% camlight(100,140);  
 % view(99,30);
% view(101,50); 
lighting phong;  material dull;

% view(13,90); 

 % view(129,34);

 %view(214,72);
% camlight(180,-60);
% %view(-189,18);
 
 % figure(4); clf; hist(double(V(:)),255);

% view(-74,-70);  