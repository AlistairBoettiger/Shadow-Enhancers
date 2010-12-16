
%% Figs_hb-embs

% Figure plotting of hb embryos for presentation
 
clear all



fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Results/';

%% Boundary vs. central region comparison for cycle 13 embryos
folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed/';

edat09 = 'MP09_22C_y_hb50_data.mat'; % cc13
edat01 = 'BAC01b_22C_y_hb38_data.mat'; % cc13
edat02 = 'BAC02_22C_y_hb08_data.mat'; % cc13  04  06  08


load([folder,edat09]); 
          age = getage(H,cent)         
     [miss_rate1,miss_rate2,miss_rate3] =  divide_regs(L2,H,pts1,pts2,ptr_nucin2,handles.In,0);
    load test2;  
    Ymax = 250;
     imReg3 = ismember(H,onReg3);
     imReg2 = ismember(H,onReg2); 
       Iz = uint8(zeros(h,w,3));
     Iz(:,:,1) = 1*uint8(Ymax*L1)+ 255*uint8(imReg3) - 255*uint8(imReg2);
     Iz(:,:,2) = 1*uint8(Ymax*L1) + .6*handles.In - 255*uint8(imReg3)   + 255*uint8(imReg2) ;
     Iz(:,:,3) = .6*handles.In - Iz(:,:,1) + 205*uint8(imReg2);
     Im_seg09 = uint8(bsxfun(@times,double(Iz)/215,double(handles.In)));
     
     Im_seg09 = imflip(imflip(Im_seg09,2),1);
     figure(1); clf; imshow(Im_seg09);

     
     
     
     
             Ns1 = ismember(H, s1s); % map of cells in region
        Ns2 = ismember(H, s2s);
        Ns3 = ismember(H, s3s);
      
      Iz = uint8(zeros(h,w,3));
     Iz(:,:,1) =255*uint8(Ns3) ;
     Iz(:,:,2) =  .5*handles.In+ 255*uint8(Ns1) + 255*uint8(Ns2);
     Iz(:,:,3) = .5*handles.In + 255*uint8(Ns2);
     reg_schema = uint8(bsxfun(@times,double(Iz)/215,double(handles.In))); 
     reg_schema = imflip(imflip(reg_schema,1),2); 
     
     figure(2); clf; imshow(reg_schema)
     
           imwrite(reg_schema,[fout,'reg_schema.tif'],'tif');
     

     edat01 = 'BAC01b_22C_y_hb38_data.mat'; % cc13
load([folder,edat01]); 
          age = getage(H,cent)         
     [miss_rate1,miss_rate2,miss_rate3] =  divide_regs(L2,H,pts1,pts2,ptr_nucin2,handles.In,0);
    load test2;  
    Ymax = 250;
     imReg3 = ismember(H,onReg3);
     imReg2 = ismember(H,onReg2); 
       Iz = uint8(zeros(h,w,3));
     Iz(:,:,1) = 1*uint8(Ymax*L1)+ 255*uint8(imReg3) - 255*uint8(imReg2);
     Iz(:,:,2) = 1*uint8(Ymax*L1) + 1*handles.In - 255*uint8(imReg3)   + 255*uint8(imReg2) ;
     Iz(:,:,3) = 1*handles.In - Iz(:,:,1) + 205*uint8(imReg2);
     Im_seg01 = uint8(bsxfun(@times,double(Iz)/55,double(handles.In)));
     Im_seg01 = imflip(Im_seg01,2); 
     figure(1); clf; imshow(Im_seg01); 


load([folder,edat02]); 
          age = getage(H,cent)         
     [miss_rate1,miss_rate2,miss_rate3] =  divide_regs(L2,H,pts1,pts2,ptr_nucin2,handles.In,0);
    load test2;  
    Ymax = 250;
     imReg3 = ismember(H,onReg3);
     imReg2 = ismember(H,onReg2); 
       Iz = uint8(zeros(h,w,3));
     Iz(:,:,1) = 1*uint8(Ymax*L1)+ 255*uint8(imReg3) - 255*uint8(imReg2);
     Iz(:,:,2) = 1*uint8(Ymax*L1) + .8*handles.In - 255*uint8(imReg3)   + 255*uint8(imReg2) ;
     Iz(:,:,3) = .8*handles.In - Iz(:,:,1) + 255*uint8(imReg2);
     Im_seg = uint8(bsxfun(@times,double(Iz)/115,double(handles.In)));

     Im_seg02 = imflip(imflip(Im_seg,2),1);
     figure(1); clf; imshow(Im_seg02); 
     
     figure(1); clf; subplot(1,3,1); imshow(Im_seg02); 
     subplot(1,3,2); imshow(Im_seg01);
     subplot(1,3,3); imshow(Im_seg09); set(gcf,'color','k'); 
     
      imwrite(Im_seg01,[fout,'Im_seg01.tif'],'tif');
        imwrite(Im_seg02,[fout,'Im_seg02.tif'],'tif');
        imwrite(Im_seg09,[fout,'Im_seg09.tif'],'tif');
      
      

%% Temperature comparison for cycle 14 boundary embryos 
% MP01 30C vs 22C
folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed/';


edat14 = 'BAC01b_30C_y_hb18_data.mat'; % cc14

edat14 = 'BAC01b_30C_y_hb15_data.mat'; % cc14
load([folder,edat14]); 
          age = getage(H,cent)
          
     [miss_rate1,miss_rate2,miss_rate3] =  divide_regs(L2,H,pts1,pts2,ptr_nucin2,handles.In,0);
    load test2;  
     imReg3 = ismember(H,onReg3);
       Iz = uint8(zeros(h,w,3));
     Iz(:,:,1) = 1*uint8(200*L1)+ 255*uint8(imReg3);
     Iz(:,:,2) = 1*uint8(200*L1) + .8*handles.In - 255*uint8(imReg3);
     Iz(:,:,3) = .8*handles.In - Iz(:,:,1);
     B1_30C = uint8(bsxfun(@times,double(Iz)/175,double(handles.In)));

        
    % edat14 = 'BAC01b_22C_y_hb27_data.mat'; % cc14
     % edat14 =   'MP01_22C_y_hb06_data.mat';
     edat14 =   'MP01_22C_y_hb07_data.mat';
     load([folder,edat14]);
          age = getage(H,cent)
       divide_regs(L2,H,pts1,pts2,ptr_nucin2,handles.In,0);
    load test2;     
     imReg3 = ismember(H,onReg3);     
       Iz = uint8(zeros(h,w,3));
     Iz(:,:,1) = 1*uint8(200*L1)+ 255*uint8(imReg3);
     Iz(:,:,2) = 1*uint8(200*L1) + .8*handles.In - 255*uint8(imReg3);
     Iz(:,:,3) = .8*handles.In - Iz(:,:,1);
     B1_22C = uint8(bsxfun(@times,double(Iz)/215,double(handles.In)));
     B1_22C = imflip(B1_22C,1); 

          figure(2); subplot(2,2,2); imshow(B1_22C);
            figure(2); subplot(2,2,4); imshow(B1_30C);
%

% MP02 

%edat14 = 'MP02_30C_LacZ_hb29_data.mat'; % cc14
edat14 = 'MP02_30C_LacZ_hb10_data.mat'; % cc14 30-good 44 43 42
load([folder,edat14]); 
     [miss_rate1,miss_rate2,miss_rate3] =  divide_regs(L2,H,pts1,pts2,ptr_nucin2,handles.In,0);
    load test2;  
     imReg3 = ismember(H,onReg3);
       Iz = uint8(zeros(h,w,3));
     Iz(:,:,1) = 1*uint8(200*L1)+ 255*uint8(imReg3);
     Iz(:,:,2) = 1*uint8(200*L1) + .8*handles.In - 255*uint8(imReg3);
     Iz(:,:,3) = .8*handles.In - Iz(:,:,1);
     B2_30C = uint8(bsxfun(@times,double(Iz)/215,double(handles.In)));
     B2_30C = imflip(B2_30C,1);
        
     
   %  edat14 = 'BAC02_22C_y_hb11_data.mat'; % cc14
     edat14 = 'MP02_22C_y_hb57_data.mat'; % cc14 % 77 81 57
     load([folder,edat14]); 
     [miss_rate1,miss_rate2,miss_rate3] =  divide_regs(L2,H,pts1,pts2,ptr_nucin2,handles.In,0);
    load test2;     
     imReg3 = ismember(H,onReg3);     
       Iz = uint8(zeros(h,w,3));
     Iz(:,:,1) = 1*uint8(200*L1)+ 255*uint8(imReg3);
     Iz(:,:,2) = 1*uint8(200*L1) + .8*handles.In - 255*uint8(imReg3);
     Iz(:,:,3) = .8*handles.In - Iz(:,:,1);
     B2_22C = uint8(bsxfun(@times,double(Iz)/165,double(handles.In)));
     

     
     
     
     
     
     
     
 %% MP09 boundary region
 
      
edat14 = 'MP09_22C_y_hb01_data.mat'; % cc14
     load([folder,edat14]); 
     [miss_rate1,miss_rate2,miss_rate3] =  divide_regs(L2,H,pts1,pts2,ptr_nucin2,handles.In,0);
    load test2;     
     imReg3 = ismember(H,onReg3);     
       Iz = uint8(zeros(h,w,3));
     Iz(:,:,1) = 1*uint8(200*L1)+ 255*uint8(imReg3);
     Iz(:,:,2) = 1*uint8(200*L1) + .8*handles.In - 255*uint8(imReg3);
     Iz(:,:,3) = .8*handles.In - Iz(:,:,1);
     BC_22C = uint8(bsxfun(@times,double(Iz)/205,double(handles.In)));
     BC_22C = imflip(BC_22C,1); 
 
edat14 = 'BAC09_30C_y_hb07_data.mat'; % cc14
load([folder,edat14]); 
     [miss_rate1,miss_rate2,miss_rate3] =  divide_regs(L2,H,pts1,pts2,ptr_nucin2,handles.In,0);
    load test2;  
     imReg3 = ismember(H,onReg3);
       Iz = uint8(zeros(h,w,3));
     Iz(:,:,1) = 1*uint8(200*L1)+ 255*uint8(imReg3);
     Iz(:,:,2) = 1*uint8(200*L1) + .8*handles.In - 255*uint8(imReg3);
     Iz(:,:,3) = .8*handles.In - Iz(:,:,1);
     BC_30C = uint8(bsxfun(@times,double(Iz)/185,double(handles.In)));
    BC_30C = imflip(imflip(BC_30C,1),2); 
        



          figure(2); subplot(2,3,3); imshow(BC_22C);
            figure(2); subplot(2,3,6); imshow(BC_30C);
 %%        Plot and write to disk   
     
       
     
     
          figure(2); clf; subplot(2,3,1); imshow(B2_22C);
     figure(2); subplot(2,3,4); imshow(B2_30C);

     
     
          figure(2); subplot(2,3,2); imshow(B1_22C);
            figure(2); subplot(2,3,5); imshow(B1_30C);
            
          figure(2); subplot(2,3,3); imshow(BC_22C);
            figure(2); subplot(2,3,6); imshow(BC_30C);
            
            set(gcf,'color','k');
     
    imwrite(B1_22C,[fout,'B1_22C_bndry.tif'],'tif');
      imwrite(B1_30C,[fout,'B1_30C_bndry.tif'],'tif');
      imwrite(B2_22C,[fout,'B2_22C_bndry.tif'],'tif');
      imwrite(B2_30C,[fout,'B2_30C_bndry.tif'],'tif');
            
       imwrite(BC_22C,[fout,'BC_22C_bndry.tif'],'tif');
      imwrite(BC_30C,[fout,'BC_30C_bndry.tif'],'tif');   
            
%% MP01 30C
folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed/';
edat14 = 'BAC01b_30C_y_hb18_data.mat'; % cc14
edat13 = 'BAC01b_30C_y_hb03_data.mat'; % cc13 
edat12 = 'BAC01b_30C_y_hb04_data.mat'; % cc12
edat11 = 'BAC01_30C_y_hb08_data.mat'; % cc11

load([folder,edat14]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) = 1*uint8(200*L1) + .8*handles.In;
     Io(:,:,3) = .8*handles.In - Io(:,:,1);
     Icc14 = uint8(bsxfun(@times,double(Io)/75,double(handles.In)));
     Icc14 = imflip(imflip(Icc14,1),2);
     age = getage(H,cent)
           

    
     
     load([folder,edat13]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) = 1*uint8(200*L1) + .8*handles.In;
     Io(:,:,3) = .8*handles.In - Io(:,:,1);
     Icc13 = uint8(bsxfun(@times,double(Io)/155,double(handles.In)));
     Icc13 = imflip(Icc13,2);
     age = getage(H,cent)
     

    
    
     
     
     load([folder,edat12]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) = 1*uint8(200*L1) + .8*handles.In;
     Io(:,:,3) = .8*handles.In - Io(:,:,1);
     Icc12 = uint8(bsxfun(@times,double(Io)/155,double(handles.In)));
     Icc12 = imflip(Icc12,2); 
     age = getage(H,cent)
     
          load([folder,edat11]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) = 1*uint8(200*L1) + .8*handles.In;
     Io(:,:,3) = .8*handles.In - Io(:,:,1);
     Icc11 = uint8(bsxfun(@times,double(Io)/155,double(handles.In)));
     age = getage(H,cent)
     
     figure(1); clf;
     subplot(4,1,1); imshow(Icc14);
     subplot(4,1,2); imshow(Icc13);
     subplot(4,1,3); imshow(Icc12);
       subplot(4,1,4); imshow(Icc11);
     set(gcf,'color','k');
     
%      imwrite(Icc14,[fout,'BAC01_30C_cc14.tif'],'tif');
%      imwrite(Icc13,[fout,'BAC01_30C_cc13.tif'],'tif');
%      imwrite(Icc12,[fout,'BAC01_30C_cc12.tif'],'tif');
%      imwrite(Icc11,[fout,'BAC01_30C_cc11.tif'],'tif');
%      

%% MP 01 22C
edat14 = 'BAC01b_22C_y_hb27_data.mat'; % cc14
edat13 = 'BAC01b_22C_y_hb38_data.mat'; % cc13
%edat11 = 'BAC01b_22C_y_hb41_data.mat'; % cc11
edat11 = 'BAC01b_22C_y_hb39_data.mat'; % cc11
edat12 = 'BAC01b_22C_y_hb35_data.mat'; % cc12
edat10 = 'BAC01b_22C_y_hb40_data.mat'; % cc10
% edat10 = 'BAC01b_22C_y_hb33_data.mat'; % cc10


load([folder,edat14]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) = 1*uint8(200*L1) + .8*handles.In;
     Io(:,:,3) = .8*handles.In - Io(:,:,1);
     Icc14 = uint8(bsxfun(@times,double(Io)/200,double(handles.In)));

     age = getage(H,cent)
               
   
     
     
     load([folder,edat13]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1) ;
     Io(:,:,2) = 1*uint8(200*L1) + handles.In;
     Io(:,:,3) = handles.In - Io(:,:,1);
     Icc13 = uint8(bsxfun(@times,double(Io)/55,double(handles.In)));
     Icc13 = imflip(Icc13,2);
     age = getage(H,cent)
     
     load([folder,edat12]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) = 1*uint8(200*L1) + handles.In;
     Io(:,:,3) = handles.In - Io(:,:,1);
     Icc12 = uint8(bsxfun(@times,double(Io)/55,double(handles.In)));
     Icc12 = imflip(imflip(Icc12,1),2);
     age = getage(H,cent)
     
          load([folder,edat11]); 
      Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) = 1*uint8(200*L1) + handles.In;
     Io(:,:,3) = handles.In - Io(:,:,1);
     Icc11 = uint8(bsxfun(@times,double(Io)/25,double(handles.In)));
     Icc11 = imflip(imflip(Icc11,1),2);
     age = getage(H,cent)
           
     load([folder,edat10]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) = 1*uint8(200*L1) + handles.In;
     Io(:,:,3) = handles.In - Io(:,:,1);
     Icc10 = uint8(bsxfun(@times,double(Io)/15,double(handles.In)));
     Icc10 = imflip(Icc10,2);
     age = getage(H,cent)
     
     
     figure(2); clf;
     subplot(5,1,1); imshow(Icc14);
     subplot(5,1,2); imshow(Icc13);
     subplot(5,1,3); imshow(Icc12);
     subplot(5,1,4); imshow(Icc11);
     subplot(5,1,5); imshow(Icc10);
 set(gcf,'color','k');

      
     imwrite(Icc14,[fout,'BAC01_22C_cc14.tif'],'tif');
     imwrite(Icc13,[fout,'BAC01_22C_cc13.tif'],'tif');
     imwrite(Icc12,[fout,'BAC01_22C_cc12.tif'],'tif');
     imwrite(Icc11,[fout,'BAC01_22C_cc11.tif'],'tif');
        imwrite(Icc10,[fout,'BAC01_22C_cc10.tif'],'tif');
 
 %% MP 02 22C
edat14 = 'BAC02_22C_y_hb11_data.mat'; % cc14
edat13 = 'BAC02_22C_y_hb10_data.mat'; % cc13
%edat13 = 'MP02_22C_y_hb39_data.mat'; % cc13
edat12 = 'MP02_22C_y_hb46_data.mat'; % cc12
edat12 = 'BAC02_22C_y_hb01_data.mat'; % cc11  % actually cc12

% edat10 = 'BAC01b_22C_y_hb33_data.mat'; % cc10


load([folder,edat14]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) = 1*uint8(200*L1) + .8*handles.In;
     Io(:,:,3) = .8*handles.In - Io(:,:,1);
     Icc14 = uint8(bsxfun(@times,double(Io)/55,double(handles.In)));

     age = getage(H,cent)
           
     load([folder,edat13]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1) ;
     Io(:,:,2) = 1*uint8(200*L1) + handles.In;
     Io(:,:,3) = handles.In - Io(:,:,1);
     Icc13 = uint8(bsxfun(@times,double(Io)/55,double(handles.In)));
          age = getage(H,cent)
     
     load([folder,edat12]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) = 1*uint8(200*L1) + handles.In;
     Io(:,:,3) = handles.In - Io(:,:,1);
     Icc12 = uint8(bsxfun(@times,double(Io)/135,double(handles.In)));
     age = getage(H,cent)
     
          load([folder,edat11]); 
      Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) = 1*uint8(200*L1) + handles.In;
     Io(:,:,3) = handles.In - Io(:,:,1);
     Icc11 = uint8(bsxfun(@times,double(Io)/125,double(handles.In)));
     age = getage(H,cent)
           
     
     
    figure(1); clf;
     subplot(4,1,1); imshow(Icc14);
     subplot(4,1,2); imshow(Icc13);
     subplot(4,1,3); imshow(Icc12);
       subplot(4,1,4); imshow(Icc11);
     set(gcf,'color','k');
     imwrite(Icc14,[fout,'BAC02_22C_cc14.tif'],'tif');
     imwrite(Icc13,[fout,'BAC02_22C_cc13.tif'],'tif');
     imwrite(Icc12,[fout,'BAC02_22C_cc12.tif'],'tif');
     imwrite(Icc11,[fout,'BAC02_22C_cc11.tif'],'tif');
      %% MP 02 30C
edat14 = 'MP02_30C_LacZ_hb29_data.mat'; % cc14
edat13 = 'MP02_30C_LacZ_hb09_data.mat'; % cc13
edat12 = 'MP02_30C_LacZ_hb22_data.mat'; % cc12
%edat12 = 'MP02_30C_LacZ_hb06_data.mat'; % cc12
edat11 = 'MP02_30C_LacZ_hb21_data.mat'; % cc11  % actually cc12

% edat10 = 'BAC01b_22C_y_hb33_data.mat'; % cc10


load([folder,edat14]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) = 1*uint8(200*L1) + .8*handles.In;
     Io(:,:,3) = .8*handles.In - Io(:,:,1);
     Icc14 = uint8(bsxfun(@times,double(Io)/155,double(handles.In)));
    Icc14 = imflip(Icc14,2); 
     age = getage(H,cent)
           
     load([folder,edat13]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1) ;
     Io(:,:,2) = 1*uint8(200*L1) + handles.In;
     Io(:,:,3) = handles.In - Io(:,:,1);
     Icc13 = uint8(bsxfun(@times,double(Io)/135,double(handles.In)));
     Icc13 = imflip(imflip(Icc13,1),2);
          age = getage(H,cent)
     
     load([folder,edat12]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) = 1*uint8(200*L1) + handles.In;
     Io(:,:,3) = handles.In - Io(:,:,1);
     Icc12 = uint8(bsxfun(@times,double(Io)/135,double(handles.In)));
     Icc12 = imflip(Icc12,1); 
     age = getage(H,cent)
     
          load([folder,edat11]); 
      Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) = 1*uint8(200*L1) + handles.In;
     Io(:,:,3) = handles.In - Io(:,:,1);
     Icc11 = uint8(bsxfun(@times,double(Io)/125,double(handles.In)));
     Icc11 = imflip(Icc11,1); 
     age = getage(H,cent)
           
     
     
    figure(1); clf;
     subplot(4,1,1); imshow(Icc14);
     subplot(4,1,2); imshow(Icc13);
     subplot(4,1,3); imshow(Icc12);
       subplot(4,1,4); imshow(Icc11);
     set(gcf,'color','k');
     
     imwrite(Icc14,[fout,'BAC02_30C_cc14.tif'],'tif');
     imwrite(Icc13,[fout,'BAC02_30C_cc13.tif'],'tif');
     imwrite(Icc12,[fout,'BAC02_30C_cc12.tif'],'tif');
     imwrite(Icc11,[fout,'BAC02_30C_cc11.tif'],'tif');
       %% MP 09 22C
edat14 = 'MP09_22C_y_hb01_data.mat'; % cc14
edat13 = 'MP09_22C_y_hb50_data.mat'; % cc13
edat12 = 'MP09_22C_y_hb51_data.mat'; % cc12
%edat12 = 'MP02_30C_LacZ_hb06_data.mat'; % cc12
edat11 = 'MP09_22C_y_hb34_data.mat'; % cc11  % actually cc12

% edat10 = 'BAC01b_22C_y_hb33_data.mat'; % cc10


load([folder,edat14]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) = 1*uint8(200*L1) + .8*handles.In;
     Io(:,:,3) = .8*handles.In - Io(:,:,1);
     Icc14 = uint8(bsxfun(@times,double(Io)/200,double(handles.In)));
    Icc14 = imflip(Icc14,1); 
     age = getage(H,cent)
           
     load([folder,edat13]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1) ;
     Io(:,:,2) = 1*uint8(200*L1) + handles.In;
     Io(:,:,3) = handles.In - Io(:,:,1);
     Icc13 = uint8(bsxfun(@times,double(Io)/200,double(handles.In)));
     Icc13 = imflip(imflip(Icc13,1),2);
          age = getage(H,cent)
     
     load([folder,edat12]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) = 1*uint8(200*L1) + handles.In;
     Io(:,:,3) = handles.In - Io(:,:,1);
     Icc12 = uint8(bsxfun(@times,double(Io)/135,double(handles.In)));
     Icc12 = imflip(Icc12,1); 
     age = getage(H,cent)
     
          load([folder,edat11]); 
      Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) = 1*uint8(200*L1) + handles.In;
     Io(:,:,3) = handles.In - Io(:,:,1);
     Icc11 = uint8(bsxfun(@times,double(Io)/155,double(handles.In)));
     Icc11 = imflip(Icc11,1); 
     age = getage(H,cent)
           
      
    figure(1); clf;
     subplot(4,1,1); imshow(Icc14);
     subplot(4,1,2); imshow(Icc13);
     subplot(4,1,3); imshow(Icc12);
       subplot(4,1,4); imshow(Icc11);
     set(gcf,'color','k');
     imwrite(Icc14,[fout,'BAC09_22C_cc14.tif'],'tif');
     imwrite(Icc13,[fout,'BAC09_22C_cc13.tif'],'tif');
     imwrite(Icc12,[fout,'BAC09_22C_cc12.tif'],'tif');
     imwrite(Icc11,[fout,'BAC09_22C_cc11.tif'],'tif');
     
     
            %% MP 09 30C
edat14 = 'BAC09_30C_y_hb07_data.mat'; % cc14
edat13 = 'BAC09_30C_y_hb71_data.mat'; % cc13
edat12 = 'BAC09_30C_y_hb65_data.mat'; % cc12
%edat12 = 'MP02_30C_LacZ_hb06_data.mat'; % cc12
edat11 = 'BAC09_30C_y_hb64_data.mat'; % cc11  % actually cc12

% edat10 = 'BAC01b_22C_y_hb33_data.mat'; % cc10


load([folder,edat14]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) = 1*uint8(200*L1) + .8*handles.In;
     Io(:,:,3) = .8*handles.In - Io(:,:,1);
     Icc14 = uint8(bsxfun(@times,double(Io)/125,double(handles.In)));
    Icc14 = imflip(imflip(Icc14,2),1); 
     age = getage(H,cent)
           
     load([folder,edat13]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1) ;
     Io(:,:,2) = 1*uint8(200*L1) + handles.In;
     Io(:,:,3) = handles.In - Io(:,:,1);
     Icc13 = uint8(bsxfun(@times,double(Io)/65,double(handles.In)));
     Icc13 =imflip(Icc13,1);
          age = getage(H,cent)
     
     load([folder,edat12]); 
   Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) = 1*uint8(200*L1) + handles.In;
     Io(:,:,3) = handles.In - Io(:,:,1);
     Icc12 = uint8(bsxfun(@times,double(Io)/55,double(handles.In)));
     age = getage(H,cent)
     
          load([folder,edat11]); 
      Io = uint8(zeros(h,w,3));
     Io(:,:,1) = 1*uint8(200*L1);
     Io(:,:,2) = 1*uint8(200*L1) + handles.In;
     Io(:,:,3) = handles.In - Io(:,:,1);
     Icc11 = uint8(bsxfun(@times,double(Io)/55,double(handles.In)));
     age = getage(H,cent)
           
      
    figure(1); clf;
     subplot(4,1,1); imshow(Icc14);
     subplot(4,1,2); imshow(Icc13);
     subplot(4,1,3); imshow(Icc12);
       subplot(4,1,4); imshow(Icc11);
     set(gcf,'color','k');
     imwrite(Icc14,[fout,'BAC09_30C_cc14.tif'],'tif');
     imwrite(Icc13,[fout,'BAC09_30C_cc13.tif'],'tif');
     imwrite(Icc12,[fout,'BAC09_30C_cc12.tif'],'tif');
     imwrite(Icc11,[fout,'BAC09_30C_cc11.tif'],'tif');
     
%% exiting mitosis MP09 at 30C


folder = '/Volumes/Data/Lab Data/Shadow_data';
emb = 'BAC09_30C_y_hb_11.tif';
name = 'BAC09_30C_y_hb_11.tif';

I = imread([folder,'/',emb]);  

C = [1,1,0;
    0,0,0;
    0,.7,.7];

T = [.13,.5;
    .17,1;
    .1,.9];

f = [1,2];

I = im_recolor(I,C,T,f) ;
figure(1); clf; imshow(I);

Io = imresize(I,1);


imwrite(Io,[fout,name]);

%% somewhat Young synchronous MP09 at 30C


folder = '/Volumes/Data/Lab Data/Shadow_data';
emb = 'BAC09_30C_y_hb_47.tif';
name = 'BAC09_30C_y_hb_47.tif';

I = imread([folder,'/',emb]);  

C = [1,1,0;
    2,0,0;
    0,.7,.7];

T = [.2,1;
    .17,1;
    .1,.9];

f = [1,2];

I = im_recolor(I,C,T,f) ;
figure(1); clf; imshow(I);

Io = imresize(I,1);


imwrite(Io,[fout,name]);


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
     Io(:,:,1) = 2*uint8(200*L1);
     Io(:,:,2) = 2*uint8(200*L1) + 1*handles.In;
     Io(:,:,3) = 1*handles.In - Io(:,:,1);
     Ired = uint8(bsxfun(@times,double(Io)/200,double(handles.In)));
     

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


