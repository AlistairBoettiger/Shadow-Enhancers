%%                              shadow_intuit.m
% Alistair Boettiger                                   Date Begun: 10/09/10
%                                                   Last Modified: 10/08/10
%

% Simulate independent enhancer function in a whole embryo
% Take an embryo nuclear map, randomly assign a fraction of cells as on,
% randomly assign all cells a second time with same probability of being
% on.  Combine the two by finding all cells on in either the first or the
% second instance. 
%
% Plot all the results and color shift to yellow-cyan.

clear all; close all;

 folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed'; 
 base_emb = 'SogP_C81_22C_LacZ_sog11_data.mat';
load([folder,'/',base_emb]); 
 
%load  hbCorrData3;% hbCorrData2 good run

N = 20; % rounds of transcription
%~~~~~~~~~~~~~~  Create null model with random probability of being on
frac_on = .5; % length(nuc_on1)/length((ptr_nucin1));

N_tot = max(H(:));

r = rand(1,N_tot); 
c = .8;

    primary = find(r<frac_on);
      P(:,:,1) = ismember(H,primary).*Reg1;  
          T =  P(:,:,1);
   figure(1); clf;
     Io = uint8(zeros(h,w,3));
     Io(:,:,1) = uint8(255*T);
     Io(:,:,3) = 30*handles.In;
     Ion = uint8(bsxfun(@times,double(Io)/100,double(handles.In)));
     Ion = imresize(Ion,.4); 
     imshow(Ion); hold on;
     pause(.01);
   
    movieframe = figure(2); clf;
       Io = uint8(zeros(h,w,3));
     Io(:,:,3) = 30*handles.In;
     Ion = uint8(bsxfun(@times,double(Io)/100,double(handles.In)));
        Ion = imresize(Ion,.3); imshow(Ion); hold on;
     axis tight
set(gca,'nextplot','replacechildren');
F = getframe; 
    

P = uint8(zeros(h,w,N));
for k=2:N
    
    r = r*c + rand(1,N_tot)*(1-c);
    primary = find(r<frac_on);
    P(:,:,k) = ismember(H,primary).*Reg1;   
    T = sum(P,3);

    movieframe = figure(2); clf;
     Io = uint8(zeros(h,w,3));
     Io(:,:,1) = uint8(255*T/(N*frac_on));
     Io(:,:,3) = 30*handles.In;
     Ion = uint8(bsxfun(@times,double(Io)/100,double(handles.In)));
        Ion = imresize(Ion,.3); imshow(Ion); hold on;
     pause(.01);
     F(k) = getframe;
end
figure(2); 
 
movie(F,1); 

     fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Results/';
     moviename = 'flicker_corr_ave.avi'

movie2avi (F, [fout,moviename]); 

%%



     
     fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Results/';
         I = uint8(zeros(h,w,3));
         I(:,:,3) = handles.In;       
     I(:,:,1) = uint8(255*P);
I = uint8(bsxfun(@times,double(I)/255,double(handles.In)));

C = [.81,.51,0;
    0,0,0;
    0,.61,.61];
T = [0,1;
    0,1;
    0,1];
f = [1,0];
I = im_recolor(I,C,T,f) ;
 figure(1); clf; imshow(I);



 
         I = uint8(zeros(h,w,3));
         I(:,:,3) = handles.In;     
     I(:,:,1) = uint8(255*S);
I = uint8(bsxfun(@times,double(I)/255,double(handles.In)));
C = [.81,.51,0;
    0,0,0;
    0,.61,.61];
T = [0,1;
    0,1;
    0,1];
f = [1,0];
I = im_recolor(I,C,T,f) ;
 figure(1); clf; imshow(I);

 

         I = uint8(zeros(h,w,3));
         I(:,:,3) = handles.In;     
     I(:,:,1) = uint8(255*(P|S));
I = uint8(bsxfun(@times,double(I)/255,double(handles.In)));
C = [.81,.51,0;
    0,0,0;
    0,.61,.61];
T = [0,1;
    0,1;
    0,1];
f = [1,0];
I = im_recolor(I,C,T,f) ;
 figure(1); clf; imshow(I);



%%
     
     
    figure(6); clf; 
   subplot(3,1,1);
     Io = uint8(zeros(h,w,3));
     Io(:,:,1) = uint8(255*P);
    %  Io(:,:,2) = uint8(255*P);
     Io(:,:,3) = 30*handles.In;
     Ion = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
     imshow(Ion); hold on;
      
     subplot(3,1,2);
     Io = uint8(zeros(h,w,3));
     Io(:,:,1) = uint8(255*S);
    %    Io(:,:,2) = uint8(255*S);
     Io(:,:,3) = 30*handles.In;
     Ion = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
     imshow(Ion); hold on;
     
      subplot(3,1,3);
     Io = uint8(zeros(h,w,3));
     Io(:,:,1) = uint8(255*(P|S)) ;
    %    Io(:,:,2) = uint8(255*(P|S));
        Io(:,:,3) = 30*handles.In;
     Ion = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
     imshow(Ion); hold on;