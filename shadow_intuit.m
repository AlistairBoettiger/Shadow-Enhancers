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



load  hbCorrData3;% hbCorrData2 good run


%~~~~~~~~~~~~~~  Create null model with random probability of being on
frac_on = .65; % length(nuc_on1)/length((ptr_nucin1));
primary = find(rand(1,N_tot)<frac_on);
P = ismember(H,primary);   

shadow = find(rand(1,N_tot)<frac_on);

S = ismember(H,shadow);


     
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
name = 'stoch_Primary.tif'; 
imwrite(I,[fout,name]);


 
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

 name ='stoch_shadow.tif'; 
imwrite(I,[fout,name]);
 

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
 name = 'stoch_2enh.tif'; 
imwrite(I,[fout,name]);


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