
%%                             ShowEctopNuc.m
% 
% Alistair Boettiger                                   Date Begun: 03/01/11
% Levine Lab                                Functionally Complete: 03/01/11
%                                                   Last Modified: 03/01/11
%

function Im_seg = ShowEctopNuc(folder,name,norm,scale,f)

% name = kr1;

    Ymax = 250;
    Nstrength = .7;

   load([folder,name]); 
    age = getage(H,cent);
    disp(age);
    
     Iz = uint8(zeros(h,w,3));
     Iz(:,:,1) = 1*uint8(Ymax*L1)+ 255*uint8(L1&(1-L2));
     Iz(:,:,2) = 1*uint8(Ymax*(L1)) + Nstrength*handles.In - 255*uint8(L1&(1-L2));
     Iz(:,:,3) = Nstrength*handles.In - Iz(:,:,1);
     Im_seg = uint8(bsxfun(@times,double(Iz)/255*norm,double(handles.In)));
     
     Im_seg = imresize(Im_seg,scale);
     Im_seg = imflip(imflip(Im_seg,f(1)),f(2));
     
     
   imshow(Im_seg);
     
     