
%%                             ShowEctopNuc.m
% 
% Alistair Boettiger                                   Date Begun: 03/01/11
% Levine Lab                                Functionally Complete: 03/01/11
%                                                   Last Modified: 03/01/11
%

function Im_seg = ShowEctopNuc(folder,name,norm,scale,f,hb)

% name = MP02; norm = 2; scale = .3; f = [1,0];

    Ymax = 250;
    Nstrength = .7;

   load([folder,name]); 
    age = getage(H,cent);
    disp(age);
    
   if hb == 1 
        if f(1) == 1
            Nnucs = max(H(:));
                  Filt = ismember(H,Nnucs-120:Nnucs);
                   % figure(11); clf; imshow(Filt);
        else
            Filt = ismember(H,1:120);
           % figure(11); clf; imshow(Filt);
        end
   else
       Filt = 1; 
   end
    
    
     Iz = uint8(zeros(h,w,3));
     Iz(:,:,1) = 1*uint8(Ymax*L1)+ 255*uint8(L1&(1-L2).*Filt);
     Iz(:,:,2) = 1*uint8(Ymax*(L1)) + Nstrength*handles.In - 255*uint8(L1&(1-L2).*Filt);
     Iz(:,:,3) = Nstrength*handles.In - Iz(:,:,1);
     Im_seg = uint8(bsxfun(@times,double(Iz)/255*norm,double(handles.In)));
     
     Im_seg = imresize(Im_seg,scale);
     Im_seg = imflip(imflip(Im_seg,f(1)),f(2));
     
     
   imshow(Im_seg);
     
     