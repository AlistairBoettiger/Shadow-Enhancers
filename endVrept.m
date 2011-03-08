
%%                             endVrept.m
% 
% Alistair Boettiger                                   Date Begun: 03/01/11
% Levine Lab                                Functionally Complete: 03/01/11
%                                                   Last Modified: 03/01/11
%

function Im_seg = endVrept(folder,name,norm,scale,f)



    Ymax = 250;
    Nstrength = .7;

   load([folder,name]); 
    age = getage(H,cent);
    disp(age);
    
    %N_y = length(unique(H.*(L1&L2)));
    N_y = length(unique(H.*L2n1)); 
    N_g = length(unique(   H.*L2    ));
    N_miss = N_y/N_g; 
    % N_miss = 1- N_y/(N_g);
    
    
   % Yellow = both;  Green = endogenous no reporter 
     Iz = uint8(zeros(h,w,3));
     Iz(:,:,1) =  1*uint8(Ymax*(L1&L2))  ;
     Iz(:,:,2) =  255*uint8(L2n1) + 1*uint8(Ymax*(L1&L2))  + Nstrength*handles.In ;
     Iz(:,:,3) = Nstrength*handles.In - Iz(:,:,1)  - 255*uint8(L2n1) - 1*uint8(Ymax*(L1&L2));
     Im_seg = uint8(bsxfun(@times,double(Iz)/255*norm,double(handles.In)));
     
     Im_seg = imresize(Im_seg,scale);
     Im_seg = imflip(imflip(Im_seg,f(1)),f(2));
     
     
   imshow(Im_seg);
     
     
   
   %  Yellow = both, red = missing
%         Iz = uint8(zeros(h,w,3));
%      Iz(:,:,1) = 1*uint8(Ymax*(L1&L2))+ 255*uint8(L2n1);
%      Iz(:,:,2) = 1*uint8(Ymax*(L1&L2)) + Nstrength*handles.In - 255*uint8(L2n1);
%      Iz(:,:,3) = Nstrength*handles.In - Iz(:,:,1);
%      Im_seg = uint8(bsxfun(@times,double(Iz)/255*norm,double(handles.In)));

% Green = both, yellow = missing

%      Iz = uint8(zeros(h,w,3));
%      Iz(:,:,1) = 255*uint8(L2n1);
%      Iz(:,:,2) = 1*uint8(Ymax*(L1&L2)) + 255*uint8(L2n1) + Nstrength*handles.In ;
%      Iz(:,:,3) = Nstrength*handles.In - Iz(:,:,1)  - 255*uint8(L2n1) - 1*uint8(Ymax*(L1&L2));
%      Im_seg = uint8(bsxfun(@times,double(Iz)/255*norm,double(handles.In)));
%   title(['Missing Frac:',num2str(N_miss,3),'  ',...
%        num2str(N_y),' rep. AND endog.  ',num2str(N_g),' endgo. NOT rep.'],'Color','w','FontSize',10);
   