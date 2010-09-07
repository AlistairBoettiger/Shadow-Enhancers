
%% anlz_young hb


% Alistair Boettiger                            Date Begun: 09/02/10
% Levine Lab                                    Last Modified: 0


% MP02_30C_LacZ
     %                                       *
young = [06,19,39,49,41,37,35,34,26,24,22,21,14,15,07,01]; % 14
% y2 11,15,08;    

% colors flipped
MP02_22C_y = [46, 79,88,82,76,39,33,25,22,16,15,11,10];

MP09 = [05,09,12,16,22,23,26,28,31,33,34,36,37,38,42,46,48,49,50,51,47,45,27,11,32 ];  % 23,34,49
%34=23

folder =  '/Volumes/Data/Lab Data/Shadow_data/';

emb_roots{z}


for k=1:length(MP09);
    emb = num2str(MP09(k));
    if MP09(k) < 10
        emb = ['0',emb];
    end
    I_09 = imread([folder,'MP09_22C_y_hb_',emb,'.tif']);
    figure(2); clf; imshow(I_09);
    pause(1); 
    
    
     load([folder,'/',emb_roots{z},emb,'_data.mat']); 
    
end

%%
I = imread([folder,'MP09_22C_y_hb_49.tif']);
figure(1); clf; imshow(I);
