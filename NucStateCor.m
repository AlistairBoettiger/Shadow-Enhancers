
%%                          Correlation Analysis


% Alistair Boettiger                                   Date Begun: 09/06/10
% Levine Lab                                This version complete: 09/08/10
%                                                   Last Modified: 09/08/10

% Method
% Find 5 neighbors out from chosen cell. 
% clasify each neighbor 
clear all;


folder =  '/Volumes/Data/Lab Data/Shadow_data/Processed';
emb_roots = {'MP09_22C_y_hb';
            'MP02_22C_y_hb';
            'MP02_30C_y_hb';
            'MP02_30C_LacZ_hb'   % yes it's actually yellow
            };  
          

names = {'2 enhancers, 22C';
         'no shadow, 22C';
         'no shadow, 30C'
         'no shadow, 30C'
         };
     
     N=90;
     Z = 4; % total number of data channels
NucCorrData = cell(N,Z);
NCDevData = cell(N,Z);  
     
for z=1:Z;

    
for n=1:N
        if n<10
            emb = ['0',num2str(n)];
        else
            emb = num2str(n);
        end

        try
        load([folder,'/',emb_roots{z},emb,'_data.mat']); 
   
%   figure(1); clf; imshow(handles.It); title(emb_roots{z});
%   pause(1); 

  C= sparse(handles.conn_map);
%%


nuc_on1 = unique(pts1); 
N_tot = length(C);
N_on = length(ptr_nucin1);

Nstatus = NaN*zeros(N_tot,5);
Ndev = NaN*zeros(N_tot,5); 
 
on1 = zeros(1,N_on);

for j=1:N_on
    i=ptr_nucin1(j);

     N1 =  find(C(i,:)>4);
        N1in = intersect(ptr_nucin1,N1); % only neighbors also inside the region count 
        on1(j) = sum(ismember(N1in,nuc_on1))/length(N1in); % fraction of layer 1 nuclei on in channel 1 
        on2 = zeros(1,length(N1)); 
     for k=1:length(N1) % loop through all neighbors 
         N2 = find(C(N1(k),:)>4);  % find all of these guys connections
          N2in = intersect(ptr_nucin1,N2); % only count the ones also inside the boundary
          on2(k) = sum(ismember(nuc_on1,N2in))/length(N2in); % fraction of layer 2 nuclei on in channel 1 
          on3 = zeros(1,length(N2)); % initialize a new entry for each of these guys to store its neighbor stats.
         for k2=1:length(N2)
             N3 = find(C(N2(k2),:)>4);
             N3in = intersect(ptr_nucin1,N3);
             on3(k2) = sum(ismember(nuc_on1,N3in))/length(N3in); % fraction of layer 3 nuclei on in channel 1  
             on4 = zeros(1,length(N3));
             for k3=1:length(N3)
                 N4 = find(C(N3(k3),:)>4);
                 N4in = intersect(ptr_nucin1,N4);
                 on4(k3) = sum(ismember(nuc_on1,N4in))/length(N4in); % fraction of layer 4 nuclei on in channel 1
                 on5 = zeros(1,length(N4));
                 for k4=1:length(N4)
                     N5 = find(C(N4(k4),:)>4);
                     N5in = intersect(ptr_nucin1,N5);
                     on5(k4) = sum(ismember(nuc_on1,N5in))/length(N5in); % fraction of layer 5 nuclei on in channel 1  
                 end
            end
         end
     end
   
     Nstatus(i,:) = [on1(j),nanmean(on2),nanmean(on3),nanmean(on4),nanmean(on5)];
     Ndev(i,:) = [0,nanstd(on2),nanstd(on3),nanstd(on4),nanstd(on5)];
     if rem(j,10)==0;
     disp(['Embryo ',emb, '  dataset ',num2str(z),  ' Progress: ',num2str(j/N_on,2)]);  
     end
end
  
  fon = 100*length(nuc_on1)/length((ptr_nucin1));
 % figure(1);

   
   off_nucs = setdiff(ptr_nucin1,nuc_on1);
   On_corr = 100*nanmean(Nstatus(nuc_on1,:))
   On_corr_err = 100*nanmean(Ndev(nuc_on1,:));
 %  figure(2); clf; errorbar(1:5,On_corr,On_corr_err);
   
   Off_corr = 100*nanmean(Nstatus(off_nucs,:))
   Off_corr_err = 100*nanmean(Ndev(off_nucs,:));
%   figure(3); clf; errorbar(1:5,Off_corr,Off_corr_err);
 
        catch
        end
   
   
NucCorrData{n,z} = [N_tot, fon, On_corr, On_corr_err, Off_corr, Off_corr_err]; 
end 
end
%%

save hbCorrData3; %hbCorrData2 good run

%%
N = 50;
z=4;
lab = cell(1,N);

figure(1); clf; figure(2); clf;
for j=1:N
    oncor = [NucCorrData{j,z}(2), NucCorrData{j,z}(3:7)] ;% /NucCorrData{j,z}(2); 
    offcor = [NucCorrData{j,z}(2), NucCorrData{j,z}(13:17)];% /NucCorrData{j,z}(2); 
    figure(1); hold on; plot(oncor,'color',[1-j/N,0,j/N],'Marker','.');
    figure(2); hold on; plot(offcor,'color',[1-j/N,0,j/N],'Marker','.');
     lab{j} = ['Nuclei: '  num2str(NucCorrData{j,z}(1))];
end
    legend(lab);
 
    z=1;
NCD = cell2mat(NucCorrData);

z=4;
cc9s = find(NCD(:,(z-1)*22+1)<500);
cc13s = find(NCD(:,(z-1)*22+1)>1200  & NCD(:,(z-1)*22+1)<1800);
cc14s = find(NCD(:,(z-1)*22+1)>1800);

embs_chosen =  cc14s;
mean(NCD(embs_chosen,(z-1)*22+2))
mean(NCD(embs_chosen,(z-1)*22+13)./NCD(embs_chosen,(z-1)*22+2))





%%  Simulate correlated activation


%~~~~~~~~~~~~~~  Create null model with random probability of being on
frac_on = .2; % length(nuc_on1)/length((ptr_nucin1));
uncorr = find(rand(1,N_tot)>frac_on);
Ons = ismember(H,uncorr);   

   figure(6); clf; 
     Io = uint8(zeros(h,w,3));
     Io(:,:,1) = uint8(255*Ons);
     Io(:,:,3) = 30*handles.In;
     Ion = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
     imshow(Ion); hold on;
      
      
frac_on = .5*length(nuc_on1)/length((ptr_nucin1));
rdraw = rand(1,N_tot);
onstate = rdraw>frac_on;
uncorr = find(onstate);

Ons = ismember(H,uncorr);         
   figure(6); clf; 
     Io = uint8(zeros(h,w,3));
     Io(:,:,1) = uint8(255*Ons);
     Io(:,:,3) = 30*handles.In;
     Ion = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
     imshow(Ion); hold on;
      
     for j=1:N_tot
         if onstate(j) == 0       
            Neibs = find(C(j,:)>4);
            draw2 = onstate(Neibs).*rand(1,length(Neibs))'; % everybody who's on gets a draw
            onstate(j)= draw2>1.75; % if the sum of the draws is greater than x, this off guy goes on
         end
     end
            
            
  corrstate = find(onstate);
  Ons = ismember(H,corrstate);  
     figure(5); clf; 
     Io = uint8(zeros(h,w,3));
     Io(:,:,1) = uint8(255*Ons);
     Io(:,:,3) = 30*handles.In;
     Ion = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
     imshow(Ion); 
  
%~~~~~~~~~~~~~~~
  
nuc_on1 =   corrstate; 
ptr_nucin1 = 1:N_tot;
  fon = 100*length(nuc_on1)/length((ptr_nucin1))

      for j=1:N_tot
    i=j;

     N1 =  find(C(i,:)>4);
        N1in = intersect(ptr_nucin1,N1); % only neighbors also inside the region count 
        on1(j) = sum(ismember(N1in,nuc_on1))/length(N1in); % fraction of layer 1 nuclei on in channel 1 
        on2 = zeros(1,length(N1)); 
     for k=1:length(N1) % loop through all neighbors 
         N2 = find(C(N1(k),:)>4);  % find all of these guys connections
          N2in = intersect(ptr_nucin1,N2); % only count the ones also inside the boundary
          on2(k) = sum(ismember(nuc_on1,N2in))/length(N2in); % fraction of layer 2 nuclei on in channel 1 
          on3 = zeros(1,length(N2)); % initialize a new entry for each of these guys to store its neighbor stats.
         for k2=1:length(N2)
             N3 = find(C(N2(k2),:)>4);
             N3in = intersect(ptr_nucin1,N3);
             on3(k2) = sum(ismember(nuc_on1,N3in))/length(N3in); % fraction of layer 3 nuclei on in channel 1  
             on4 = zeros(1,length(N3));
             for k3=1:length(N3)
                 N4 = find(C(N3(k3),:)>4);
                 N4in = intersect(ptr_nucin1,N4);
                 on4(k3) = sum(ismember(nuc_on1,N4in))/length(N4in); % fraction of layer 4 nuclei on in channel 1
                 on5 = zeros(1,length(N4));
                 for k4=1:length(N4)
                     N5 = find(C(N4(k4),:)>4);
                     N5in = intersect(ptr_nucin1,N5);
                     on5(k4) = sum(ismember(nuc_on1,N5in))/length(N5in); % fraction of layer 5 nuclei on in channel 1  
                 end
            end
         end
     end
   
     Nstatus(i,:) = [on1(j),nanmean(on2),nanmean(on3),nanmean(on4),nanmean(on5)];
     Ndev(i,:) = [0,nanstd(on2),nanstd(on3),nanstd(on4),nanstd(on5)];
     if rem(j,10)==0;
     disp(['Embryo ',emb, '  dataset ',num2str(z),  ' Progress: ',num2str(j/N_on,2)]);  
     end
end
  
  fon = 100*length(nuc_on1)/length((ptr_nucin1));
 % figure(1);

   
   off_nucs = setdiff(ptr_nucin1,nuc_on1);
   On_corr = 100*nanmean(Nstatus(nuc_on1,:))
   On_corr_err = 100*nanmean(Ndev(nuc_on1,:));
 %  figure(2); clf; errorbar(1:5,On_corr,On_corr_err);
   
   Off_corr = 100*nanmean(Nstatus(off_nucs,:))
   Off_corr_err = 100*nanmean(Ndev(off_nucs,:));
%   figure(3); clf; errorbar(1:5,Off_corr,Off_corr_err);

