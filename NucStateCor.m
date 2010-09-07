
% Correlation Analysis


% Method
% Find 5-8 neighbors out from chosen cell. 
% clasify each neighbor 


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
     
     E=90;
NucCorrData = cell(E,4);
NCDevData = cell(E,4);  
     
for z=1:4;

    
    found = 0; n = 0;  
    while found ==0;  
    n=n+1;
        if n<10
            emb = ['0',num2str(n)];
        else
            emb = num2str(n);
        end

        try
        load([folder,'/',emb_roots{z},emb,'_data.mat']); 
        catch
        end
        
        found = 1; 
    end
  
  figure(1); clf; imshow(handles.It); title(emb_roots{z});
  pause(1); 

  
  clear C; 
  C= sparse(handles.conn_map);
%%


nuc_on1 = unique(pts1); 
N_tot = length(C);
N_on = length(ptr_nucin1);

Nstatus = zeros(N_tot,5);
Ndev = zeros(N_tot,5); 
 
on1 = zeros(1,N_on);

for j=1:N_on
    i=ptr_nucin1(j);

     N1 =  find(C(i,:)>4);
        N1in = intersect(ptr_nucin1,N1); % only neighbors also inside the region count 
        on1(j) = sum(ismember(N1in,nuc_on1))/length(N1); % fraction of layer 1 nuclei on in channel 1 
        on2 = zeros(1,length(N1)); 
     for k=1:length(N1) % loop through all neighbors 
         N2 = find(C(N1(k),:)>4);  % find all of these guys connections
          N2in = intersect(ptr_nucin1,N2); % only count the ones also inside the boundary
          on2(k) = sum(ismember(nuc_on1,N2in))/length(N2); % fraction of layer 2 nuclei on in channel 1 
          on3 = zeros(1,length(N2)); % initialize a new entry for each of these guys to store its neighbor stats.
         for k2=1:length(N2)
             N3 = find(C(N2(k2),:)>4);
             N3in = intersect(ptr_nucin1,N3);
             on3(k2) = sum(ismember(nuc_on1,N3in))/length(N3); % fraction of layer 3 nuclei on in channel 1  
             on4 = zeros(1,length(N3));
             for k3=1:length(N3)
                 N4 = find(C(N3(k3),:)>4);
                 N4in = intersect(ptr_nucin1,N4);
                 on4(k3) = sum(ismember(nuc_on1,N4in))/length(N4); % fraction of layer 4 nuclei on in channel 1
                 on5 = zeros(1,length(N4));
                 for k4=1:length(N4)
                     N5 = find(C(N4(k4),:)>4);
                     N5in = intersect(ptr_nucin1,N5);
                     on5(k4) = sum(ismember(nuc_on1,N5in))/length(N5); % fraction of layer 5 nuclei on in channel 1  
                 end
            end
         end
     end
   
     Nstatus(i,:) = [on1(j),nanmean(on2),nanmean(on3),nanmean(on4),nanmean(on5)];
     Ndev(i,:) = [0,nanstd(on2),nanstd(on3),nanstd(on4),nanstd(on5)];
     if rem(j,10)==0;
     disp(['Embryo ',emb, '  dataset ',num2str(z),  ' Progress: ',num2str(j/N_tot,2)]);  
     end
end
  
 length(nuc_on1)/length(unique(ptr_nucin1))
  figure(1);
   clf;
   subplot(3,1,1); hist(Nstatus(:,1),100);
   subplot(3,1,2); hist(Nstatus(:,2),100);
   subplot(3,1,3); hist(Nstatus(:,3),100);
   
   off_nucs = setdiff(1:length(C),nuc_on1);
   On_corr = nanmean(Nstatus(nuc_on1,:))
   figure(2); clf; plot(On_corr);
   
   Off_corr = nanmean(Nstatus(off_nucs,:)) 
   figure(3); clf; plot(Off_corr);
  
   
NucCorrData{z} = Nstatus; 
NCDevData{z} = Ndev;
end 
%    
%    Nstatus = zeros(N_tot,5); 
% for j=1:N_tot; % j=pts1(100)
%      N1{j} =  find(C(j,:)>4);
%         on1 = sum(ismember([N1{j}],nuc_on1))/length(N1{j}); % fraction of layer 2 nuclei on in channel 1 
%      for k=1:length(N1{j})
%          N2{k} = find(C(N1{j}(k),:)>4);
%           on2 = sum(ismember(nuc_on1,N2{k}))/length(N2{k}); % fraction of layer 2 nuclei on in channel 1 
%          for k2=1:length(N2{k})
%              N3{k2} = find(C(N2{k}(k2),:)>4);
%              on3 = sum(ismember(nuc_on1,N3{k2}))/length(N3{k2}); % fraction of layer 3 nuclei on in channel 1  
%              for k3=1:length(N3{k})
%                  N4{k3} = find(C(N3{k}(k3),:)>4);
%                  on4 = sum(ismember(nuc_on1,N4{k3}))/length(N4{k3}); % fraction of layer 3 nuclei on in channel 1  
%                  for k4=1:length(N4{k})
%                      N5{k4} = find(C(N4{k}(k4),:)>4);
%                      on5 = sum(ismember(nuc_on1,N5{k4}))/length(N5{k4}); % fraction of layer 3 nuclei on in channel 1  
%                  end
%             end
%          end
%      end
%      Nstatus(j,:) = [on1,on2,on3,on4,on5];
%      disp(['Embryo ',emb, '  dataset ',num2str(z),  ' Progress: ',num2str(j/N_tot,2)]);  
% end
   