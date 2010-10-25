function [mc, ptr_nucin2] = anlz_major_reg(folder,emb_roots,emb)

    load([folder,'/',emb_roots,emb,'_data.mat']);   
    
Reg_size1 = zeros(1,handles.Rs1);
  for k=1:handles.Rs1
Reg_size1(k) = sum(diff(bndrys1{k}(:,1)).^2 +diff(bndrys1{k}(:,2)).^2);% perimiter of regions for mRNA1  
  end

  [jnk,k1] = max(Reg_size1);
  
  Reg_size2 = zeros(1,handles.Rs2);
  for k=1:handles.Rs2
Reg_size2(k) = sum(diff(bndrys2{k}(:,1)).^2 +diff(bndrys2{k}(:,2)).^2)
  end
  
  [jnk,k2] = max(Reg_size2);
  
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~mRNA 1 Stuff~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
     % Nuclei and RNA in each region 
      nin1 = false; RNA_in1 = false;
    for k=k1
       nin1 = logical( inpolygon(cent(:,1),cent(:,2),bndrys1{k}(:,1),bndrys1{k}(:,2)) + nin1 ) ;  
       RNA_in1 = logical( inpolygon(centRNA1(:,1),centRNA1(:,2),bndrys1{k}(:,1),bndrys1{k}(:,2)) + RNA_in1 ); 
    end
        rna_x1 = centRNA1(RNA_in1(:,1),1);       
        rna_y1 = centRNA1(RNA_in1(:,1),2);
%        rna_x1 = centRNA1(:,1); 
%        rna_y1 = centRNA1(:,2); 
    
    % get pointers (indices) of the nuclei in each region
         all_nucs = 1:max(max(H));
         ptr_nucin1 = all_nucs(nin1);
  % Identify all pixels of the nuclei which contain a nasenct transcript
    inds1 = round(rna_y1)+round(rna_x1)*h;   % convert x-y indexing to linear indexing
    pts1 = nonzeros(H(inds1)); % Indices of nuclei who contain an mRNA_1 transcript
    L1a = ismember(H,pts1); 
    Reg1 = ismember(H,ptr_nucin1); 
     
     
         [mRNAs_per_nucF1,mRNAs_per_nuc1]  = mRNA_border(H,cellbords,bnd_buff,inds1,pars); 
          def_on1 = all_nucs(mRNAs_per_nuc1 > 0);
          ons1 = all_nucs(mRNAs_per_nucF1 > 0);
            L1 = ismember(H,ons1);

 %~~~~~~~~~~~~~~~~~~~~~~~~~~end mRNA 1 stuff~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
     
     
 %~~~~~~~~~~~~~~~~~~~~~~~~~~~mRNA 2 Stuff~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%   
      nin2 = false; 
    for k=k2
       nin2 = logical( inpolygon(cent(:,1),cent(:,2),bndrys2{k}(:,1),bndrys2{k}(:,2)) + nin2);
    end        
 ptr_nucin2 = all_nucs(nin2); % pointers of nuclei in boundary2;     
 Reg2 = ismember(H,ptr_nucin2);   
 
 % if region 2 is just exonic transcript, there are no nascent transcripts
 % to process.  We just use the expression region.  
    %-----------------------------------------%
     if isempty(centRNA2)==0   
            RNA_in2 = false; 
            for k=1:handles.Rs2
               RNA_in2 = logical(inpolygon(centRNA2(:,1),centRNA2(:,2),bndrys2{k}(:,1),bndrys2{k}(:,2)) + RNA_in2);
            end            
               rna_x2 = centRNA2(RNA_in2(:,1),1); 
               rna_y2 = centRNA2(RNA_in2(:,1),2);            
        %        rna_x2 = centRNA2(:,1); 
        %        rna_y2 = centRNA2(:,2); 

            inds2 = round(rna_y2)+round(rna_x2)*h;  % convert x-y indexing to linear indexing    
            pts2 = nonzeros(H(inds2));% Indices of nuclei who contain an mRNA_2 transcript
            L2a = ismember(H,pts2);    
             
          [mRNAs_per_nucF2,mRNAs_per_nuc2]  = mRNA_border(H,cellbords,bnd_buff,inds2,pars); 
          ons2 = all_nucs(mRNAs_per_nucF2 > 0);
          def_on2 = all_nucs(mRNAs_per_nuc2 > 0);
          L2 = ismember(H,ons2); 
            
     else
         pts2 = ptr_nucin2;
         L2a = Reg2;
         L2 = L2a; 
         Chns = 1; 
         rna_x2 = NaN;
         rna_y2 = NaN;
     end   
 %~~~~~~~~~~~~~~~~~~~~~~~~~~end mRNA 2 stuff~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~% 
 
%length(unique(ons1)) - length(unique(pts1))
Cell_bnd = false(h,w);
Cell_bnd(cellbords) = 1;
Cell_bnd = imdilate(Cell_bnd,strel('disk',1)); 
Nuc_line = H; Nuc_line(Nuc_line>0)=1; Nuc_line(Cell_bnd == 1)=0;
Nuc_cntr = H; Nuc_cntr(Nuc_cntr>0)=1; Nuc_cntr(imdilate(Cell_bnd,strel('disk',1))==1)=0;
Cell_bnd = logical(Nuc_line-Nuc_cntr); 



% % show effect of disambiguation
%        L2n1 = (L2-L1) > 0; 
%       L2n1a = (L2a - L1a) > 0;
%  figure(1); clf;
%     Io = uint8(zeros(h,w,3));
%     Io(:,:,1) = imadd(uint8(0*L2n1a.*Cell_bnd),1*handles.Im1); 
%     Io(:,:,2) =  imadd(uint8(255*L2n1.*Cell_bnd),1.2*handles.Im2); % imadd(uint8(255*Cell_bnd),1*handles.Im2);  %
%     Io(:,:,3) =   imadd(uint8(255*L2n1a.*Cell_bnd),.8*handles.In); % uint8(255*Cell_bnd); %
%     % DI = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
%   clf; imshow(Io); hold on;
% plot(0,0,'b','LineWidth',3); plot(0,0,'g','LineWidth',3); plot(0,0,'c','LineWidth',3);
% legend('was off','newly off','still off');
%  %  plot(rna_x1,rna_y1,'r.','MarkerSize',5);
% % plot(rna_x2,rna_y2,'g.','MarkerSize',5);
%      title('disambiguation of mRNA 1'); 
     
  
% % Outline missing nuclei  load temp_dat;
%        L2n1 = (L2-L1) > 0; 
%       L2n1a = (L2a - L1a) > 0;
%  figure(1); clf;
%     Io = uint8(zeros(h,w,3));
%     Io(:,:,1) = imadd(uint8(0*L2n1a.*Cell_bnd),1*handles.Im1); 
%     Io(:,:,2) =  imadd(uint8(255*Reg1.*L2n1a.*Cell_bnd),1*handles.Im2); % imadd(uint8(255*Cell_bnd),1*handles.Im2);  %
%     Io(:,:,3) =   imadd(uint8(255*L2n1a.*Cell_bnd),.8*handles.In); % uint8(255*Cell_bnd); %
%     % DI = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
%   clf; imshow(Io); hold on;
% plot(0,0,'b','LineWidth',3); plot(0,0,'c','LineWidth',3); 
% legend('off, not in region','off, in region','still off');
%  %  plot(rna_x1,rna_y1,'r.','MarkerSize',5);
% % plot(rna_x2,rna_y2,'g.','MarkerSize',5);
%      title('transcription activity difference'); 


% %   
%   figure(2); clf; 
%     Io = uint8(zeros(h,w,3));
%     Io(:,:,1) = uint8(255*L2n1.*Reg1);
%      Io(:,:,2) = uint8(255*L2n1.*Reg1);
%      Io(:,:,3) = 2*handles.In;
%      Idif = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
%      imshow(Idif); hold on;
% 
%  %    
%     figure(6); clf; 
%     set(gcf,'color','k'); 
%  subplot(2,1,1); 
%      Io = uint8(zeros(h,w,3));
%      Io(:,:,1) = uint8(255*L1);
%      Io(:,:,3) = 30*handles.In;
%      Ired = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
%      imshow(Ired); hold on;
%       %plot(rna_x1,rna_y1,'r.');  
%  subplot(2,1,2); 
%      Io = uint8(zeros(h,w,3));
%      Io(:,:,2) = uint8(255*L2);
%      Io(:,:,3) = 30*handles.In;
%      Igreen = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
%      imshow(Igreen); hold on;
% % plot(rna_x2,rna_y2,'g.');
%  
%  
%  figure(2); clf;
%      Io = uint8(zeros(h,w,3));
%     Io(:,:,1) = uint8(255*L1);
%     Io(:,:,2) = uint8(255*L2);
%     Io(:,:,3) = 30*handles.In;
%      DI = uint8(bsxfun(@times,double(Io)/255,double(handles.In)));
%   clf; imshow(DI);
%     for k=1:handles.Rs1
%  hold on; plot(bndrys1{k}(:,1),bndrys1{k}(:,2),'r','LineWidth',3);
%     end
%     for k=1:handles.Rs2
%  hold on; plot(bndrys2{k}(:,1),bndrys2{k}(:,2),'g','LineWidth',3);
%     end
%  hold on;
%  plot(rna_x1,rna_y1,'r.');
%  plot(rna_x2,rna_y2,'g.');
%   

mc = length(intersect(setdiff(pts2,pts1), ptr_nucin2));  

