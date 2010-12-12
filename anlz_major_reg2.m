% defines region using the im_nucdots_exon.m exported regions.  Chooses
% only major region from these.  (region based on disjoint nuclei looses
% too many detached nuclei from major region).  


function [foff, ptr_inhb] = anlz_major_reg2(folder,emb_roots,emb)


load([folder,'/',emb_roots,emb,'_data.mat']);      
 % load([folder,'/',emb_roots{z},'02','_data.mat']);
 
  Reg_size2 = zeros(1,handles.Rs2);
  for k=1:handles.Rs2
Reg_size2(k) = sum(diff(bndrys2{k}(:,1)).^2 +diff(bndrys2{k}(:,2)).^2);
  end
  
  [jnk,k] = max(Reg_size2);
       nin2 = logical( inpolygon(cent(:,1),cent(:,2),bndrys2{k}(:,1),bndrys2{k}(:,2)));
      
maj_reg = all_nucs(nin2); % pointers of nuclei in boundary2;     
ptr_inhb = intersect(ptr_nucin2,maj_reg); % pointer of nuclei in endogenous boundaries AND in major region of endogenous expression


hb_dots2 = intersect(pts2,ptr_inhb); % endogenous channel 
hb_dots1 = intersect(pts1,ptr_inhb); % reporter channel
noff = setdiff( hb_dots2,hb_dots1); % number not on in reporter, on in endogenous
foff = length(noff) ./ length(intersect(pts1,ptr_inhb)); % fraction not on it reporter

% % troubleshooting / validation code 
% T = ismember(H,maj_reg);
% 
% figure(1); clf; 
% imshow(T); pause(.001);

