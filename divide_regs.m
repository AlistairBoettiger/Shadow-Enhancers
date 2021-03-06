
%%                         Function divide_regs.m                        %%
%
%
%
% Alistair Boettiger                                   Date Begun: 11/16/10
% Levine Lab                                     Functional Since: 11/17/10
%                                                   Last Modified: 03/21/11

%% Description
% 11/17/10
% Divide major expression region into 3 AP sections, oriented from anterior
% to posterior (1:3).  
% 
% Rewritten 03/21/11 to divide expression based on rows of cells from
% posterior border of the anterior pattern.  
%% Updates



function layer_miss_rate = divide_regs(Reg1,Reg2,L2n1,H,cellbords,dispim)
%[f1,f2,f3] = divide_regs(L2,H,pts1,pts2,ptr_nucin2,In,dispim)

fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Shadow Enhancers/Code_Data/';

%%
% Get major region
    % Compute region sizes
    
% dispim = 1;
   
  fname = [folder,'/',emb_roots{z},emb,'_data.mat'];
            load(fname);   

    [hbreg,N] = bwlabel(Reg2);  
    
    
    
    [hbreg_props] = regionprops(hbreg);
    [area,major_reg] = max([hbreg_props.Area]);  % find major region
    disp(['Regions: ', num2str(N),'  Area =', num2str(area)]);

    % still need to keep small off-shoots near major region.  
    
    % create a mask that has only the major region
    [h,w] = size(H); 
    reg_mask = zeros([h,w]); 
    reg_mask(hbreg == major_reg) = 1;
    
    % figure(3); clf; imshow(reg_mask);

% Get indices of nuclei in anterior region
    hb_reg_labels = H.*reg_mask;  % indices of all nuclei in region
    nmin = min(hb_reg_labels(hb_reg_labels>0)); % nucleus with smallest index
    nmax = max(hb_reg_labels(hb_reg_labels>0)); % nucleus with largest index
    
    nin = nmax-nmin + 1; % number of nuclei in region 
    thrd = round(nin/3); % how many nucs in one third of region

% Divide region into thirds 
        n1s = nmin:nmin+thrd;       % indices of cells in region
        n2s =  nmin+thrd+1:nmin+2*thrd;
        n3s = nmin+2*thrd+1:nmax;
 

%     % Troubleshooting: test visualization
%         Ns1 = ismember(H, n1s); % map of cells in region
%         Ns2 = ismember(H, n2s);
%         Ns3 = ismember(H, n3s);
%         Ns = Ns1 + 2*Ns2 + 3*Ns3;
%         figure(1); clf; imagesc(Ns); 

% Figure out which is anterior and which is posterior
% 1 is for anterior posterior
Ntot = max(H(:));

 % relabel embryo
         % Create new bw map (old one was not saved)
            B = logical(H);   % needs to be on fullsize H because cellbords are on fullsize H; 
            B(cellbords) = 0;      
           %  figure(1); clf; imagesc(B); colormap hot;
    
    if  nmax < 3*Ntot/4;   % nmax < 7*Ntot/9; %  nmax < 5*Ntot/9; % 
        orientation = 0;               
        % Flip embryo
            B = fliplr(B); 
            Reg1 = fliplr(Reg1);
            Reg2 = fliplr(Reg2); 
            L2n1 = fliplr(L2n1);
            Cell_bnd = fliplr(Cell_bnd); 
            
        % Relabel nuclei
            NucLabeled = bwlabel(B,4);
        
    else
        orientation = 1;
        NucLabeled = H.*B; % keep current labels

    end
   
       Reg = ismember(NucLabeled,Ntot-nin+1:Ntot);
    %   figure(1); clf; imagesc(Reg); colormap hot;

       RegLabeled = bwlabel(Reg,4);
       R = regionprops(RegLabeled,'Centroid'); 
       reg_cents = reshape([R.Centroid],2,length(R)); 
       r_inds = sub2ind([h,w],floor(reg_cents(2,:)),floor(reg_cents(1,:)));
       d = reg_cents(1,:);
       



hbL =  floor(sqrt(nin)); % length of hb pattern in cells 


    % Sort by distance from upper left corner (min bcd).  
    nuc_order = NucLabeled(r_inds);
    [b,m,n] = unique(nuc_order);
    dists = d(m);
    dists = dists - min(dists); 
    dists = round(dists/max(dists)*hbL); % divide into 30
    
    if dispim == 1
        figure(1); clf; imagesc(L2n1); colormap hot;
    end
        
    RegD = RegLabeled;
    ns = min(nin,length(dists));
        for n=1:ns
            if orientation == 1; 
                RegD(RegLabeled==n) = dists(n);
            else
                 RegD(RegLabeled==n) = dists(ns-n+1);
            end
        end
        
        % figure(1); clf; imagesc(RegD); colormap(lines);colorbar;
        
        C = [1,1,1; colormap(jet(15)); 0,0,0]; 
        
        % C = [1,1,1; colormap(copper(15)); 0,0,0]; 
        MissedD = RegD.*L2n1.*Reg1.*Reg2; % must be inside both regs 
     if dispim == 1
        figure(2); clf; imagesc(flipud(fliplr(MissedD + 5*Cell_bnd )));  colormap(C);%  lines;
        axis off; colorbar;
     end
        
        layer_miss_cnt = hist(nonzeros(MissedD(r_inds)),1:hbL);
        layer_length = hist(nonzeros(RegD(r_inds)),1:hbL);
        layer_miss_rate = layer_miss_cnt./layer_length;
        
    if dispim == 1
        figure(3); clf; bar(layer_miss_rate);
    end
    
   
%     
%     
% % Assign Anterior as Sect1, and boundary region as Sect3
%     if orientation == 1  
% %   Map only for trouble shooting
% %         Sect1 = Ns1; % map of cells in region, AP oriented
% %         Sect2 = Ns2;
% %         Sect3 = Ns3;
%         s1s = n1s; % indices of cells in region, AP oriented
%         s2s = n2s;
%         s3s = n3s;
%     else
% %         Sect1 = Ns3; % map of cells in region, AP oriented
% %         Sect2 = Ns2;
% %         Sect3 = Ns1;
%         s1s = n3s; % indices of cells in region, AP oriented
%         s2s = n2s;
%         s3s = n1s; 
%     end
%     onReg1 = intersect(s1s,setdiff(pts2,pts1));
%     onReg2 = intersect(s2s,setdiff(pts2,pts1));
%     onReg3 = intersect(s3s,setdiff(pts2,pts1));
%     
%     
% %     figure(2); clf; subplot(3,1,1); imshow(Ired);
% %     subplot(3,1,2); imshow(Igreen);
% %     subplot(3,1,3); imshow(Idif);
%  
% if dispim == 1
%     figure(1); clf; % In = handles.In;
%      imReg1 = ismember(H,onReg1);
%      imReg2 = ismember(H,onReg2);
%      imReg3 = ismember(H,onReg3);
%     imReg = double(In)/255+(imReg1 + 2*imReg2 + 3*imReg3);
%     imagesc(imReg); colormap(jet);
%     pause(.1);
% end 
% 
% % compute fraction of INactive cells in region
% f1 = length(onReg1)/length(s1s); % anterior most region
% f2 = length(onReg2)/length(s2s); % middle region
% f3 = length(onReg3)/length(s3s); % boundary
% 
% save([fout,'test2']);

