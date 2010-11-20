
%%                         Function divide_regs.m                        %%
%
%
%
% Alistair Boettiger                                   Date Begun: 11/16/10
% Levine Lab                                     Functional Since: 11/17/10
%                                                   Last Modified: 11/17/10

%% Description
% Divide major expression region into 3 AP sections, oriented from anterior
% to posterior (1:3).  
%
%
%% Updates



function [f1,f2,f3] = divide_regs(L2,H,pts1,pts2,ptr_nucin2,In)
%%
% Get major region
    % Compute region sizes
    endog_reg = ismember(H,ptr_nucin2);
    
    [hbreg,N] = bwlabel(endog_reg); 

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

    if nmax > 3*Ntot/4;  
        orientation = 0;
    else
        orientation = 1;
    end
% Assign Anterior as Sect1, and boundary region as Sect3
    if orientation == 1  
%   Map only for trouble shooting
%         Sect1 = Ns1; % map of cells in region, AP oriented
%         Sect2 = Ns2;
%         Sect3 = Ns3;
        s1s = n1s; % indices of cells in region, AP oriented
        s2s = n2s;
        s3s = n3s;
    else
%         Sect1 = Ns3; % map of cells in region, AP oriented
%         Sect2 = Ns2;
%         Sect3 = Ns1;
        s1s = n3s; % indices of cells in region, AP oriented
        s2s = n2s;
        s3s = n1s; 
    end
    onReg1 = intersect(s1s,setdiff(pts2,pts1));
    onReg2 = intersect(s2s,setdiff(pts2,pts1));
    onReg3 = intersect(s3s,setdiff(pts2,pts1));
    
    
%     figure(2); clf; subplot(3,1,1); imshow(Ired);
%     subplot(3,1,2); imshow(Igreen);
%     subplot(3,1,3); imshow(Idif);
 
    figure(1); clf; % In = handles.In;
     imReg1 = ismember(H,onReg1);
     imReg2 = ismember(H,onReg2);
     imReg3 = ismember(H,onReg3);
    imReg = double(In)/255+(imReg1 + 2*imReg2 + 3*imReg3);
    imagesc(imReg); colormap(jet);
    pause(.1);
    
% compute fraction of INactive cells in region
f1 = length(onReg1)/length(s1s); % anterior most region
f2 = length(onReg2)/length(s2s); % middle region
f3 = length(onReg3)/length(s3s); % boundary

