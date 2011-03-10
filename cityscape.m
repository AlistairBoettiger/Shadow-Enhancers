%%                              cityscape.m

% Alistair Boettiger                                 Date Begun: 02/28/2011
% Levine Lab                                   Functional Since: 02/-/2010
%                                                 Last Modified: 02/28/2011

%% Description:
% function for taking cell data and making normalized cityscape plots of
% the data. These are simply just the tops of normalized histograms.  
%
%
% Called by Plot_DotComps family of analysis code.  
% Originally written to look at distribution of shadow data fraction of
% active transcription sites and variability in total transcript.  

%% Requires functions
% hist2dist.m 

%% Update Log:
% developed fro CompDist.m, 

%% Source Code:


function cityscape(data,names,xlab,F)
% 
% x = linspace(0,1,15);
%   xx = linspace(0,1,100); 
%  method = 'pcubic';
% sigma = .1; 
% 

% data = plot_miss;
G = length(data);
X = cell(1,G); 
Y = cell(1,G);
for k=1:G
    Es = length(nonzeros(data{k}));
    X{k} = linspace(0,1,10*log2(Es));
end


var_hist = cellfun(@nonzeros,data,'UniformOutput',false);
var_hist = cellfun(@hist,var_hist,X,'UniformOutput',false);
N_var = cellfun(@sum,var_hist,'UniformOutput',false);




% figure(1); clf; 
C = hsv(G);
legend_labels = cell(1,G);
for k=1:G  
     xx = imresize(X{k},2,'nearest');
     XX{k} = [xx(1,2:end), 1];
     y = imresize(var_hist{k},2,'nearest'); 
     Y{k} = y(1,:);  
    legend_labels{k} = [names{k}, ' N=', num2str(N_var{k} )]; 
    plot(XX{k},Y{k},'color',C(k,:),'LineWidth',3); hold on;
end
legend(legend_labels);
xlabel(xlab,'FontSize',F+1); 
set(gca,'FontSize',F);