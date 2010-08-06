%%                              CompDist.m

% Alistair Boettiger                                 Date Begun: 04/04/2010
% Levine Lab                                   Functional Since: 04/04/2010
%                                                 Last Modified: 04/04/2010

%% Description:
% function for taking cell data and making normalized frequency
% distributions of the data.  
% Called by Plot_DotComps family of analysis code.  
% Originally written to look at distribution of shadow data fraction of
% active transcription sites and variability in total transcript.  

%% Requires functions
% hist2dist.m 

%% Update Log:

%% Source Code:


function CompDist(data,x,xx,method,sigma,names,xlab)
% 
%  x = linspace(0,1,12);
%   xx = linspace(0,1,100); 
%  method = 'pcubic';
% sigma = .1; 
% 
G = length(data);
X = cell(1,G); XX = cell(1,G); M = cell(1,G); S = cell(1,G);


for k=1:G
    X{k} = x;
    XX{k} =xx;
    M{k} = method;
    S{k} = sigma;
end


var_hist = cellfun(@nonzeros,data,'UniformOutput',false);
var_hist = cellfun(@hist,var_hist,X,'UniformOutput',false);
N_var = cellfun(@sum,var_hist,'UniformOutput',false);
dist_var = cellfun(@hist2dist,var_hist,X,XX,M,S,'UniformOutput',false);


C = hsv(G);
legend_labels = cell(1,G);
for k=1:G
    legend_labels{k} = [names{k}, ' N=', num2str(N_var{k} )]; 
    plot(xx,dist_var{k},'color',C(k,:),'LineWidth',3); hold on;
end
legend(legend_labels);
xlabel(xlab,'FontSize',16); 
set(gca,'FontSize',15);