%%                          merge_data.m 
%
% Alistair Boettiger                                   Date Begun: 06/28/10
% Levine Lab                                        Last Modified: 06/28/10

%% Description
% Written to merge dot data from different tracks, a and b, into data for
% track a.  Combined nonzeros lengths must be less than max length of data
% in analysis.  

%% Called from:
% Plot_DotComs2_sna.m

% has bugs
% disabled function call until bugs can be fixed.

function disabled 10/21/10 [miss_cnt,miss_rate,nd,lowon] = merge_data(a,b,N,miss_cnt,miss_rate,nd,lowon)

% merge data seemlessly
d_cnt = [nonzeros(miss_cnt{a}); nonzeros(miss_cnt{b})];
d_rate = [nonzeros(miss_rate{a}); nonzeros(miss_rate{b})];
d_nd =  [nonzeros(nd{a}); nonzeros(nd{b})];
d_lowon = [nonzeros(lowon{a}); nonzeros(lowon{b})];

% update data structures, respecting size (for later manipulations).
miss_cnt{a} = [d_cnt; zeros(N-length(d_cnt),1) ] ; 
miss_rate{a} = [d_rate; zeros(N-length(d_rate),1) ] ; 
nd{a} =[d_nd; zeros(N-length(d_nd),1) ] ; 
lowon{a} = [d_lowon; zeros(N-length(d_lowon),1) ] ; 

K = length(miss_cnt); 

% Remove track b
miss_cnt(b:K) = {miss_cnt{b+1:K},[]};
miss_rate(b:K) = {miss_rate{b+1:K},[]};
nd(b:K) = {nd{b+1:K},[]};
lowon(b:K) = {lowon{b+1:K},[]};


% % old command
% miss_cnt{1} = [miss_cnt{1}(1:N/2); miss_cnt{5}(1:N/2)]; 
% miss_rate{1} = [miss_rate{1}(1:N/2); miss_rate{5}(1:N/2)]; 
% nd{1} = [nd{1}(1:N/2); nd{5}(1:N/2)]; 
% lowon{1} = [lowon{1}(1:N/2); lowon{5}(1:N/2)]; 