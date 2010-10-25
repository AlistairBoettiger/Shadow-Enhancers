
plot_miss = cell(1,G); 
 for k=1:5:10;     plot_miss{k} = miss_rate{k}(cc13{k}); end

 plot_miss{2} = [plot_miss{1};plot_miss{6}];

 
  [mc,mr,nd2,er] = merge_data(1,6,N,miss_cnt,miss_rate,nd,ectop_rate);
  
  ND = cell2mat(nd2); age_offset = 4.8;
emb_cycle = age_offset + log2( nonzeros( sort(ND(:)) ) );
cc14 =cell(1,G); cc13 = cell(1,G); cc12 = cell(1,G); cc11 = cell(1,G); cc10 = cell(1,G); cc9 = cell(1,G); 
for z=1:G-1
    logage =   age_offset + log2( ND(:,z) );
    cc14{z} = logage >14;
    cc13{z} = logage <14  & logage> 13;
    cc12{z}  = logage <13 & logage > 12;
    cc11{z} = logage <12 ;
end

  
  
  for k=1;     plot_miss{k+2} = mr{k}(cc13{k}); end
 
figure(30); clf;

x = linspace(0,1,20);  % range and number of bins for histogram
xx = linspace(0,1,pts); % range a number of bins for interpolated distribution
 method = 'pcubic'; % method for interpolation
sigma = .1;  % smoothing factor for interpolation
CompDist(plot_miss,x,xx,method,sigma,names,xlab,F)
ylim([0,ymax]); title('cc13');


