%%                              BoxDist.m

% Alistair Boettiger                                 Date Begun: 04/04/2010
% Levine Lab                                   Functional Since: 04/04/2010
%                                                 Last Modified: 04/04/2010

%% Description:


%% Requires functions


%% Update Log:

%% Source Code:


function BoxDist(data,names,xlab)
% 
% 
% data = plot_lowon;
% data = plot_miss; 

G = length(data); 
% C1 = hsv(G); C = 1:G;
C1 = [1,0,.4;
     0,.4,1];
  
 C = zeros(1,G); 
for k=1:G/2;     
C(1+(k-1)*2:k*2) = [1,2];
end

B = zeros(100,G); 
for k=1:G
 B(1:length(data{k}),k) = data{k};
end
Ns = sum(B>0);
B(B==0) = NaN; 

legend_labels = cell(1,G); 
for k=1:G
    legend_labels{k} = [names{k}, ' N=', num2str(Ns(k) )]; 
end

B = fliplr(B);
legend_labels = fliplr(legend_labels); 

boxplot(B,'orientation','horizontal','whisker',1,'colorgroup',C,'colors',C1,'labels',legend_labels); % ,'medianstyle','target'
hold on;
plot(-1,-1,'color',C1(1,:)); 
plot(-1,-1,'color',C1(2,:));
legend('30C','22C'); 
set(gca,'FontSize',15);
xlabel(xlab,'FontSize', 16);


%legend(legend_labels,'Location','SouthOutside'); 
% set(gca,'XTick',[],'XTickLabel',{' '},'FontSize',15);

