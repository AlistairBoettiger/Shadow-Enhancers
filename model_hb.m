

function model_hb
clear all;

N  = 100;  T = 1500;
x = linspace(0,1,N);

tau = 1/5; % from Gregor 2007
bcd = exp(-x/tau);

% figure(1); clf; plot(x,bcd,'g');


hbM = .1;



hb = hbM*ones(N,1); 
P = ones(N,1); 
PB = zeros(N,1);
PBH = zeros(N,1); 
PBH2 = zeros(N,1); 

R0 = [hb;P;PB;PBH;PBH2];
[t,R] = ode15s(@sys,[0,T],R0,[],[N,bcd]);
hb = R(:,1:N); 
P=R(:,N+1:2*N);  % 
PB=R(:,2*N+1:3*N);  % 
PBH=R(:,3*N+1:4*N);  % 
PBH2=R(:,4*N+1:5*N); 
save test;

figure(1); set(gcf,'color','w');
for i=1:length(t)
    figure(1); clf; subplot(2,1,1); 
    plot(x,hb(i,:),'b'); hold on; 
    plot(x,bcd,'g');
   
   figure(1); subplot(2,1,2);
   plot(x,P(i,:),'k'); hold on;
   plot(x,PB(i,:),'c'); 
    plot(x,PBH(i,:),'m');
    plot(x,PBH2(i,:),'r');
end


function dRdt = sys(t,R,Pars)
N = Pars(1);
bcd = Pars(2:end)';


kd = 10;
k_bon = 2;
k_boff = .5;
k_hon = 1;
k_hoff = .5;
k_hon2 = 6;
k_hoff2 = .1; 
k_bs = 1;
k_hbs = 5;

hb = R(1:N); %
P=R(N+1:2*N);  % 
PB=R(2*N+1:3*N);  % 
PBH=R(3*N+1:4*N);  % 
PBH2=R(4*N+1:5*N); 

% promoter binding
dP = - k_bon.*bcd.*P + k_boff.*PB;
dPB = k_bon.*bcd.*P - k_boff.*PB - k_hon.*PB.*hb + k_hoff.*PBH;
dPBH = k_hon.*PB.*hb - k_hoff.*PBH - k_hon2.*PBH.*hb + k_hoff2.*PBH2;
dPBH2 =  k_hon2.*PBH.*hb - k_hoff2.*PBH2;

% hb synthesis and decay
dhb = k_bs*PB + k_hbs*PBH2 - kd*hb - k_hon.*PB.*hb - k_hon2.*PBH.*hb;
dRdt = [dhb; dP; dPB; dPBH; dPBH2];



%%