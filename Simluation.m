%Simularion of Networks Control project
clear all
close all
%%% Initizliation 
%% changing parameters 

% Set of optimal gain for different alpha1/alpha2 ratio:
AirCireff = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%hi=10;  
%ho=5; 

% 100
Ki = 0.008;
Kp = 9.8121;

% 50
Kp =7.0009;
Ki = 0.0053;

% 10
Kp =3.2298;
Ki = 0.0016;

% 7.5
Kp =   2.8369 ; 
Ki=  0.0011;

% 5
Kp = 2.4439;
Ki = 2.8620e-07;

Ki = Ki*10;

% 1
%Kp = 1.0189;
%Ki = 2.2394e-8;



[Ts,Bp,P,n,A,L,Desired_Temp,dis1,dis2,dis3,Temp_hvac]=Model();


% Siimulation time (seconds)
Tt=3600*30;



%% variables define

X = zeros(length(L),Tt);
U= zeros(1,Tt);

Disturbance = zeros(length(L),Tt);

%% initial conditions


C = zeros(Tt,7);


i_err = 0;

initail_temp=32.22;
initail_temp=21.11;
initail_temp=30;
 X(1:15,1) = initail_temp;
 X(1:15,2) = initail_temp; % Assign initial temperature to the rooms
%X(1:7,1) = [60,30,40,90,80,70,50]';
%X(1:7,2) = [60,30,40,90,80,70,50]'; % Assign initial temperature to the rooms


%initial location of observer
OLS = ceil(n*rand(1)); % select a random starting room for resident
                        % OLS position of resident
OLS = 4;
                        
res_temp= zeros(1,Tt);
res_temp(1,1:2)=[initail_temp,initail_temp];

time= zeros(1,Tt);
%%% system equations



for k=3:Tt

 %%%   Disturbance generator
   
%for solar radiation, we consider simus disturbance with density: with
%period 12 hours
radiation_ave=0;
radiation_peak=2;
%ds1= dis1*[2*(radiation_ave+radiation_peak*sin(k/43200)),2*(radiation_ave+radiation_peak*sin(k/43200)),2*(radiation_ave+radiation_peak*sin(k/43200+pi/2)),(radiation_ave+radiation_peak*sin(k/43200+pi/2)),3*(radiation_ave+radiation_peak*sin(k/43200+pi/2))];
ds1= dis1*2*(radiation_ave+radiation_peak*sin(k/43200))* [2,2,2,1,3];

%Number of wall are not included, however, sum does not shine on all of
%them directly

for nn=1:5
    if ds1(nn)<0
       ds1(nn)=0; 
    end
end
da1=zeros(1,5);
Dis1= [0,0,0,0,0,  ds1  ,0,0,0,0,0]';

% for leakage disturbance, we consider simun this period 24 hours:
outtemp_ave=35;
outtemp_peak=7;
    
%ds2=dis2*[2*(outtemp_ave+outtemp_peak*sin(2*pi*k/86400)),2*(outtemp_ave+outtemp_peak*sin(2*pi*k/86400)),2*(outtemp_ave+outtemp_peak*sin(2*pi*k/86400)),(outtemp_ave+outtemp_peak*sin(2*pi*k/86400)),3*(outtemp_ave+outtemp_peak*sin(2*pi*k/86400))];
ds2=dis2*[(outtemp_ave+outtemp_peak*sin(2*pi*k/86400)),(outtemp_ave+outtemp_peak*sin(2*pi*k/86400)),(outtemp_ave+outtemp_peak*sin(2*pi*k/86400)),(outtemp_ave+outtemp_peak*sin(2*pi*k/86400)),(outtemp_ave+outtemp_peak*sin(2*pi*k/86400))];


Dis2= [0,0,0,0,0,  ds2  ,0,0,0,0,0]';

 % for  unmeasured heat gain, we use norman  white rand noise
% Noise change every 5 min 
if rem(k,300)==3
    ds3= dis3*(randn(1,5)+1);
end
%ds3=zeros(1,5);
Dis3= [ds3 ,0,0,0,0,0  ,0,0,0,0,0]';
 

%%
OLS = find(rand(1)<cumsum(P(OLS,:)), 1,'first');  % Select a randpm room according P matrix
C(k,OLS)= 1;

%Normal case 

i_err = i_err+(Desired_Temp - X(OLS,k-1))*Ts; %error intergral
Pro = Kp*(Desired_Temp-X(OLS,k-1)); % Propotional part of controlller

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% For stationary case when theromostat in room 4 %%%
%i_err = i_err+(Desired_Temp - X(4,k-1))*Ts; %error intergral
%pid = Kp*(Desired_Temp-X(4,k-1)); % Propotional part of controlller

%%%%%%
Int= Ki*(i_err);

pid = Pro + Int; %Add Integral part of Controller

% For debugging, (it is not a part of simulation)
if rem(k,1800)==0
   AAA=0; 
end


%Saturation
if pid>0
    pid=0;
elseif pid<-1/AirCireff
    pid=-1/AirCireff;
end
%%%%%%

%Make actuator non-linear
B=Bp;
B(1)=-B(1)*(Temp_hvac - X(1,k-1) );
B(2)=-B(2)*(Temp_hvac - X(2,k-1));
B(3)=-B(3)*(Temp_hvac - X(3,k-1));
B(4)=-B(4)*(Temp_hvac - X(4,k-1));
B(5)=-B(5)*(Temp_hvac - X(5,k-1));
%B=B*10;
X(:,k) = A*X(:,k-1) + B*pid +Dis2+Dis1+Dis3; %calculate next state we not miltiplied by A yet
U(k)=pid;

Disturbance(:,k)= Dis2+Dis1+Dis3; 

%B*pid+
res_temp(k)=X(OLS,k);

time(k)=k/3600;
end
 
%% Plotting
figure
hold off
plot(time,X(1:5,:)','linewidth',1)
%legend('room 1','room 2','room 3','room 4','room 5')
%legend({'room 1','room 2','room 3','room 4','room 5'},'Position',[0.75 0.19 0.1 0.2],'FontSize', 12)

%plot(X')
hold on 
dashline = Desired_Temp*ones(1,Tt);


plot(time,dashline,'color','b' ,'linestyle',':','linewidth',0.1)

plot(time,res_temp,'color','k' ,'linestyle','--','linewidth',1)

xlabel('time (hours)','FontSize', 13)
ylabel('Temperature (^{\circ}C)','FontSize', 13)


set(gca,'FontSize',13)

%legend({'room 1','room 2','room 3','room 4','room 5'},'Position',[0.75 0.19 0.1 0.2],'FontSize', 12)
legend({'room 1','room 2','room 3','room 4','room 5'},'Position',[0.81 0.76 0.1 0.2],'FontSize', 12)

title('(a) Temperature of Rooms ','FontSize', 13)

grid on
%ylim([-2 2])




figure('Renderer', 'painters', 'Position', [400 400 560 240])
plot(time,-U','linewidth',1)

xlabel('time (hours)','FontSize', 15)
ylabel('Mass flow rate (m^3/s)','FontSize', 15)

title('(b) Control Input ','FontSize', 13)
grid on


%%%% Generatre Disturbance moments for H1 & H2 analysis

%Disturbance = [Disturbance;zeros(5,Tt)];  %***** % I comment this line for PI***8
Disturbance=Disturbance(:,3:Tt);


ED=mean(Disturbance');
ED=ED';
Covv=cov(Disturbance');
Covv=Covv(:);
ED2=Covv+kron(ED,ED);

Disturbance = [Disturbance;zeros(1,Tt-2)];

EDi=mean(Disturbance');
EDi=EDi';
Covv=cov(Disturbance');
Covv=Covv(:);
EDi2=Covv+kron(EDi,EDi);

%save ED.mat ED
%save ED2.mat ED2
%save EDi.mat EDi
%save EDi2.mat EDi2


