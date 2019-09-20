function [Ts,Bp,P,n,A,L,Desired_Temp,dis1,dis2,dis3,Temp_hvac]=Model()

Ts = 1;      


%B=[0,1,0,0,0,0]';

c=.99975;


n = 5; %Number of rooms
P = (n*c-1)/(n-1)*eye(n)+((1-c)/(n-1))*ones(n,n); %Matkov chain Matrix
                                                  % Person moves according
                                                 
                                                  
%%% Physical parameters  
Li=0.12; % Thickness of wall


Cpa=1007; % Air specific termal capacity

Ma=1.1839*5*5*3; %Air mass, Density of air in 25 centigrad 1.1839, Check it!
Ca=Cpa*Ma; % Air termal capacity

Mw=5*3*Li*500;  %Density: Brick 2000 
                %Wall mass, 2370 Density of conrete
                %Plywood 10.4
                %wood 400-700 ==> we select wood
                
Cw=1200*Mw;      %thermal capacity of each wall  
                 %check the value  %brick 800
                             %concrete 1000
                             %timber 1200-2000

Ar=3*5; %area of each wall
hi=10;  %check the value  % 5 to 37 % if we choose smaller values (0.1)
ho=5;  %check the value
k=0.14;  %thermal conductivity of wall  % concrete 1.13
                                        % timber 0.14                 
                                        %brick 0.73
 
hi=10;  
ho=5;                                       
                                        
alpha=0.7;

Desired_Temp = 21;
Temp_hvac=12;
%Temp_hvac=-5;
%Temp_hvac=53.3;

Rcov=1/(hi*Ar);
Rw=Li/(k*Ar);

Rcov_out=1/(ho*Ar);
Ro = Rw/2 + Rcov_out;

                                          
%%%% Construct Laplacian matrix (L)

% a to e belongs to air temperature
a2=1/(Ca*0.5*(Rcov+Rw/2));
a3=1/(Ca*(Rcov+Rw/2));
a4=1/(Ca*(Rcov+Rw/2));
a1=-(a2+a3+a4);

b2=1/(Ca*0.5*(Rcov+Rw/2));
b3=1/(Ca*(Rcov+Rw/2));
b4=1/(Ca*(Rcov+Rw/2));
b1=-(b2+b3+b4);

c2=1/(Ca*0.5*(Rcov+Rw/2));
c3=1/(Ca*(Rcov+Rw/2));
c4=1/(Ca*(Rcov+Rw/2));
c1=-(c2+c3+c4);

d2=1/(Ca*(Rcov+Rw/2));
d3=1/(Ca*(Rcov+Rw/2));
d4=1/(Ca*(Rcov+Rw/2));
d5=1/(Ca*(Rcov+Rw/2));
d1=-(d2+d3+d4+d5);


e2=1/(Ca*(Rcov+Rw/2)/3);
e3=1/(Ca*(Rcov+Rw/2));
e1=-(e2+e3);

% modified here from old version Cw must be number_of_walls*Cw
% f to j belongs to exterior walls
f1=1/(Cw*(Rcov+Rw/2));
f2=-f1-1/(Cw*(Ro));

g1=1/(Cw*(Rcov+Rw/2));
g2=-g1-1/(Cw*(Ro));

h1=1/(Cw*(Rcov+Rw/2));
h2=-h1-1/(Cw*(Ro));

i1=1/(Cw*(Rcov+Rw/2));
i2=-i1-1/(Cw*(Ro));

j1=1/(Cw*(Rcov+Rw/2));
j2=-j1-1/(Cw*(Ro));

% k to p belongs to interior walls

k1=1/(Cw*(Rcov+Rw/2));
k2=1/(Cw*(Rcov+Rw/2));
k3=-(k1+k2);

m1=1/(Cw*(Rcov+Rw/2));
m2=1/(Cw*(Rcov+Rw/2));
m3=-(m1+m2);

n1=1/(Cw*(Rcov+Rw/2));
n2=1/(Cw*(Rcov+Rw/2));
n3=-(n1+n2);

o1=1/(Cw*(Rcov+Rw/2));
o2=1/(Cw*(Rcov+Rw/2));
o3=-(o1+o2);

p1=1/(Cw*(Rcov+Rw/2));
p2=1/(Cw*(Rcov+Rw/2));
p3=-(p1+p2);

L = [a1 0 0 0 0   a2 0 0 0 0   a3 0 0 a4 0
     0 b1 0 0 0   0 b2 0 0 0   b3 b4 0 0 0
     0 0 c1 0 0   0 0 c2 0 0   0 c3 c4 0 0
     0 0 0 d1 0   0 0 0 d2 0   0 0 d3 d4 d5
     0 0 0 0 e1   0 0 0 0 e2   0 0 0  0  e3
     
     f1 0 0 0 0   f2 0 0 0 0   0 0 0 0 0
     0 g1 0 0 0   0 g2 0 0 0   0 0 0 0 0
     0 0 h1 0 0   0 0 h2 0 0   0 0 0 0 0
     0 0 0 i1 0   0 0 0 i2 0   0 0 0 0 0
     0 0 0 0 j1   0 0 0 0 j2   0 0 0 0 0
     
     k1 k2 0 0 0   0 0 0 0 0   k3 0 0 0 0
     0 m1 m2 0 0   0 0 0 0 0   0 m3 0 0 0
     0 0 n1 n2 0   0 0 0 0 0   0 0 n3 0 0
     o1 0 0 o2 0   0 0 0 0 0   0 0 0 o3 0
     0 0 0 p1 p2   0 0 0 0 0   0 0 0 0 p3];
 


%%%%

% to this model
                                                  

A = expm (Ts*L); % find discrete model

%Simulate actuator in non-linear way!
%q1=1007*(Temp_hvac- Desired_Temp )/Ca;
q1=Cpa/Ca;

% equal distribution
%Bp=[q1/5 q1/5 q1/5 q1/5 q1/5   0 0 0 0 0   0 0 0 0 0]';


%base
Bp=[q1/8 q1/8 q1/8 q1/2 q1/8   0 0 0 0 0   0 0 0 0 0]';

% 20% more in room 4
%Bp=[q1/10 q1/10 q1/10 6*q1/10 q1/10   0 0 0 0 0   0 0 0 0 0]';

% 50% more in  room 4
%Bp=[q1/16 q1/16 q1/16 12*q1/16 q1/16   0 0 0 0 0   0 0 0 0 0]';

% 100% more in room 4
%Bp=[0 0 0 q1 0   0 0 0 0 0   0 0 0 0 0]';

%Bp=[q1/5 q1/5 q1/5 q1/5 q1/5   0 0 0 0 0   0 0 0 0 0]';
%Bp=[3*q1/20 3*q1/20 3*q1/20 2*q1/5 3*q1/20   0 0 0 0 0   0 0 0 0 0]';
%Bp=[9*q1/40 9*q1/40 9*q1/40 q1/10 9*q1/40   0 0 0 0 0   0 0 0 0 0]';


%%% Calculate disturbances
dis1=alpha*Ar/Cw;
dis2=1/(Cw*(Ro));
dis3=1/Ca;

end
