clc
clear all



% For all equal

 %Pgain = -1.5:0.05*4*2:12;%-0.1:0.05/2:1.9;
 %Dgain = 0.01;
 %Igain = -0.39:0.0001*5*4*10*10*2*2:32;
 
%   Pgain = -1.5:0.05*4*10:12;%-0.1:0.05/2:1.9;
%  Dgain = 0.01;
%  Igain = -0.39:0.0001*5*4*10*10*2*10:32;
%  

% For half to central
 Pgain = -4*3.5+2:3.5:375;
  Igain = -4*6.5+2.8620e-07:6.5:700;
  
  
  
  Kp = 2.4439;
Ki = 2.8620e-07;

Igain = Ki*10;
Pgain = Kp;

  %50 percent to central room, maximum stable Kp
%Pgain =354.75   
%Igain=0.0012

  %75 percent to central room, maximum stable Kp
%Pgain =236.75
%Igain=0.0012
  
% Equally Distributed
%Pgain =850+25+6+3+1.5+.75
%Igain=0.0012

%%%
%Pgain = 2.5297;
%Igain = 2.8620e-07*10;

alpha_1=[5];
alpha_2=1;
  

aa = length(Pgain);
bb = length(Igain);
cc= length(alpha_1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    [Ts,B,P,n,A,L,Desired_Temp,dis1,dis2,dis3]=Model(); 
%[Ts,heaterroom,P,n,A,L,outside_temp]=Model();


Pss=null(P-eye(length(P)));
Pss=Pss/sum(Pss);
%Disturbance generator
%%%%
r=length(P);
EDis1=[zeros(r,1);dis1*20*ones(r,1);zeros(r,1)];
EDis2=[zeros(r,1);dis2*20*ones(r,1);zeros(r,1)];
EDis3=[0.002*ones(r,1);zeros(10,1)];

ED=EDis1+EDis2+EDis3;
ED=[ED;zeros(r,1)];


EDi=[ED;0];

ED2=kron(ED,ED);
EDi2=kron(EDi,EDi);

 for i4 = 1:10
   ED2((i4-1)*length(ED)+i4)= ED2((i4-1)*length(ED)+i4)+ 0.2;
   EDi2((i4-1)*length(EDi)+i4)= EDi2((i4-1)*length(EDi)+i4)+ 0.2;
 end

ED= struct2cell(load('ED.mat'));
ED2= struct2cell(load('ED2.mat'));
EDi= struct2cell(load('EDi.mat'));
EDi2= struct2cell(load('EDi2.mat'));

ED=cell2mat(ED);
ED2=cell2mat(ED2);
EDi=cell2mat(EDi);
EDi2=cell2mat(EDi2);
%%%%
F_1=zeros(aa,bb);
F_2=zeros(aa,bb);
J_Matrix=zeros(aa,bb,cc);
StabilityM = zeros(aa,bb);

Counter = 0;

for i1 = 1:aa
    for i2 = 1:bb
Counter = Counter +1;
 

 if Igain(i2)==0
     EDs=ED;
     ED2s=ED2;  
 else
     EDs=EDi;
     ED2s=EDi2;
 end


%[HH1,HH2,r,GSz,H12]=PIH1H2Func3(A,P,B,Pgain(i1),Igain(i2),EDs);
[HH1,HH2,r,GSz,H12]=PIH1H2Func3Stationary(A,P,B,Pgain(i1),Igain(i2),EDs);



%E1 =max(eig(HH1));

%E2 =max(eig(HH2));


EE1 = eig(HH1);
EE2 = eig(HH2);

Max_Eig_H1 = max(EE1);
Max_Eig_H2 = max(EE2);


%%% To display during run
%Disp('attemp' Counter 'from',i1*i2)
disp(['attemp ' num2str(Counter) ' from ' num2str(aa*bb)])

disp(['Pgain = ' num2str(Pgain(i1))])
disp(['Igain = ' num2str(Igain(i2))])

%%%%%

if (Max_Eig_H1> 1 || Max_Eig_H2> 1 )
    
    %J_Matrix(i1,i2)=0;
    StabilityM(i1,i2)=1;
    
    %%% To display during run time
    disp('Not Stable for these gains')
    %%%%
else
 
 %%% To display during run time
    disp('System is Stable for these gains')
 %%%%
 


 
 %H12= kron(EDs,HH1) + kron(HH1,EDs);
 
 HBig=[HH1,zeros(length(HH1),length(HH2));H12,HH2];
 Pss=ones(length(P),1);
 tildeD=[kron(Pss,EDs);kron(Pss,ED2s)];
 
 
 N2=inv(eye(length(HBig))-HBig)*tildeD;

if Igain(i2)==0
    n=GSz-r;
else
    n=GSz-r-1;
end

J = 0;
INPUTCOST=0;
STATECOST=0;
for m = 1:r

    
    %Input Cost
    EX2 = N2(GSz*r+((m-1)*(GSz)^2+(m-1)*(GSz)+m),1);
    StateCost = EX2;
    
    %InputCost
    %EX2= N2(GSz*r+((m-1)*(GSz)^2+(m-1)*(GSz)+m),1);
    
    %
    %EXold2 = N2(GSz*r+((m-1)*(GSz)^2+(n+m-1)*(GSz)+(n+m)),1);
    
    Eacc2= N2(GSz*r+((m)*(GSz)^2),1);
    
    %
    %EXXold = N2(GSz*r+((m-1)*(GSz)^2+(m-1)*(GSz)+m+n),1);
    
    %EXXold2=  N2(GSz*r+((m-1)*(GSz)^2+(n+m-1)*(GSz)+m),1); %To check
    
    
    EXacc= N2(GSz*r+((m-1)*(GSz)^2+m*GSz),1);
    
    EXacc2= N2(GSz*r + (m-1)*GSz^2 +(GSz-1)*GSz+m);
    
    %EXoldacc= N2(GSz*r+((m-1)*(GSz)^2+(n+m)*(GSz)),1);
    
    %EXoldacc2= N2(GSz*r+((m-1)*(GSz)^2+(GSz-1)*(GSz)+m+n));
   
    
    if (Igain(i2)~=0)
        F1=EX2*(Pgain(i1)+Igain(i2)*Ts)^2;
        
        %F1= EX2*(Pgain(i1))^2;
      
    
        %F2 = EXold2*(Dgain/Ts)^2;
    
        F3 = Eacc2*(Igain(i2)*Ts)^2;
    
        %F4 = -2*EXXold*(Dgain/Ts)*(Pgain(i1)+Dgain/Ts+Igain(i2)*Ts);
    
        F5 = 2*EXacc*(Pgain(i1)+Igain(i2)*Ts)*(Igain(i2)*Ts);
        
        %F5 = 2*EXacc*Pgain(i1)*Igain(i2)*Ts;
    
        %F6 = -2*EXoldacc*(Dgain/Ts)*(Igain(i2)*Ts);
    
        %F4 = -(EXXold+EXXold2)*(Dgain/Ts)*(Pgain(i1)+Dgain/Ts);
    
        %F5 = (EXacc+EXacc2)*(Pgain(i1)+Dgain/Ts)*(Igain(i2)*Ts);
    
        %F6 = -(EXoldacc+EXoldacc2)*(Dgain/Ts)*(Igain(i2)*Ts);
        
    InputCost=F1+F3+F5;%+F4+F5+F6;
    KKK=1;
    else
        F1 = EX2*(Pgain(i1))^2;
        %F1=EX2*(Pgain(i1)^2+Dgain^2*(1/Ts)^2+2*Pgain(i1)*Dgain*(1/Ts));
    
        %F2 = EXold2*(Dgain/Ts)^2;
    
    
        %F4 = -2*EXXold*(Dgain/Ts)*(Pgain(i1)+Dgain/Ts);
    
        InputCost=F1;
    end
 
    if (InputCost<0 || StateCost<0)
        KKK=0;
    end
    
    INPUTCOST = INPUTCOST+InputCost;
    STATECOST= STATECOST+StateCost;
    %J = J+alpha_2*InputCost+alpha_1*StateCost;


 
end


F_1(i1,i2)= STATECOST;
F_2(i1,i2)=INPUTCOST;
%J_Matrix(i1,i2)=J;
disp(['State cost is equal to ' num2str(STATECOST) 'Input cost is equal to' num2str(INPUTCOST)])

end

    end
end

for i3=1:cc
   J_Matrix(:,:,i3)=alpha_1(i3)*F_1+alpha_2*F_2 ;



Max_Cost =  max(max(J_Matrix(:,:,i3)));

for i1 = 1:aa
    for i2 = 1:bb
        if (StabilityM(i1,i2)==1)
            J_Matrix(i1,i2,i3)=1.0001*Max_Cost;
        end
    end
end

%%% Draw the plot %%%%%%%%%%%%
y = [Pgain(1);Pgain(aa)];
x = [Igain(1),Igain(bb)];

%%% Delete negative cost function (those are related to unstable gains!)
% [a,b]=find(J_Matrix<0);
% 
% qqq=max(max(J_Matrix));
%  for i=1:length(a)
%      J_Matrix(a(i),b(i))=qqq;
%  end

%%% Plot log10 J matrix to be more visible
%filename = 'J_Matrix.xlsx';
%xlswrite(filename,J_Matrix,1,'A1')

figure
JLog=log10(J_Matrix(:,:,i3));

imagesc(x,y,JLog)
colormap(copper)
colormap(flipud(colormap))

xlabel('Intergration Gain (K_i)','FontSize', 13) % x-axis label
ylabel('Proporinal Gain (K_p)','FontSize', 13) % y-axis label

colorbar

set(gca,'FontSize',13);


min_Cost = min(min(J_Matrix(:,:,i3)));
[a1,a2]=find(J_Matrix(:,:,i3) ==min_Cost);

disp(['Optimal gains are Pgain= ' num2str(Pgain(a1)) ' and Igain = ' num2str(Igain(a2)) ])
disp(['Those Gain give cost fuction' num2str(J_Matrix(a1,a2)) 'for = ' num2str(alpha_1(i3)) ])

end

%J_Matrix = xlsread('J_Matrix.xlsx','A1:AY106')
%[q1,q2]=find(J_Matrix<0)


%Optimal gains are Pgain= 2.5 and Igain = 7.81 
%Those Gain give cost fuction65454.9266for = 1000

%Optimal gains are Pgain= 5.5 and Igain = 0.01
%Those Gain give cost fuction64686.3674for = 10