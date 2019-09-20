function J = CostFunct(KGAIN)
Pgain=KGAIN(1)
Igain=KGAIN(2)

alpha_1=1;
alpha_2=1;

[Ts,B,P,n,A,L,Desired_Temp,dis1,dis2,dis3]=Model(); 


ED= struct2cell(load('ED.mat'));
ED2= struct2cell(load('ED2.mat'));
EDi= struct2cell(load('EDi.mat'));
EDi2= struct2cell(load('EDi2.mat'));

ED=cell2mat(ED);
ED2=cell2mat(ED2);
EDi=cell2mat(EDi);
EDi2=cell2mat(EDi2);

 if Igain==0
     EDs=ED;
     ED2s=ED2;  
 else
     EDs=EDi;
     ED2s=EDi2;
 end


[HH1,HH2,r,GSz,H12]=PIH1H2Func3(A,P,B,Pgain,Igain,EDs);

EE1 = eig(HH1);
EE2 = eig(HH2);

Max_Eig_H1 = max(EE1);
Max_Eig_H2 = max(EE2);


if (Max_Eig_H1> 1 || Max_Eig_H2> 1 )
    
    %J_Matrix(i1,i2)=0;
    StabilityM=1;
    
    %%% To display during run time
    disp('Not Stable for these gains')
    %%%%
    J=10^4;
else
 
 %%% To display during run time
    disp('System is Stable for these gains')
    StabilityM=0;
 %%%%
    HBig=[HH1,zeros(length(HH1),length(HH2));H12,HH2];
 Pss=ones(length(P),1);
 tildeD=[kron(Pss,EDs);kron(Pss,ED2s)];
 
 
 N2=inv(eye(length(HBig))-HBig)*tildeD;

if Igain==0
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
   
    
    if (Igain~=0)
        F1=EX2*(Pgain+Igain*Ts)^2;
        
        %F1= EX2*(Pgain(i1))^2;
      
    
        %F2 = EXold2*(Dgain/Ts)^2;
    
        F3 = Eacc2*(Igain*Ts)^2;
    
        %F4 = -2*EXXold*(Dgain/Ts)*(Pgain(i1)+Dgain/Ts+Igain(i2)*Ts);
    
        F5 = 2*EXacc*(Pgain+Igain*Ts)*(Igain*Ts);
        
        %F5 = 2*EXacc*Pgain(i1)*Igain(i2)*Ts;
    
        %F6 = -2*EXoldacc*(Dgain/Ts)*(Igain(i2)*Ts);
    
        %F4 = -(EXXold+EXXold2)*(Dgain/Ts)*(Pgain(i1)+Dgain/Ts);
    
        %F5 = (EXacc+EXacc2)*(Pgain(i1)+Dgain/Ts)*(Igain(i2)*Ts);
    
        %F6 = -(EXoldacc+EXoldacc2)*(Dgain/Ts)*(Igain(i2)*Ts);
        
    InputCost=F1+F3+F5;%+F4+F5+F6;
    KKK=1;
    else
        F1 = EX2*(Pgain)^2;
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

J = alpha_2*INPUTCOST+alpha_1*STATECOST 
 
 
end


















end
