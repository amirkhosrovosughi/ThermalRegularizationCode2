function [HH1,HH2,r,GSz,H12]=PIDH1H2Func3(A,P,heaterroom,Pgain,Igain,EDs);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
r = length(P); 
n = length (A);


%%%% Added for Handling Igain=0 issue %%%%%%%%%%%%%
if (Igain==0)
    GSz =n;
else
    GSz =n+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


O = zeros (r);
I = eye(r);
Znr=zeros(n,r);

%StabilityM = 0;


 T = [];
 TH2 = [];
 
 for ii  = 1:r
 % Controller P gain
%        a = zeros(1,r);
%        if (ii<=r)
%      a(ii)= 1;
%        end
%      S2 = a;
%     for k = 1:r-1
%         S2 = [S2;a];    
%     end
    S2=zeros(n);
    S2(:,ii)=ones(n,1);
    
    
    S1 = zeros(n);
    %S1(heaterroom,heaterroom)= 1;
    for nn=1:r
        S1(nn,nn)= heaterroom(nn);
    end 
   
    C2 = -Pgain*S1*S2; %%%%%%%%%%%%%%%%%%%%%%%%%Sign Pgain??
   % fr
      %C22 =[C21,zeros(r-1,1)];
      %C2 = [C22;-1*A(r,:)];
      
      
       % Big Matrix    
       Ab = A; 
     %Ab(2*r,r)=0;
     
%%%% Added for Handling Igain=0 issue %%%%%%%%%%%%%
if (Igain~=0)
     Ab = [Ab,zeros(n,1)];
     Ab = [Ab;zeros(1,n+1)];
     Ab((n+1),ii)= 1;
     Ab(n+1,n+1)=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     %Propertional controller
     Cpb = C2;%[C2,Znr,;Znr',O];

     
%%%% Added for Handling Igain=0 issue %%%%%%%%%%%%%
if (Igain~=0)
     Cpb = [Cpb,zeros(n,1)];
     Cpb = [Cpb;zeros(1,n+1)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     
%       % Controller D
%       
%       S1b = [S1,Znr;Znr',O];
%       K1=[S2(1:r,1:r);zeros(n-r,r)];
%       K2=[S2,-K1];
%       K3=[Znr',O];
%       S2d = [K2;K3];
%       Cdb = -Dgain*S1b*S2d;
% 
% %%%% Added for Handling Igain=0 issue %%%%%%%%%%%%%
% if (Igain~=0)
%      Cdb = [Cdb,zeros(n+r,1)];
%      Cdb = [Cdb;zeros(1,n+r+1)];
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %Controller I
%%%% Added for Handling Igain=0 issue %%%%%%%%%%%%%
if (Igain==0)
    %Cib = [O,O;O,O];
    Cib=zeros(n);
else
      S1bi = [S1,zeros(n,1)];
     S1bi = [S1bi;zeros(1,n+1)];
      
      S2i1 = zeros(n+1,n);
      S2i1 = [S2i1,ones(n+1,1)];
      
      S2i2 = S2;%[S2,Znr;Znr',O];
      S2i2 = [S2i2,zeros(n,1)];
      S2i2 = [S2i2;zeros(1,n+1)];
      
      Cib = -Igain*S1bi*(S2i1+S2i2);
      
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
      %clear S2 for next round
      %clear S2
        
    
    T1 = Ab+Cpb+Cib;
    
     
  
    T = [T;T1];
  
    q2 = size(T1);
    z = q2(1);
    %--------This part is added for H2 Calculation
        T2 = [];
    for j = 1:z
       Hc = [];
        for i = 1:z
        
            Hc1 = T1(i,j)*T1;
            Hc = [Hc;Hc1];    
        end
       T2 = [T2,Hc]; 
      % clear Hc
    end
    
    TH2 = [TH2;T2];
  
    %------- end of calculation for obtaining H2
    
 end

%q=q1;

HH1 = [];
HH2 = [];
H12 = [];

 a2 = size(P);
 a1 = a2(1);
for j = 1:a1
    Hc = [];
    HCC = [];
    HH12 = [];
    for i = 1:a1
        
   %    Hc1 = P(i,j)*T((i-1)*z+1:(i-1)*z+z,:);
   %   Hc2 = P(i,j)*TH2((i-1)*z*z+1:(i-1)*z*z+z*z,:);

       Hc1 = P(i,j)*T((j-1)*z+1:(j-1)*z+z,:);
       Hc2 = P(i,j)*TH2((j-1)*z*z+1:(j-1)*z*z+z*z,:);
       
       Hc3=kron(EDs,Hc1)+kron(Hc1,EDs);
       
       Hc = [Hc;Hc1];
       HCC = [HCC;Hc2];
       HH12 = [HH12;Hc3];
    end

       HH1 = [HH1,Hc]; 
       HH2 = [HH2,HCC];
       H12=[H12,HH12];
    
end






end

