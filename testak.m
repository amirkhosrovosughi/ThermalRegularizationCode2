Pgain = 2.5297;
%Igain =  1.8723e-10;
Igain = 0.0000541;

KGAIN=[Pgain,Igain]
x = fminsearch(@CostFunct,KGAIN);
x