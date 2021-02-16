function [Nt,L,nu,c2,d,t,sigma,alfa,u0,r0,stop] = normal_heart()
Nt=300; %Time steps

L=1; %Square size
nu=31; %Constant nu
c2 = 3; %Constant delta
d = 0.5; %Constant d
u0 = 1; %Initial value for potential
r0 = 0; %Initial value for recovery variable
stop = 200; %At this step pulses are not created anymore
sigma = 0.01; %Conductivity constant
alfa = 0.1; %Constant alpha
sigma = @(x,y) sigma; %Conductivity function

t=100; %Period for pulses
end