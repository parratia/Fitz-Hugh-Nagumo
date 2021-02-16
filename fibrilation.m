function [Nt,L,nu,c2,d,t,sigma,alfa,u0,r0,stop] = fibrilation()
Nt = 1000; %Time steps
L=2; %Square size
nu = 31; %Constant nu
c2 = 3; %Constant delta
d = 0.5; %Constant d
u0 = 1; %Initial value for potential
r0 = 0; %Initial value for recovery variable
stop = Nt; %At this step pulses are not created anymore
sigma=0.02; %Conductivity constant
alfa=0.05; %Constant alpha
sigma = @(x,y) sigma*exp(-(x-L/2)^2-(y-L/2)^2); %Conductivity function

t=20; %Period for pulses
end