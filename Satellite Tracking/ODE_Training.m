%% MAE 4410_ODE_Training

clc
clear all 

sym y 
ODEFUN = (-2*y) + 2 - (exp(-4*t)); 
t = [0 5];
x0 = 1;
ODE45(ODEFUN, t,x0)
