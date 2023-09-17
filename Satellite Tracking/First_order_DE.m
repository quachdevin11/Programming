%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program will solve an example first order differential equation  % 
% using two methods:  Euler's integration method and Matlab's ODE 45    %    
% function.  The differential equation is:                              %
%                          y'=2-2y-exp(-4t)                             %
%                              y(1) = 1
%                             0 <= t <= 5 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Written by Lynnane George 1/26/20                      %
%                   For MAE 4410, Astrodynamics                         %
%                           Spring 2020                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

% Set initial conditions
y(1) = 1;
t(1) = 0;


% Euler's Method
% set a small time step 
dt = 0.01;
for n = 1:500
    t(n+1)= t(n)+dt;
ydot=2-2*y(n)-exp(-4*t(n+1));
y(n+1)=y(n)+dt*ydot;
end
figure(1);
plot(t,y); title('Eulers Method');

% ODE45 Method
time=[0,5];
initz=y(1);
[time,z]=ode45(@f2, time, initz);
figure(2);
plot(time,z(:,1),'x'); title('ODE 45 Method');

% Actual response
ya(1) = 1;
ta(1) = 0;
for n = 1:500
    ta(n+1)= ta(n)+dt;
ya(n+1)=1+.5*exp(-4*ta(n+1))-.5*exp(-2*ta(n+1));
end
figure(3);
plot(ta,ya,'+'); title('Actual Response')

% ODE45 function
function dz=f2(time,z)
dz=[2-2*z-exp(-4*time)];
end


