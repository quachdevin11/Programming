% Euler's Method with given a and e to plot R vector of an orbit
% Devin Quach 

function [R_x, R_y] = Euler(a,e)

% Euler's Method
mu = 3.986005*(10^5);               %Constants 
Rp = a*(1-e);                       %Constants
T = 2*pi*sqrt((a^3)/mu);            %Constants

Rdot_x(1) = 0;                      %Initial Condition for X Velocity 
Rdot_y(1) = sqrt((2*mu/Rp)-(mu/a)); %Initial Condition for Y Velocity 
R_x(1) = Rp;                        %Initial Condition for X distance   
R_y(1) = 0;                         %Initial Condition for Y distance

dt = T/100000;                      %Small time step to get one period 
for n = 1:100000 
    R = sqrt((R_x(n))^2 + (R_y(n))^2);     %Equ. for R Vector
    Rddot_x = (-mu/(R^3))*(R_x(n));        %Second ODE of X Comp. of R with indexing
    Rddot_y = (-mu/(R^3))*(R_y(n));        %Second ODE of Y Comp. of R with indexing
  
    Rdot_x(n+1) = Rdot_x(n) + dt*Rddot_x;   %First ODE of X Comp. of R with indexing
    Rdot_y(n+1) = Rdot_y(n) + dt*Rddot_y;   %First ODE of Y Comp. of R with indexing
    
    R_x(n+1) = R_x(n) + dt*Rdot_x(n);       %Distance equ of X component with indexing
    R_y(n+1) = R_y(n) + dt*Rdot_y(n);       %Distance equ of Y component with indexing
    
end

end
