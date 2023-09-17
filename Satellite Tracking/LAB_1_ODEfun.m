%% ODE Function

function xdot = LAB_1_ODEfun(t,x)
    R = sqrt((x(1))^2 + (x(2))^2);          %R vector
    mu = 3.986005*(10^5);                   %Constant 
    xdot(1,1) = x(3);                       %Labeling of X state vector, xdot
    xdot(2,1)= x(4);                        %Labeling of X state vector, ydot
    xdot(3,1)= (-mu/(R^3))*x(1) ;           %Labeling of X state vector, xdoubledot
    xdot(4,1)= (-mu/(R^3))*x(2);            %Labeling of X state vector, ydoubledot
end

