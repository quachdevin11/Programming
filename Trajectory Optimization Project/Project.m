%% Variation of Extremals with PMP for Low Thrust from LEO to GEO
% Devin Quach
clc 
clear all
close all

mu = 398600.5;                          % In km^3 /s^2
% === Desired Final State === %
TOF = 30002; %9*60*60;                          % 9 Hours In Seconds
r_final = 35786+6378;                   % Final Radius: kMeters
u_final = 0 ;                           % Final Radial Velocity:km/s 
v_final = sqrt(mu/r_final);             % Final Horizontal Velocity: km/s

% === Physical Constants === %
mdot = 0.1;                          % In kg/s Assumed to be Constant 
mo = 3000; %11600;                  % In kg (Mass of Hubble Telescope)


% === Computational Variables ==== %
NumTSteps = TOF;                         % # of TimeSteps                              
Deltat = TOF/NumTSteps;                  % Time Increments for # of time steps till TOF
t = 0:Deltat:TOF;
tol = 0.001;                              % End Condition
NumXVars = 3;

% === Weighting/Constants === % 
tau = 2;                              % Weighting on New Control
del = 1e-3;                              % Pertubation Size for Inital Costates 

% === Preallocating === %
X = ones(NumTSteps+1,NumXVars);         % Preallocating State Vector
a = zeros(NumTSteps+1,NumXVars);        % Preallocating EOM Vector
P = ones(NumTSteps+1,NumXVars)*(10e3);   % Initial Guess for P & Preallocating
Pdot = zeros(NumTSteps+1,NumXVars);     % Preallocating Costate Dot Vector
P(1,:) = [5-4.6833462165128,-857.693225185878,-4804.96591508175];  %From Dr.T's max radius problem 
% P(1,:) = [-1.56004780670084,2461.87889706418,-3223.51891396842];

% Assumption Givens for in a Vacumm (Merlin 1D) 
% === Physical Givens === %
r = 160+6378 ;                          % Initial Radius: kMeters
u = 0 ;                                 % Initial Radial Velocity: km/s
v = sqrt(mu/r);                         % Initial Horizontal Velocity: km/s
X(1,:) = [r, u, v];                     % Initial State Vector 

Thigh =0.001;% 981/mo;%981.5;                          % Higher Bound Constraint of Thrust Acceleration in km/s^2
Tlow = 0;                          % Lower Bound Constraint of Thrust in kN

iter = 1;                               % Start of iterations
while norm( X(end,:) - [r_final, u_final, v_final] ) > tol 
        % === Reduced Diff. Eqs & Propogate Costate Guess Forward to Find Final States=== %
        for k = 1:1:NumTSteps                   % Uses initial guess for P to perturb EOM vector and Pdot, then use Euler's Iteration to solve for a new X and P            
                
                
            u(k,2) = atan( P(k,2)/P(k,3) );                                                                % Control for Phi    
            
            if ((P(k,2)*sin(u(k,2)))/(mo-(mdot*t(k)))) + ((P(k,3)*cos(u(k,2)))/(mo-(mdot*t(k)))) > 1     % PMP Piecewise Control for T
               u(k,1) = Tlow;
            elseif  ((P(k,2)*sin(u(k,2)))/(mo-(mdot*t(k)))) + ((P(k,3)*cos(u(k,2)))/(mo-(mdot*t(k)))) < -1 
               u(k,1) = Thigh;
            else %((P(k,2)*sin(u(k,2)))/(mo-(mdot*t(k)))) + ((P(k,3)*cos(u(k,2)))/(mo-(mdot*t(k)))) < 1 && ((P(k,2)*sin(u(k,2)))/(mo-(mdot*t(k)))) + ((P(k,3)*cos(u(k,2)))/(mo-(mdot*t(k)))) > -1
               u(k,1) = 0; 
            end
                
            a(k,:) = [ X(k,2), ...
                (X(k,3)^2/X(k,1)) - (mu/(X(k,1))^2) + ((u(k,1)*sin(u(k,2)))/(mo-(mdot*t(k)))),...
                (-(X(k,3)*X(k,2))/(X(k,1))) + ((u(k,1)*cos(u(k,2)))/(mo-(mdot*t(k)))) ];                     % EOM Vector
            
            X(k+1,:) = a(k,:)*Deltat + X(k,:);                                                              % New State Vector using Forward Integration
            
            Pdot(k,:) = [ ((P(k,2)*(X(k,3))^2)/(X(k,1)^2)) - ((P(k,3)*X(k,3)*X(k,2))/(X(k,1)^2)) - ((2*P(k,2)*mu)/(X(k,1)^3)),...
                -P(k,1) + ((P(k,3)*X(k,3))/(X(k,1))),...
                -((2*P(k,2)*X(k,3))/X(k,1)) + ((P(k,3)*X(k,2))/(X(k,1))) ];                                 % CoState Dot Vector
            
            P(k+1,:) = Pdot(k,:)*Deltat + P(k,:);                                                           % New CoState Vector using Forward Integration
        end
        
        % === Part C: Influence Fcns === % Using Finite Difference
        
        for m = 1:1:3                                               % Column Wise  
            Ppert(1,:) = P(1,:);
            Xpert(1,:) = X(1,:);    
            Ppert(1,m) = (P(1,m)+del); 
                for k = 1:1:NumTSteps     % Uses initial guess for P to make perturb a vector, then use Euler's Iteration to solve for a new X 
          
                    u(k,2) = atan( Ppert(k,2)/Ppert(k,3) );                                                              % Control for Phi                
                    
                    if ((Ppert(k,2)*sin(u(k,2)))/(mo-(mdot*t(k)))) + ((Ppert(k,3)*cos(u(k,2)))/(mo-(mdot*t(k)))) > 1     % PMP Piecewise Control for T
                        u(k,1) = Tlow;
                    elseif  ((Ppert(k,2)*sin(u(k,2)))/(mo-(mdot*t(k)))) + ((Ppert(k,3)*cos(u(k,2)))/(mo-(mdot*t(k)))) < -1 
                        u(k,1) = Thigh;
                    else
                        u(k,1) = 0;
                    end
                    
            apert(k,:) = [ Xpert(k,2), ...
                (Xpert(k,3)^2/Xpert(k,1)) - (mu/(Xpert(k,1))^2) + ((u(k,1)*sin(u(k,2)))/(mo-(mdot*t(k)))),...
                (-(Xpert(k,3)*Xpert(k,2))/(Xpert(k,1))) + ((u(k,1)*cos(u(k,2)))/(mo-(mdot*t(k)))) ];                     % EOM Vector
            
            Xpert(k+1,:) = apert(k,:)*Deltat + Xpert(k,:);                                                              % New State Vector using Forward Integration
            
            Pdotpert(k,:) = [ ((Ppert(k,2)*(Xpert(k,3))^2)/(Xpert(k,1)^2)) - ((Ppert(k,3)*Xpert(k,3)*Xpert(k,2))/(Xpert(k,1)^2)) - ((2*Ppert(k,2)*mu)/(Xpert(k,1)^3)),...
                -Ppert(k,1) + ((Ppert(k,3)*Xpert(k,3))/(Xpert(k,1))),...
                -((2*Ppert(k,2)*Xpert(k,3))/Xpert(k,1)) + ((Ppert(k,3)*Xpert(k,2))/(Xpert(k,1))) ];                                 % CoState Dot Vector                                              
                    
            Ppert(k+1,:) = Pdotpert(k,:)*Deltat + Ppert(k,:);                                                               % Perturbed CoState Vector using Forward Integration
                end
            Px(:,m) =   (Xpert(end,:)-X(end,:))' / del ;                     % Column of Px                                       
        end

        P(1,:) = P(1,:) - tau*( (Px^-1)*( X(end,:) - [r_final, u_final, v_final] )')';        % Propogating through to get New Costate Vector
iter = iter+1;
end
u(end+1,:) = [0 0];

z(:,1) = X(:,1); 
z(:,2) = X(:,2); 
z(:,3) = X(:,3); 
z(:,4) = u(:,1); 
z(:,5) = u(:,2)*(180/pi);
% === Plotting Results === %

figure(1)
subplot(3,1,1)
plot(t,X(:,1))
title('States: Radial Distance')
xlabel('t (s)')
ylabel('Radius (km)')
axis([0 TOF r-10 r_final+10])
grid on

subplot(3,1,2)
plot(t,X(:,2))
title('States: Radial Velocity')
xlabel('t (s)')
ylabel('Radial Velocity (km/s)')
axis([0 TOF min(X(:,2)) max(X(:,2)) ])

subplot(3,1,3)
plot(t,X(:,3))
title('States: Horizontal Velocity')
xlabel('t (s)')
ylabel('Horizontal Velocity (km/s)')
axis([0 TOF min(X(:,3)) max(X(:,3))] )

figure(2)
subplot(2,1,1)
plot(t(1:end),u(:,1))
title('Controls: Thrust')
xlabel('t (s)')
ylabel('Thrust (kN)')
grid on
axis([0 TOF -Thigh Thigh ])

subplot(2,1,2)
plot(t(1:end),u(:,2)*180/pi)
title('Controls: Phi')
xlabel('t (s)')
ylabel('Phi (Degrees)')
axis([0 TOF -100 100] )
