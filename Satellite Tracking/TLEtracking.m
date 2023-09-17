%% Takes TLE elements and gives initial & final satelite COES 
% [initialCOE,finalCOE]= TLEtracking(ei,ii,omegai,wi,Mi,ni,n_doti,TOF)

function [initialCOE,finalCOE] = TLEtracking(ei,ii,omegai,wi,Mi,ni,n_doti,TOF,LWf,LEf,Lf,Hf,GST)
% Constants
mu = 398600.5;                          % km^3/s^2

ai = (((mu)./(ni.^2))).^(1/3);          % Semi-major Axis

while Mi >= 2*pi
    Mi = Mi - 2*pi;
end

Ei = 0;
while abs(Ei - (Ei + (Mi - Ei + ei*sin(Ei))/(1-(ei*cos(Ei)))) ) > .00001
    Ei = Ei + ((Mi - Ei + ei*sin(Ei))/(1-(ei*cos(Ei))));
end

nui = acosd( (cos(Ei)-ei)/(1-(ei*cos(Ei))) );         % True Anomaly Initial
if Ei > pi 
    nui = 360-abs(nui);
end

% Inital COES
        fprintf('\n\n SATELITEs INITIAL COEs \n')
        fprintf('a = Semi-major axis = %f \n',ai);
        fprintf('e = Eccentricity = %f \n',ei);
        fprintf('i = Inclination = %f \n',ii);
        fprintf('Omega = RAAN = %f \n',omegai);
        fprintf('w = Argument of Perigee = %f \n',wi);
        fprintf('nu = True Anomaly = %f \n',nui);
        initialCOE = [ai; ei; ii; omegai; wi; nui];
        
 % Final COEs
J2 = 1082.64*(10^(-6));                 % 
Re = 6378;                              % km
Po = ai*((1-ei)^2);
n_bar =  (1+( (3/2)*(J2)*((Re/Po)^2)*((1-(ei^2))^.5)*(1 - (3/2)*((sind(ii))^2)) ))*ni;

 a_dot = (-2*n_doti*ai)/(3*ni); 
 e_dot = (-2/3)*(1-ei)*(n_doti/ni); 
 omega_dot = ((-(3/2)*J2*((Re/Po)^2)*cosd(ii))*n_bar)*(180/pi); 
 w_dot = ( (3/2)*J2*((Re/Po)^2)*(2-((5/2)*(sind(ii))^2)) )*n_bar*(180/pi);

 af = ai + (a_dot*TOF);
 ef = ei + (e_dot*TOF);
 omegaf = omegai + (omega_dot*TOF);
 wf = wi + (w_dot*TOF);
 [Ei,Mi,Mf,Ef,nuf] = KEPLER2(nui,TOF,ei,ai);
 
 fprintf('\n SATELITLEs FINAL COEs \n')
        fprintf('a = Semi-major axis = %f \n',af);
        fprintf('e = Eccentricity = %f \n',ef);
        fprintf('i = Inclination = %f \n',ii);
        fprintf('Omega = RAAN = %f \n',omegaf);
        fprintf('w = Argument of Perigee = %f \n',wf);
        fprintf('nu = True Anomaly = %f \n',nuf);
        finalCOE = [af; ef; ii; omegaf; wf; nuf];
 
% Tracking Info From Station
[R_ijk, V_ijk] = Nuf2RV_ijk(nuf,ef,af,ii,omegaf,wf);
        
[rof,Azf,Elf] = IJK_2_SEZ(Lf,GST,Hf,R_ijk,TOF,LEf);
        
w_dot = w_dot*86400;
omega_dot = omega_dot*86400;
        fprintf('\nw_dot = Change in Argument of Perigee = %f degrees/day\n',w_dot);
        fprintf('Omega_dot = Change in RAAN = %f degrees/day \n',omega_dot);

end




