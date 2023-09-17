%% Nuf2RV_ijk
% Takes Nuf to RV in ijk
% Inputs: Nuf, e, a, i, omega, w
% Outputs: R_ijk, V_ijk

function [R_ijk, V_ijk] = Nuf2RV_ijk(Nuf,e,a,i,omega,w)

mu = 398600.5;

% 
R = (a*((1-(e^2)))) / (1 + (e*cosd(Nuf)));
R_pqw = [R*cosd(Nuf) ; R*sind(Nuf); 0];

P = a*(1-(e)^2);
V_pqw = sqrt(mu/P).*[-sind(Nuf); (e+(cosd(Nuf))); 0 ];

Rot3_N_omega = [cosd(-omega) sind(-omega) 0 ;...
                -sind(-omega) cosd(-omega) 0; ... 
                0 0 1];
            
Rot1_N_i = [1 0 0; ... 
            0 cosd(-i) sind(-i);...
            0 -sind(-i) cosd(-i)];
        
Rot3_N_w = [cosd(-w) sind(-w) 0 ;...
                -sind(-w) cosd(-w) 0; ... 
                0 0 1];

R_ijk = Rot3_N_omega*Rot1_N_i*Rot3_N_w*R_pqw;
V_ijk = Rot3_N_omega*Rot1_N_i*Rot3_N_w*V_pqw;

fprintf('\n FINAL RANGE TRACKING SITE LOCATION:\n') 
fprintf('R_ijk = [%f I + %f J + %f K] km \n',R_ijk(1),R_ijk(2),R_ijk(3))
fprintf('V_ijk = [%f I + %f J + %f K] km/s',V_ijk(1),V_ijk(2),V_ijk(3))




