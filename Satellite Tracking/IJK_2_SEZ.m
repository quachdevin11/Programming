%% IJK to SEZ Position Vector's 
% Inputs: [rof,Azf,Elf] = IJK_2_SEZ(Lf,GST,Hf,R_ijk,TOF,LEf)
% Outputs: rof, Azf, Elf

function [rof,Azf,Elf] = IJK_2_SEZ(Lf,GST,Hf,R_ijk,TOF,LEf)

GSTf = GST*15 + (TOF*15/3600);      %(hr)*(degrees/hour) + sec*(hr/sec)*(15deg/hr)
LSTf = GSTf + LEf;                  %degrees

while LSTf >= 360
    LSTf = LSTf - 360;
end

% Constants
a_not = 6378.137;
e = 0.0818;

% Position
x = (((a_not)/sqrt(1 - (e^2)*(sind(Lf))^2))+Hf)*cosd(Lf);
z = (((a_not*(1-(e^2)))/sqrt(1 - (e^2)*(sind(Lf))^2))+Hf)*sind(Lf);
R_site = [x*cosd(LSTf); x*sind(LSTf); z ]; 

ro_ijk = R_ijk - R_site;

COLATf = 90 - Lf; 

ROT2_N_COLATf = [cosd(-COLATf) 0 -sind(-COLATf);...
                 0 1 0;...
                 -sind(COLATf) 0 cosd(-COLATf)];
ROT3_N_LSTf = [cosd(-LSTf) sind(-LSTf) 0; ...
                -sind(-LSTf) cosd(-LSTf) 0;...
                0 0 1];
            
ro_sez = (ROT3_N_LSTf*ROT2_N_COLATf)'*ro_ijk;

rof = norm(ro_sez);
Elf = asind(ro_sez(3)/rof);
Azf = acosd( (-ro_sez(1)) / (rof*cosd(Elf)) );

if ro_sez(2) < 0
    Azf = 360-Azf;    
end

fprintf('\n Range = %f km \n Azimuth = %f degrees\n Elevation = %f degrees \n',rof,Azf,Elf)
fprintf('LST = Local Sidereal Time =  %f degrees\n',LSTf);
