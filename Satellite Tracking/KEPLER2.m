%% Keplers Type II
% Inputs: Nui, TOF, e, a
% Outputs: Ei, Mi, Mf, Ef, Vf

function [Ei,Mi,Mf,Ef,Nuf] = KEPLER2(Nui,TOF,e,a)

mu = 398600.5;

Ei = acos( (e + cosd(Nui))/(1+(e*cosd(Nui))) );

if Nui > 180 
    Ei = (2*pi)-Ei;
end

Mi = Ei - (e*sin(Ei));  %Radians

n = sqrt(mu/(a^3));

Mf = Mi + (n*TOF);      

while Mf >= 2*pi
    Mf = Mf - 2*pi;
end

Ef = Mf;
while abs(Ef - (Ef + (Mf - Ef + e*sin(Ef))/(1-(e*cos(Ef)))) ) > .00001
    Ef = Ef + ((Mf - Ef + e*sin(Ef))/(1-(e*cos(Ef))));
end

Nuf = acosd( (cos(Ef)-e)/(1-(e*cos(Ef))) ); 
if Ef > pi 
    Nuf = 360-abs(Nuf);
end
