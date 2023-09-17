%% COE: Finds Classical Orbital Elements Given V and R Vectors
% Inputs: (R_Bar,V_Bar)
% Outputs: [V,R,E,a,H_Bar,h,e_Bar,e,nu,i,omega,w] Gives orbit shape as
% 'O'

function [V,R,E,a,h_Bar,h,e_Bar,e,nu,i,n_Bar,n,omega,w] = COE(R_Bar,V_Bar)

        mu = 398600.5;                                      %Constant mu
        V= norm(V_Bar);                                     %Mag of V
        R= norm(R_Bar);                                     %Mag of R
        E = ((V^2)/2) - mu/R;                               %Specific Energy
        a = -mu/(2*E);                                      %Semi-major Axis
        h_Bar = cross(R_Bar,V_Bar);                         %H-Vector
        h = norm(h_Bar);                                    %H Magnitude
        e_Bar = ((cross(V_Bar,h_Bar))./mu) - (R_Bar./R);    %Eccentricity Vector
        e = norm(e_Bar);                                    %Eccentricity
        nu = acosd((dot(e_Bar,R_Bar))/(e*R));               %True Anomoly
        i = acosd(h_Bar(3)/h);                              %Inclination
        n_Bar = cross([0;0;1],h_Bar);                       %Normal Vector
        n = norm(n_Bar);                                    %Magnitude of Normal Vector
        omega = acosd(n_Bar(1)/n);                          %RAAN
        w = acosd(dot(n_Bar,e_Bar)/(n*e));                  %Argument of Perigee 
        
        if dot(R_Bar,V_Bar) < 0                             %QUAD CHECKING NU
            nu = 360 - nu;
        end
        
        if n_Bar(2) < 0                                     %QUAD CHECKING OMEGA
            omega = 360 - omega;
        end
        
        if e_Bar(3) < 0                                     %QUAD CHECKING Arg. Perigee
            w = 360-w;
        end
        
        
        % Assigning Orbits
        
        if e < .001  
            O = 'The Orbit is Circular';
        elseif (.001<e)&&(e<0.999)
            O = 'The Orbit is Elliptical';
        elseif abs((1-e)) < .001
            O = 'The Orbit is Parabolic';
            a = inf;
        else
            O = 'The Orbit is Hyperbolic';
        end
       
%         fprintf('R = [%f %f %f] \n',R_Bar(1),R_Bar(2),R_Bar(3))
%         fprintf('V = [%f %f %f] \n',V_Bar(1),V_Bar(2),V_Bar(3))
%         fprintf('\n-----------------------------------------------\n') 
        
            if e < .001 && i > .001                 %alternate COE argument of latitude
                w = 0;
                nu = 0;
                u = acosd(dot(n_Bar,R_Bar)/(n*R));
                    if R_Bar(3) < 0
                        u = 360 - u;                    %quad check
                    end
            else
                u = 0;
            end
            
            if e > .001 && i < .001                 %alternate COE longitude of perigee
                w = 0;
                omega = 0;
                Pi = acosd(e_Bar(1)/e);
                    if ev(2) < 0                        %quad check
                        Pi = 360 - Pi;
                    end
            else
                Pi = 0;
            end
            
            if e < .001 && i < .001                 %alternate COE true longitude
                omega = 0;
                w = 0;
                nu = 0; 
                l = acosd(R_Bar(1)/R);
                    if R_Bar(2) < 0                 %quad check
                        l = 360 - l;
                    end
            else
                l = 0;
            end
            
%         fprintf('\n-----------------------------------------------\n ALTERNATE COEs \n\n')    
%         fprintf('u = Arg. of Latitude = %f \n\n',u);
%         fprintf('pI = Longitude of Perigee = %f \n\n',Pi);
%         fprintf('l = True Longitude = %f ',l);
%         fprintf('\n-----------------------------------------------\n') 
%         
        fprintf('\n\n INITIAL COEs \n')
        fprintf('a = Semi-major axis = %f \n',a);
        fprintf('e = Eccentricity = %f \n',e);
        fprintf('i = Inclination = %f \n',i);
        fprintf('Omega = RAAN = %f \n',omega);
        fprintf('w = Argument of Perigee = %f \n',w);
        fprintf('nu = True Anomaly = %f ',nu);
        fprintf('\n\n ORBIT \n')
        disp(O)


  

        
