% This function calculates ThetaBN by prolonging B vector from the SC to
% the average bow shock position specified by L, e, x0 parameters in
% cylindrical coordinates; then ThetaBN is calcutated as the angle between
% B and Bow Shock normal angle in place of intersection
%
%
% =====SYNTAX=====
% INPUT:
% 1) MSO spacetraft position in Rm [1x3]
% 2) Magnetic field in MSO [1x3]
%
% OUTPUT:
% 1) ThetaBN in degrees [1x1]
% 2) Distance to intersection in Rm [1x1] (optional)
%
%           Author: Sergey Shuvalov (shuvalovsergei@gmail.com,
%                                    +79153481805)
%           Last change date: 16.11.2017 13:42


function [thetabn, dist] = ThetaBN_func(sc_mso, B)
n = zeros(1, 3);

% Bow Shock parameters (Trotignon, 2006)
L = 2.081;
e = 1.026;
x0 = 0.6;

% Parameters of B line in cylindrical coordinates
k = sqrt(B(2)^2+B(3)^2)/B(1);
b = sqrt(sc_mso(2)^2+sc_mso(3)^2) - sc_mso(1)*sqrt(B(2)^2+B(3)^2)/B(1);

% Calculation of B and BS intersection
A_eq = k^2 - e^2 + 1;
B_eq = 2*k*b + 2*e*L - 2*x0 + 2*e^2*x0;
C_eq = b^2 - L^2 - 2*e*L*x0 - e^2*x0^2 + x0^2;

D_eq = B_eq^2 - 4*A_eq*C_eq;
if(D_eq>0)
        % 2 solutions of 2-order equation
    x_cand = [(-B_eq+sqrt(D_eq))/(2*A_eq) (-B_eq-sqrt(D_eq))/(2*A_eq)];
        % calculating distance between SC and intersection for each solution
    distance = sqrt((x_cand-sc_mso(1)).^2 + (k*x_cand+b-sqrt(sc_mso(2)^2+sc_mso(3)^2)).^2);
        % choosing the closest to the SC intetsection and writing the
        % distance
    xn = x_cand(distance == min(distance));
    dist = min(distance);
        % in case that D==0 we need to delete the second xn
    if(length(xn)>1)
        xn(2)=[];
    end
    
    % Parameters of the equation of normal in cylindrical coordinates
    kn = sqrt((L - e*(xn-x0)).^2 - (xn-x0).^2)/(e*(L-(xn-x0)*e)+xn-x0);
    %bn = sqrt((L - e*(xn-x0)).^2 - (xn-x0).^2) - xn*sqrt((L - e*(xn-x0)).^2 - (xn-x0).^2)/(e*(L-(xn-x0)*e)+(xn-x0));
    
    % Calculation of the components of the normal in MSO
    n(1) = 1/sqrt(kn^2 + 1);
    n(2) = sign(sc_mso(2))*kn/sqrt(2*(kn^2+1));
    n(3) = sign(sc_mso(3))*kn/sqrt(2*(kn^2+1));
    
    thetabn = acos(sum(n.*B)/sqrt(sum(B.^2)))*180/pi;
    if(thetabn>90)
        thetabn = 180-thetabn;
    end
else
        % if there is no intersection
    thetabn = nan;
    dist = nan;
end
end