function [rVeh, rOut, thetaOut, vOut, psiOut, outStatus, Tflight, r_mn] = InBetween(r0, theta0, v0, psi0, dt_min, T_sim_max, MnTinit)

G = 6.67427*(10^(-11));
M_E = 5.9736*(10^24); %[kg]
R_E = 6370000; %[m]
u_E = G*M_E;
M_mn = 7.3477*(10^22); %[kg]
R_mn = 1730000; %[m]
u_mn = G*M_mn;
outStatus = 'error';

%Moon's orbit characteristics:
T_mn = 2.361*(10)^6; %[sec]
e_mn = 0.0549;
a_mn = 384.4*(10)^6;
lts_mn = a_mn*(1-e_mn^2);

%Zone B boundaries:
Re99 = 1.71*(10^8); %[m] 1.723*(10^8)-13000km for histeresis
Rm99 = 3.7*(10^6); %[m] 4*(10^6)-300km for histeresis

% R = (r,thetha) (position is described in angular coordinate system with 
% the origin in the centre of the Earth)
% velocity V = (v, xsi) (xsi - an angle btw the V vector and a circle of 
% radius R centred in the origin; -90 to 90; positive outwards)

rVeh(1,1) = r0*cos(theta0); %x
rVeh(1,2) = r0*sin(theta0); %y
if (v0 < 0)
    psiAbs = theta0 - pi/2 + psi0; 
else
    psiAbs = theta0 + pi/2 - psi0;
end
v0Abs = abs(v0);
vVeh(1,1) = v0Abs*cos(psiAbs); %Vx
vVeh(1,2) = v0Abs*sin(psiAbs); %Vy [m/sec]

dt = dt_min;
T_sim = 0;
i = 1;
%number of iterations is not predetermined
    
while (T_sim < T_sim_max),
    
    MeanAnom = 2*pi*(T_sim+MnTinit)/T_mn;
    %true anomaly as calculated by Kepler:
    phi = TrueAnom(MeanAnom,e_mn,4);
    r_mn_polar(i) = lts_mn/(1+e_mn*cos(phi));
    r_mn(i,1) = r_mn_polar(i)*cos(phi);
    r_mn(i,2) = r_mn_polar(i)*sin(phi); 
    
    r_E_veh(i,1) = rVeh(i,1); %Earth to vehicle distance along x 
    r_E_veh(i,2) = rVeh(i,2); %Earth to vehicle distance along y
    r_mn_veh(i,1) = rVeh(i,1) - r_mn(i,1); %Moon to vehicle distance along x
    r_mn_veh(i,2) = rVeh(i,2) - r_mn(i,2); %Moon to vehicle distance along y
    
    r_E_veh2 = (r_E_veh(i,1)^2+r_E_veh(i,2)^2)^(1/2); %Earth to vehicle 2D distance
    r_mn_veh2 = (r_mn_veh(i,1)^2+r_mn_veh(i,2)^2)^(1/2); %Moon to vehicle 2D distance
       
    r_E_veh3 = r_E_veh2^3;
    r_mn_veh3 = r_mn_veh2^3;
    
    ax(i) = - u_mn*r_mn_veh(i,1)/r_mn_veh3 - u_E*r_E_veh(i,1)/r_E_veh3;
    ay(i) = - u_mn*r_mn_veh(i,2)/r_mn_veh3 - u_E*r_E_veh(i,2)/r_E_veh3;
    %calculate the acceleration and speed, then make an aprox of the next position
    
    vVeh(i+1,1) = vVeh(i,1) + ax(i)*dt;
    vVeh(i+1,2) = vVeh(i,2) + ay(i)*dt;
    
    rVeh(i+1,1) = rVeh(i,1) + (vVeh(i,1) + vVeh(i+1,1))*dt/2; %x
    rVeh(i+1,2) = rVeh(i,2) + (vVeh(i,2) + vVeh(i+1,2))*dt/2; %y

    if (r_E_veh2 < Re99)
        outStatus = 'A';
        break;    
    end
    
    if (r_mn_veh2 < Rm99)
        outStatus = 'C';
        break;    
    end

    %k = |a x V|/|V^3| - curvature of a 2 D line:
    k = 10^7*abs(-ax(i)*vVeh(i,2)+ay(i)*vVeh(i,1))/[(vVeh(i,1)^2+vVeh(i,2)^2)^(3/2)];
    %time increment dt used is inversely proportional to trajectory's
    %curvature k
    
    %maximum gravitational acceleration experienced in zone B has an order
    %of a = 0.31 m/(s^2) just by the zone C boundary
        
    %new dt is calculated each iteration so that bigger
    %increments were used when trajectory's curvature is smaller and vice-versa
    
    if (k > 1/dt_min)
        dt = dt_min;
    else 
        if (k < 0.0083)
            dt = 120; %sec
        else
            dt = 1/k;
        end
    end
    
    T_sim = T_sim + dt;
    i = i+1;
end

Tflight = T_sim;

if (T_sim >= T_sim_max) outStatus = 'timeout'; end

%-----------------------Output values-----------------------------
rOut = (rVeh(i,1)^2 + rVeh(i,2)^2)^0.5;
thetaOut = atan2(rVeh(i,2),rVeh(i,1)); %-pi <= ATAN2(Y,X) <= pi

vOut = (vVeh(i,1)^2 + vVeh(i,2)^2)^0.5;
psiOutAbs = atan2(vVeh(i,2),vVeh(i,1)); %-pi <= ATAN2(Y,X) <= pi
if (psiOutAbs < 0) psiOutAbs = psiOutAbs + 2*pi; end

if (v0 < 0) 
    vOut = -vOut;
    rNormal = thetaOut - pi/2;
    if (rNormal < 0) rNormal = rNormal + 2*pi; end
    psiOut = psiOutAbs - rNormal;
else
    psiOut = thetaOut + pi/2 - psiOutAbs;
end
%-----------------------------------------------------------------

