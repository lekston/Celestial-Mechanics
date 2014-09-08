function [rVeh2nd, thetaVeh2nd, vVeh2nd, psiVeh2nd] = AngCoordConv(MnTpres, conv, rVeh1st, thetaVeh1st, vVeh1st, psiVeh1st)

%rVeh1st = 3.615*10^8;
%thetaVeh1st =0;
%vVeh1st = 1076;
%psiVeh1st = -pi/2;
%MnTpres = 0;
%conv = 'EM';

G = 6.67427*(10)^(-11);
M_E = 5.9736*(10^24); %for Moon movement calculations
u = G*M_E;

%Moon's orbit characteristics:
T_mn = 2.361*(10)^6; %[sec]
e_mn = 0.0549;
a_mn = 384.4*(10)^6;
lts_mn = a_mn*(1-e_mn^2);
h_mn = sqrt(lts_mn*u);

% conv = 'EM' | 'ME'

%present moon position:
MeanAnomMn = 2*pi*(MnTpres)/T_mn;
thetaMn = TrueAnom(MeanAnomMn,e_mn,4);
rMn = lts_mn/(1+e_mn*cos(thetaMn));
vMn = sqrt(u*(2/rMn-1/a_mn));
psiMn = real(acos(h_mn/(vMn*rMn)));
if (thetaMn > pi) psiMn = -psiMn; end;
psiAbsMn = thetaMn + pi/2 - psiMn; %as moon's orbit is direct
        
%-----------Angular coordinates in the other system-----
%position of the moving object must be expressed in coordinates of
%a angular system related to the moon i.e: rVeh2nd, thetaVeh2nd, vVeh2nd, psiVeh2nd

if (conv == 'EM')
    alpha = thetaVeh1st - thetaMn;
    vMn = -vMn; %since we will be substracting its value
else
    thetaE = thetaMn - pi;
    if (thetaE < 0) thetaE = thetaE + 2*pi; end
    alpha = thetaVeh1st - thetaE;
end

rVeh2nd = sqrt(rMn^2 + rVeh1st^2 - 2*rMn*rVeh1st*cos(alpha)); %law of cosines
beta = asin(rVeh1st*sin(alpha)/rVeh2nd);%law of sines, ambiguous outcome
if (conv == 'EM')
    if (rVeh1st^2 > (rMn^2 + rVeh2nd^2))
        if (beta<0)
            beta = -pi - beta;
        else
            beta = pi - beta;
        end
    end
    
    thetaVeh2nd = pi + thetaMn - beta;
else
    thetaVeh2nd = pi + thetaE - beta;
end
        
if (thetaVeh2nd < 0) thetaVeh2nd = thetaVeh2nd + 2*pi; end
if (thetaVeh2nd > 2*pi) thetaVeh2nd = thetaVeh2nd - 2*pi; end
        
%absolute orientation of velocity vector in earths's coordinates: 
if (vVeh1st < 0)
    psiAbsVeh1st = thetaVeh1st - pi/2 + psiVeh1st; 
else
    psiAbsVeh1st = thetaVeh1st + pi/2 - psiVeh1st;
end

        
%absolute orientation and value of velocity vector in 2nd body's coordinates:        
vVeh1stAbs = abs(vVeh1st);
vVeh2ndX = vVeh1stAbs*cos(psiAbsVeh1st) + vMn*cos(psiAbsMn);
vVeh2ndY = vVeh1stAbs*sin(psiAbsVeh1st) + vMn*sin(psiAbsMn);
vVeh2nd = sqrt(vVeh2ndX^2 + vVeh2ndY^2);
psiAbsVeh2nd = atan2(vVeh2ndY,vVeh2ndX);
if (psiAbsVeh2nd < 0) psiAbsVeh2nd = psiAbsVeh2nd + 2*pi; end
        
psiTmp = thetaVeh2nd + pi/2 - psiAbsVeh2nd;
if (psiTmp >= pi) psiTmp = psiTmp - 2*pi; end
if (psiTmp < -pi) psiTmp = psiTmp + 2*pi; end
        
if (abs(psiTmp)<(pi/2))
    psiVeh2nd = psiTmp;
else
    vVeh2nd = -vVeh2nd;
    rNormal = thetaVeh2nd - pi/2;
    if (rNormal < 0) rNormal = rNormal + 2*pi; end
    psiVeh2nd = psiAbsVeh2nd - rNormal;
end
        
%--------------------------End of Function---------------------------------
%--------------------------------------------------------------------------