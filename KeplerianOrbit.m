function [lts, e, periapsisArg, retro, iMax, tFlIter, varargout1, outStatus] = KeplerianOrbit(r0, theta0, v0, psi0, zone, T_sim)

           
%r0 = 4*10^6;
%theta0 = 0; 
%v0 = 3000;
%psi0 = -pi/4;
%zone = 'C';
            
%------------physical constants & boundaries of the system-----------------
%--------------------------------------------------------------------------

G = 6.67427*(10)^(-11);
ReAtm = 6.6*(10^6);
Re99 = 1.723*(10^8);
Rm99 = 4*(10^6);
RmL = 1.74*(10^6);

if (zone == 'A')
    %for Terrestrial neighborhood:
    M_E = 5.9736*(10^24); %[kg]
    u = G*M_E; %[(m^3)/(kg*(s^2)); u is the gravitational parameter 
    RlimitH = Re99;
    RlimitL = ReAtm;
elseif (zone == 'C')
    %for Lunar neighborhood
    M_M = 7.3477*(10^22); %[kg]
    u = G*M_M;
    RlimitH = Rm99;
    RlimitL = RmL;
end

% R = (r,thetha) (position is described in angular coordinate system with 
% the origin in the centre of the Main Gravitational Body)
% velocity V = (v, xsi) (xsi - an angle btw the V vector and a circle of 
% radius R centred in the origin; -90 to 90; positive outwards)

h = abs(v0*r0*cos(psi0));

a = u*r0/(-(v0^2)*r0 + 2*u); %from vis viva equation

e = sqrt(1-(h^2)/(u*a)); %eccentricity

lts = (1-e^2)*a; %semi-latus rectum

%present true anomaly of the orbiting object is given by:
TrueAnom0 = real(acos( (1/e)*((lts/r0) - 1))); 
%the outcome (in radians) is between 0 and pi (true value is in the range of -pi to pi)

if (psi0 < 0) %object on the "inbound" leg
    %the object TrueAnomaly from it's trayectory's periapsis is btw -pi and 0
    TrueAnom0 = -TrueAnom0;
end


if (v0 > 0) 
    %trajectory angular speed vector points counterclockwise
    %TrueAnomaly orbital coordinate increases counterclockwise
    periapsisArg = theta0 - TrueAnom0;
    retro = 0; %direct orbit
elseif (v0 < 0) 
    %trajectory angular speed vector points clockwise
    %TrueAnomaly orbital coordinate increases clockwise
    periapsisArg = theta0 + TrueAnom0;
    retro = 1; %retrograde orbit
else
    outStatus = 'V'; %vertical trajectory    
end

%varargout1 = [lts e periapsisArg retro 0 0 0 0 0];
varargout1 = [0 0 0 0 0 0 0 0];

if (periapsisArg < 0) periapsisArg = periapsisArg + 2*pi; end
if (periapsisArg > 2*pi) periapsisArg = periapsisArg - 2*pi; end

r_periapsis = a*(1-e);
v_periapsis = h/r_periapsis;

%--------------------------------------------------------------------------
%----------------------------selection algorithm:--------------------------
option = 0;
if (e < 1) 
    option = 8;
    r_apoapsis = a*(1+e);
    v_apoapsis = h/r_apoapsis;
else
    r_apoapsis = 1e20; %arbitrary big number - in reality it's infinity
end
if (r_periapsis > RlimitL) option = option + 4; end
if (r_apoapsis > RlimitH) option = option + 2; end
if (TrueAnom0 > 0) option = option + 1; end
 

switch option
    case {11,14,15} 
        %-----------------------alpha scenario-----------------------------
        %leaving Planet's region via elliptic trajectory
        outStatus = 'B';
        
        TrueAnomOut = acos( (1/e)*( ((lts/RlimitH) - 1) )); %allways > 0 on the outbound leg
        if (retro == 0) 
            thetaOut = periapsisArg + TrueAnomOut;
        else 
            thetaOut = periapsisArg - TrueAnomOut;
        end
        if (thetaOut > 2*pi) thetaOut = thetaOut - 2*pi; end
        if (thetaOut < 0) thetaOut = thetaOut + 2*pi; end
        
        vOut = sqrt(u*(2/RlimitH-1/a)); %from vis-viva
        if (retro) vOut = -vOut; end
        psiOut = acos(abs(h/(vOut*RlimitH))); %should be smaller than pi/2 but positive (outbound)
        
        %-----------------------time of flight-----------------------------
        T_orbit = 2*pi*sqrt((a^3)/u);
        EccAnom0 = atan2(sqrt(1-e^2)*sin(TrueAnom0),e+cos(TrueAnom0)); %EccAnom = arg(Y,X)
        MeanAnom0 = EccAnom0 - e*sin(EccAnom0);
        %time from last theoretical periapsis:
        tLstPeriapsis = T_orbit*MeanAnom0/(2*pi);
        %when tLstPeriapsis is < 0 option 14nth is realized
        
        EccAnomOut = atan2(sqrt(1-e^2)*sin(TrueAnomOut),e+cos(TrueAnomOut)); %EccAnom = arg(Y,X)
        MeanAnomOut = EccAnomOut - e*sin(EccAnomOut);
        tOut = T_orbit*MeanAnomOut/(2*pi);
        
        tFlight = tOut - tLstPeriapsis;
              
        %-----------------------plotting limits----------------------------
        d_theta = TrueAnomOut - TrueAnom0;        
        if (d_theta < 0) d_theta = d_theta + 2*pi; end
        if (d_theta > 2*pi) d_theta = d_theta - 2*pi; end

        iMax = ceil(1000*(d_theta)/(2*pi)); %limits the plotting loop
        %-------------------------output data------------------------------
        varargout1 = [RlimitH thetaOut vOut psiOut tFlight option r_periapsis v_periapsis]; 
        %--------------------iterated time variable------------------------
        if (retro) 
            iMax = 1000 - iMax;
            for i = iMax:1:1000
                fiVehMn = 2*pi*i/1000 + theta0;
                
                TrueAnomOutIt = periapsisArg - fiVehMn;
                
                EccAnomOut = atan2(sqrt(1-e^2)*sin(TrueAnomOutIt),e+cos(TrueAnomOutIt)); %EccAnom = arg(Y,X)
                MeanAnomOutIt(1001 - i) = EccAnomOut - e*sin(EccAnomOut);
                tOutIt = T_orbit*MeanAnomOutIt(1001 - i)/(2*pi);
                if ((tOutIt < 0) && (tLstPeriapsis > 0)) tOutIt = - tOutIt; end %i.e not if option == 14
                tFlIter(1001 - i) = tOutIt - tLstPeriapsis;
                
                %-------------upgraded output------------------------------
                if ((tFlight > T_sim) && (i~=iMax))
                    if ((tFlIter(1001-(i-1)) > T_sim) && (tFlIter(1001-i) <= T_sim)) 
                        %calculate output at (1001-i)
                        thetaOut = fiVehMn;
                        rOut = lts/(e*cos(TrueAnomOutIt) + 1);
                        vOut = sqrt(u*(2/rOut-1/a));
                        psiOut = acos(abs(h/(vOut*rOut)));
                        if ((option == 14) && (tOutIt < 0)) psiOut = -psiOut; end
                        varargout1 = [rOut thetaOut vOut psiOut tFlIter(1001-i) option r_periapsis v_periapsis r_apoapsis v_apoapsis T_orbit];
                    end
                end
            end
        else
            for i = 1:1:iMax+1
                fiVehMn = 2*pi*(i-1)/1000 + theta0;
                
                TrueAnomOutIt = fiVehMn - periapsisArg;
                
                EccAnomOut = atan2(sqrt(1-e^2)*sin(TrueAnomOutIt),e+cos(TrueAnomOutIt)); %EccAnom = arg(Y,X)
                MeanAnomOut = EccAnomOut - e*sin(EccAnomOut);
                tOutIt = T_orbit*MeanAnomOut/(2*pi);
                if ((tOutIt < 0) && (tLstPeriapsis > 0)) tOutIt = - tOutIt; end %i.e not if option == 14
                tFlIter(i) = tOutIt - tLstPeriapsis;
                
                %-------------upgraded output------------------------------
                if ((tFlight > T_sim) && (i~=1))
                    if ((tFlIter(i-1) <= T_sim) && (tFlIter(i) > T_sim))
                        %calculate output at i
                        thetaOut = fiVehMn;
                        rOut = lts/(e*cos(TrueAnomOutIt) + 1);
                        vOut = sqrt(u*(2/rOut-1/a));
                        psiOut = acos(abs(h/(vOut*rOut)));
                        if ((option == 14) && (tOutIt < 0)) psiOut = -psiOut; end
                        varargout1 = [rOut thetaOut vOut psiOut tFlIter(i) option r_periapsis v_periapsis r_apoapsis v_apoapsis T_orbit];
                    end
                end
            end
        end
        %------------------------------------------------------------------

        
    case {8,9,10} 
        %--------------------------beta scenario---------------------------
        %entering Planet's atmosphere along elliptic trajectory
        outStatus = 'Z'; %periapsis within Atmosphere (or dangerously close to the surface)
        
        TrueAnomLout = -acos( (1/e)*( ((lts/RlimitL) - 1) )); %allways < 0 on the inbound leg
        if (retro == 0) 
            thetaLout = periapsisArg + TrueAnomLout;
        else 
            thetaLout = periapsisArg - TrueAnomLout;
        end
        if (thetaLout > 2*pi) thetaLout = thetaLout - 2*pi; end
        if (thetaLout < 0) thetaLout = thetaLout + 2*pi; end
        
        vLout = sqrt(u*(2/RlimitL-1/a)); %from vis-viva
        if (retro) vLout = -vLout; end
        psiLout = -acos(abs(h/(vLout*RlimitL))); %should be bigger than -pi/2 but negative (inbound leg)       
        
        %-----------------------time of flight-----------------------------
        T_orbit = 2*pi*sqrt((a^3)/u); %theoretical (assuming 0 Planet's radius)
        EccAnom0 = atan2(sqrt(1-e^2)*sin(TrueAnom0),e+cos(TrueAnom0)); %EccAnom = arg(Y,X)
        MeanAnom0 = EccAnom0 - e*sin(EccAnom0);
        %time from last theoretical periapsis:
        tLstPeriapsis = T_orbit*MeanAnom0/(2*pi);
        if (tLstPeriapsis < 0) tLstPeriapsis = T_orbit + tLstPeriapsis; end

        EccAnomLout = atan2(sqrt(1-e^2)*sin(TrueAnomLout),e+cos(TrueAnomLout)); %EccAnom = arg(Y,X)
        MeanAnomLout = EccAnomLout - e*sin(EccAnomLout);
        %time from theoretical periapsis to atmosphere entry (or surface impact)
        tLout = T_orbit*MeanAnomLout/(2*pi);
        if (tLout < 0) tLout = T_orbit + tLout; end
        
        tFlight = tLout - tLstPeriapsis;
              
        %-----------------------plotting limits----------------------------
        d_theta = TrueAnomLout - TrueAnom0;
        if (d_theta < 0) d_theta = d_theta + 2*pi; end
        if (d_theta > 2*pi) d_theta = d_theta - 2*pi; end

        iMax = ceil(1000*(d_theta)/(2*pi)); %limits the plotting loop
        
        %-------------------------output data------------------------------
        varargout1 = [RlimitL thetaLout vLout psiLout tFlight option r_periapsis v_periapsis]; 
        %--------------------iterated time variable------------------------
        if (retro) 
            iMax = 1000 - iMax;
            for i = iMax:1:1000
                fiVehMn = 2*pi*i/1000 + theta0;
                
                TrueAnomOutIt = periapsisArg - fiVehMn;
                
                EccAnomOut = atan2(sqrt(1-e^2)*sin(TrueAnomOutIt),e+cos(TrueAnomOutIt)); %EccAnom = arg(Y,X)
                MeanAnomOut = EccAnomOut - e*sin(EccAnomOut);
                tOutIt = T_orbit*MeanAnomOut/(2*pi);
                if ((tOutIt < 0) && (tLstPeriapsis > 0)) tOutIt = T_orbit + tOutIt; end %i.e not if option == 9
                tFlIter(1001 - i) = tOutIt - tLstPeriapsis;
                
                %-------------upgraded output------------------------------
                if ((tFlight > T_sim) && (i~=iMax))
                    if ((tFlIter(1001-(i-1)) > T_sim) && (tFlIter(1001-i) <= T_sim)) 
                        %calculate output at (1001-i)
                        thetaLout = fiVehMn;
                        rLout = lts/(e*cos(TrueAnomOutIt) + 1);
                        vLout = sqrt(u*(2/rLout-1/a));
                        psiLout = -acos(abs(h/(vLout*rLout)));
                        if ((option == 9) && (tFlIter(1001-i) < (T_orbit/2 - tLstPeriapsis))) psiLout = -psiLout; end
                        varargout1 = [rLout thetaLout vLout psiLout tFlIter(1001-i) option r_periapsis v_periapsis r_apoapsis v_apoapsis T_orbit];
                    end
                end
            end
        else
            for i = 1:1:iMax+1
                fiVehMn = 2*pi*(i-1)/1000 + theta0;
                
                TrueAnomOutIt = fiVehMn - periapsisArg;
                
                EccAnomOut = atan2(sqrt(1-e^2)*sin(TrueAnomOutIt),e+cos(TrueAnomOutIt)); %EccAnom = arg(Y,X)
                MeanAnomOut = EccAnomOut - e*sin(EccAnomOut);
                tOutIt = T_orbit*MeanAnomOut/(2*pi);
                if ((tOutIt < 0) && (tLstPeriapsis > 0)) tOutIt = T_orbit + tOutIt; end %i.e not if option == 9
                tFlIter(i) = tOutIt - tLstPeriapsis;
                
                %-------------upgraded output------------------------------
                if ((tFlight > T_sim) && (i~=1))
                    if ((tFlIter(i-1) <= T_sim) && (tFlIter(i) > T_sim))
                        %calculate output at i
                        thetaLout = fiVehMn;
                        rLout = lts/(e*cos(TrueAnomOutIt) + 1);
                        vLout = sqrt(u*(2/rLout-1/a));
                        psiLout = -acos(abs(h/(vLout*rLout)));
                        if ((option == 9) && (tFlIter(i) < (T_orbit/2 - tLstPeriapsis))) psiLout = -psiLout; end
                        varargout1 = [rLout thetaLout vLout psiLout tFlIter(i) option r_periapsis v_periapsis r_apoapsis v_apoapsis T_orbit];
                    end
                end
            end
        end
        %------------------------------------------------------------------ 
        
    case {12,13} 
        %--------------------------gamma scenario--------------------------
        %maintaining an orbit around Earth or Moon
        outStatus = 'O';
        
        T_orbit = 2*pi*sqrt((a^3)/u);
        EccAnom0 = atan2(sqrt(1-e^2)*sin(TrueAnom0),e+cos(TrueAnom0)); %EccAnom = arg(Y,X)
        MeanAnom0 = EccAnom0 - e*sin(EccAnom0);
        %time from last theoretical periapsis:
        tLstPeriapsis = T_orbit*MeanAnom0/(2*pi);
        if (tLstPeriapsis < 0) tLstPeriapsis = T_orbit + tLstPeriapsis; end
        
        T_simComplement = T_sim;
        k = 0;
        while (T_simComplement > T_orbit)
            T_simComplement = T_simComplement - T_orbit;
            k = k+1;
        end
        
        %--------------------iterated time variable------------------------
        if (retro) 
            iMax = 1;
            for i = 1:1:1000
                fiVehMn = 2*pi*(i-1)/1000 + theta0;
                
                TrueAnomOutIt = periapsisArg - fiVehMn;
                
                EccAnomOut = atan2(sqrt(1-e^2)*sin(TrueAnomOutIt),e+cos(TrueAnomOutIt)); %EccAnom = arg(Y,X)
                MeanAnomOut = EccAnomOut - e*sin(EccAnomOut);
                tOutIt = T_orbit*MeanAnomOut/(2*pi);
                if (tOutIt < 0) tOutIt = T_orbit + tOutIt; end
                tFlIter(1001 - i) = tOutIt - tLstPeriapsis;
                if (tFlIter(1001 - i) < 0) tFlIter(1001 - i) = T_orbit + tFlIter(1001 - i); end
                
                %-------------upgraded output------------------------------
                if ((i~=1) && (tFlIter(1001-(i-1)) > T_simComplement) && (tFlIter(1001-i) <= T_simComplement)) 
                    %calculate output at (1001-i)
                    thetaOut = fiVehMn;
                    rOut = lts/(e*cos(TrueAnomOutIt) + 1);
                    vOut = sqrt(u*(2/rOut-1/a));
                    psiOut = acos(abs(h/(vOut*rOut)));
                    if ((option == 14) && (tOutIt < 0)) psiOut = -psiOut; end
                    %-------------------------output data--------------
                    varargout1 = [rOut thetaOut vOut psiOut tFlIter(1001-i) option r_periapsis v_periapsis r_apoapsis v_apoapsis T_orbit k T_simComplement];
                end
            end
            tFlIter(1000) = tFlIter(999);
        else
            iMax = 1000;
            for i = 1:1:1000
                fiVehMn = 2*pi*i/1000 + theta0;
                
                TrueAnomOutIt = fiVehMn - periapsisArg;
                
                EccAnomOut = atan2(sqrt(1-e^2)*sin(TrueAnomOutIt),e+cos(TrueAnomOutIt)); %EccAnom = arg(Y,X)
                MeanAnomOut = EccAnomOut - e*sin(EccAnomOut);
                tOutIt = T_orbit*MeanAnomOut/(2*pi);
                if (tOutIt < 0) tOutIt = T_orbit + tOutIt; end
                tFlIter(i) = tOutIt - tLstPeriapsis;
                if (tFlIter(i) < 0) tFlIter(i) = T_orbit + tFlIter(i); end
                
                %-------------upgraded output------------------------------
                if ((i~=1) && (tFlIter(i-1) <= T_simComplement) && (tFlIter(i) > T_simComplement))
                    %calculate output at i
                    thetaOut = fiVehMn;
                    rOut = lts/(e*cos(TrueAnomOutIt) + 1);
                    vOut = sqrt(u*(2/rOut-1/a));
                    psiOut = acos(abs(h/(vOut*rOut)));
                    if ((option == 14) && (tOutIt < 0)) psiOut = -psiOut; end
                    %-------------------------output data--------------
                    varargout1 = [rOut thetaOut vOut psiOut tFlIter(i) option r_periapsis v_periapsis r_apoapsis v_apoapsis T_orbit k T_simComplement];
                end
            end
            tFlIter(1000) = tFlIter(999);
        end
        %------------------------------------------------------------------
        
    case {3,6,7} 
        %--------------------------delta scenario--------------------------
        %leaving Planet's region via hyper-/para-bolic trajectory
        outStatus = 'Bh';
        
        TrueAnomOut = acos( (1/e)*( ((lts/RlimitH) - 1) )); %allways > 0 on the outbound leg
        if (retro == 0) 
            thetaOut = periapsisArg + TrueAnomOut;
        else 
            thetaOut = periapsisArg - TrueAnomOut;
        end
        if (thetaOut > 2*pi) thetaOut = thetaOut - 2*pi; end
        if (thetaOut < 0) thetaOut = thetaOut + 2*pi; end
        
        vOut = sqrt(u*(2/RlimitH-1/a)); %from vis-viva
        if (retro) vOut = -vOut; end
        psiOut = acos(abs(h/(vOut*RlimitH))); %should be smaller than pi/2 but positive (outbound)
        
        %option 6 has a negative initial TrueAnom0 which yields a negative
        %value of Ainbound
        
        %calculating time of flight using constant area speed
        Ainbound = (lts^2)/(2*sqrt(e^2-1))*(e*sin(TrueAnom0)/(e*cos(TrueAnom0)+1) - 2*atanh((e-1)*tan(TrueAnom0/2)/sqrt(e^2-1)));        
        Aoutbound = (lts^2)/(2*sqrt(e^2-1))*(e*sin(TrueAnomOut)/(e*cos(TrueAnomOut)+1) - 2*atanh((e-1)*tan(TrueAnomOut/2)/sqrt(e^2-1)));
        
        Atot = Aoutbound - Ainbound;
        tFlight = abs(2*Atot/h);
        
        %-----------------------plotting limits----------------------------
        d_theta = TrueAnomOut - TrueAnom0;
        if (d_theta < 0) d_theta = d_theta + 2*pi; end
        if (d_theta > 2*pi) d_theta = d_theta - 2*pi; end

        iMax = ceil(1000*(d_theta)/(2*pi)); %limits the plotting loop
        
        %-------------------------output data------------------------------
        varargout1 = [RlimitH thetaOut vOut psiOut tFlight option r_periapsis v_periapsis]; 
        %--------------------iterated time variable------------------------
        if (retro) 
            iMax = 1000 - iMax;
            for i = iMax:1:1000
                fiVehMn = 2*pi*i/1000 + theta0;
                
                TrueAnomOutIt = periapsisArg - fiVehMn;
                
                AoutboundIt = (lts^2)/(2*sqrt(e^2-1))*(e*sin(TrueAnomOutIt)/(e*cos(TrueAnomOutIt)+1) - 2*atanh((e-1)*tan(TrueAnomOutIt/2)/sqrt(e^2-1)));
                AtotIt = AoutboundIt - Ainbound;
                tFlIter(1001-i) = abs(2*AtotIt/h);
                %-------------upgraded output------------------------------
                if ((tFlight > T_sim) && (i~=iMax))
                    if ((tFlIter(1001-(i-1)) > T_sim) && (tFlIter(1001-i) <= T_sim)) 
                        %calculate output at (1001-i)
                        thetaOut = fiVehMn;
                        rOut = lts/(e*cos(TrueAnomOutIt) + 1);
                        vOut = sqrt(u*(2/rOut-1/a));
                        psiOut = acos(abs(h/(vOut*rOut)));
                        if ((option == 6) && (AoutboundIt < 0)) psiOut = -psiOut; end
                        varargout1 = [rOut thetaOut vOut psiOut tFlIter(1001-i) option r_periapsis v_periapsis];
                    end
                end
            end
        else
            for i = 1:1:iMax+1
                fiVehMn = 2*pi*(i-1)/1000 + theta0;
                
                TrueAnomOutIt = fiVehMn - periapsisArg;
                
                AoutboundIt(i) = (lts^2)/(2*sqrt(e^2-1))*(e*sin(TrueAnomOutIt)/(e*cos(TrueAnomOutIt)+1) - 2*atanh((e-1)*tan(TrueAnomOutIt/2)/sqrt(e^2-1)));
                AtotIt(i) = AoutboundIt(i) - Ainbound;
                tFlIter(i) = abs(2*AtotIt(i)/h);
                
                %-------------upgraded output------------------------------
                if ((tFlight > T_sim) && (i~=1))
                    if ((tFlIter(i-1) <= T_sim) && (tFlIter(i) > T_sim))
                        %calculate output at i
                        thetaOut = fiVehMn;
                        rOut = lts/(e*cos(TrueAnomOutIt) + 1);
                        vOut = sqrt(u*(2/rOut-1/a));
                        psiOut = acos(abs(h/(vOut*rOut)));
                        if ((option == 14) && (AoutboundIt < 0)) psiOut = -psiOut; end
                        varargout1 = [rOut thetaOut vOut psiOut tFlIter(i) option r_periapsis v_periapsis];
                    end
                end
            end
        end
        %------------------------------------------------------------------
    
    case {2} 
        %---------------------------zeta scenario--------------------------
        %entering Planet's atmosphere hyper-/para-bolic elliptic trajectory
        outStatus = 'Zh';
        
        TrueAnomLout = -acos( (1/e)*( ((lts/RlimitL) - 1) )); %allways < 0 on the inbound leg
        if (retro == 0) 
            thetaLout = periapsisArg + TrueAnomLout;
        else 
            thetaLout = periapsisArg - TrueAnomLout;
        end
        if (thetaLout > 2*pi) thetaLout = thetaLout - 2*pi; end
        if (thetaLout < 0) thetaLout = thetaLout + 2*pi; end
        
        vLout = sqrt(u*(2/RlimitL-1/a)); %from vis-viva
        if (retro) vLout = -vLout; end
        psiLout = -acos(abs(h/(vLout*RlimitL))); %should be bigger than -pi/2 but negative (entry)
        
        %calculating time of flight using constant area speed
        Ainbound = (lts^2)/(2*sqrt(e^2-1))*(e*sin(TrueAnom0)/(e*cos(TrueAnom0)+1) - 2*atanh((e-1)*tan(TrueAnom0/2)/sqrt(e^2-1)));        
        Aoutbound = (lts^2)/(2*sqrt(e^2-1))*(e*sin(TrueAnomLout)/(e*cos(TrueAnomLout)+1) - 2*atanh((e-1)*tan(TrueAnomLout/2)/sqrt(e^2-1)));
        
        Atot = Aoutbound - Ainbound;
        tFlight = abs(2*Atot/h);
        
        %-----------------------plotting limits----------------------------
        d_theta = TrueAnomLout - TrueAnom0;        
        if (d_theta < 0) d_theta = d_theta + 2*pi; end
        if (d_theta > 2*pi) d_theta = d_theta - 2*pi; end
        
        iMax = ceil(1000*(d_theta)/(2*pi)); %limits the plotting loop 
        
        %-------------------------output data------------------------------
        varargout1 = [RlimitL thetaLout vLout psiLout tFlight option r_periapsis v_periapsis]; 
        %--------------------iterated time variable------------------------
        if (retro) 
            iMax = 1000 - iMax;
            for i = iMax:1:1000
                fiVehMn = 2*pi*i/1000 + theta0;
                
                TrueAnomOutIt = periapsisArg - fiVehMn;
                
                AoutboundIt = (lts^2)/(2*sqrt(e^2-1))*(e*sin(TrueAnomOutIt)/(e*cos(TrueAnomOutIt)+1) - 2*atanh((e-1)*tan(TrueAnomOutIt/2)/sqrt(e^2-1)));
                AtotIt = AoutboundIt - Ainbound;
                tFlIter(1001-i) = abs(2*AtotIt/h);
                
                %-------------upgraded output------------------------------
                if ((tFlight > T_sim) && (i~=iMax))
                    if ((tFlIter(1001-(i-1)) > T_sim) && (tFlIter(1001-i) <= T_sim)) 
                        %calculate output at (1001-i)
                        thetaLout = fiVehMn;
                        rLout = lts/(e*cos(TrueAnomOutIt) + 1);
                        vLout = sqrt(u*(2/rLout-1/a));
                        psiLout = -acos(abs(h/(vLout*rLout)));
                        varargout1 = [rLout thetaLout vLout psiLout tFlIter(1001-i) option r_periapsis v_periapsis];
                    end
                end
            end
        else
            for i = 1:1:iMax+1
                fiVehMn = 2*pi*(i-1)/1000 + theta0;
                
                TrueAnomOutIt = fiVehMn - periapsisArg;
                
                AoutboundIt = (lts^2)/(2*sqrt(e^2-1))*(e*sin(TrueAnomOutIt)/(e*cos(TrueAnomOutIt)+1) - 2*atanh((e-1)*tan(TrueAnomOutIt/2)/sqrt(e^2-1)));
                AtotIt = AoutboundIt - Ainbound;
                tFlIter(i) = abs(2*AtotIt/h);
                
                %-------------upgraded output------------------------------
                if ((tFlight > T_sim) && (i~=1))
                    if ((tFlIter(i-1) <= T_sim) && (tFlIter(i) > T_sim))
                        %calculate output at i
                        thetaLout = fiVehMn;
                        rLout = lts/(e*cos(TrueAnomOutIt) + 1);
                        vLout = sqrt(u*(2/rLout-1/a));
                        psiLout = -acos(abs(h/(vLout*rLout)));
                        varargout1 = [rLout thetaLout vLout psiLout tFlIter(i) option r_periapsis v_periapsis];
                    end
                end
            end
        end
        %------------------------------------------------------------------
end

%--------------------------------------------------------------------------
%--------------------------Visualisation-----------------------------------

if (~retro)
    for i=0:1:iMax %correct for normal trajectories
        fi = 2*pi*i/1000 + theta0; %starting point
        trayectory(i+1,1) = lts*cos(fi)/(e*cos(fi - periapsisArg) + 1);
        trayectory(i+1,2) = lts*sin(fi)/(e*cos(fi - periapsisArg) + 1);
    end
else
    for i=iMax:1:1000 %correct for retrograde trajectories
        fi = 2*pi*i/1000 + theta0; %starting point
        trayectory(1001-i,1) = lts*cos(fi)/(e*cos(fi - periapsisArg) + 1);
        trayectory(1001-i,2) = lts*sin(fi)/(e*cos(fi - periapsisArg) + 1);
    end
end

R_E = 6370000;
R_mn = 1730000;

for i=0:1:200 %to show earth's and moon's surfaces
    E(i+1,1)= R_E*cos(2*pi*i/200);
    E(i+1,2)= R_E*sin(2*pi*i/200);
    M(i+1,1)= R_mn*cos(2*pi*i/200);
    M(i+1,2)= R_mn*sin(2*pi*i/200);    
end

if (zone == 'A')
    figure(3);
    plot(trayectory(:,1), trayectory(:,2), 'b', E(:,1), E(:,2), 'g');
    title('Part of the flight in the neighborhood of Earth');
    legend('Covered path(Object)','Earth',2);
    axis equal;
    xlabel('metres');
    ylabel('metres');
end
