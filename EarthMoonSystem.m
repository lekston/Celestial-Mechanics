% Position is described in angular coordinate system with 
% the origin in the centre of the Earth: R = (r,thetha)
% Verocity V = (v, psi) is described by:
% * v - value in m/s,
%   SIGN of v determines if objects trajectory will follow:
%   -> DIRECT orbit (counterclockwise), when v is POSITIVE;
%   -> RETROGRADE orbit (clockwise), when v is NEGATIVE;
% * psi - an angle between the V vector and a circle of radius R centred in the origin
%   psi ranges from -90 (object in a nose down attitude) to 90 (object in a nose up attitude)


%-------ENTER THE INITIAL CONDITION OF THE SIMULATION HERE-------%

%DISTANCE FROM EARTH'S CENTRE
r0              = 10^7;
%ANGULAR COORDINATE (0,2*pi)
theta0          = (12.5/180)*pi;
%VELOCITY IN m/s
v0              = 8850;
%PITCH ANGLE (-pi,pi)
psi0            = 0;
%MINIMUM TIME INCREAMENT (0s,120s)
dt_min          = 5;
%PERIOD OF SIMULATION (0s, inf)
T_sim_max       = 400000;
%MOON POSITION AT THE START OF SIMULATION (0,2*pi)
MnPos0          = 2.5393;


%---------ONLY THE VARIABLES ABOVE ARE TO BE EDITED--------------%
%                                                                %
%                                                                %
%                                                                %
%----ALL OF THE FOLLOWING DATA ARE EITHER PHISICAL CONSTANTS-----%
%---OR VARIABLES REQUIRED BY THE OPERATIONS OF THIS ALGORITHM----%
%------------ANY CHANGES MAY CAUSE PROGRAM FAILURE---------------%

%The presented below EarthMoonSystem is devided into following zones:
%A - Earth's neighborhood - where Earth's graviation comprises at least 99% of all considered forces 
    %Zone A has a radius of:
    Re99 = 1.723*(10^8); %[m]
%B - both Earth's and Moon's influences are taken into account
    %Zone B is simply outside of zone A (except the small circle defining zone C)
%C - Moon's neighborhood - where Moon's gravitation comprises at least 99% of all considered forces
    %Zone C surrounds the Moon and has a radius of:
    Rm99 = 4*(10^6); %[m]

%Boundaries limiting the simulation:
% -> Object is assumed to be subjected to significant atmospheric drag 
    %(which is not included in this model) when its distance to the 
    %centre of the Earth (origin) is less than:
    ReAtm = 6.6*(10^6); %[m]   
% -> Object is assumed to crash on the surface of the Moon when flying
    %lower that 10km over its surface:
    RmnC = 1.74*(10^6); %[m]

%Simulation stops if:
% -> Either of above mentioned boundaries is reached
% -> Stated maximum simulation time T_sim_max is reached

% This script makes calls to all files stated below:
% * InBetween.m
% * KeplerianOrbit.m
% * TrueAnom.m
% * AngCoordConv.m

close all;

G = 6.67427*(10)^(-11);
M_M = 7.3477*(10^22); %[kg]
u = G*M_M;

R_E = 6370000; %[m]

R_mn = 1730000; %[m]
%Moon's orbit characteristics:
T_mn = 2.361*(10)^6; %[sec]
e_mn = 0.0549;
a_mn = 384.4*(10)^6;
lts_mn = a_mn*(1-e_mn^2);

%Calculating time of moon's flight corresponding to its initial position:
EccAnomMn0 = atan2(sqrt(1-e_mn^2)*sin(MnPos0),e_mn+cos(MnPos0)); %EccAnom = arg(Y,X)
MeanAnomMn0 = EccAnomMn0 - e_mn*sin(EccAnomMn0);
%time from last theoretical periapsis:
MnTinit = T_mn*MeanAnomMn0/(2*pi);
if (MnTinit < 0) MnTinit = T_mn + MnTinit; end

Bvisited = 0;
Cvisited = 0;
Orb = 0;
pValid = 0;
pValidMn = 0;
T_sim = 0;
MnTpres = MnTinit;
zone = 'B';
if (r0 < Re99) zone = 'A'; end

%all data about the trajectory will be added to the trajectory(:,1)
% and trajectory(:,2) variables
iVeh0 = 0;
history = 'H: ';

while (T_sim < T_sim_max),
    
    if (zone == '/')
        % This part of code is entered when the object passes from one of
        % the zones (A,B,C) to another
        zone = outStatus;
        
        %------------------Standarization of variables---------------------
        while (abs(thetaOut) - pi > pi),
            if (thetaOut > 2*pi) thetaOut = thetaOut - 2*pi; end
            if (thetaOut < 0) thetaOut = thetaOut + 2*pi; end
        end
        
        while (abs(psiOut) > (pi/2)),
            if (psiOut > pi) 
                psiOut = psiOut - 2*pi; 
            elseif (psiOut > pi/2)
                psiOut = pi - psiOut;
                vOut = -vOut;
            end
            if (psiOut < -pi) 
                psiOut = psiOut + 2*pi; 
            elseif (psiOut < -pi/2)
                psiOut = -(psiOut + pi);
                vOut = -vOut;
            end
        end
        %------------------------------------------------------------------
        
        r0 = rOut;
        theta0 = thetaOut;
        v0 = vOut;
        psi0 = psiOut;
    end
    
    if  (zone == 'A')
        
        % Zone A, below code corresponds to part of the flight conducted in 
        % the neighborhood of the Earth
        
        T_simAhead = T_sim_max - T_sim;
        % Calculating flight trajectory characteristics according to
        % clasical celestial mechanics 
        % (algorith mcontained in KeplerianOrbit):
        [ltsOrbE, eOrbE, periapsisArgOrbE, retroOrbE, iMaxOrbE, tFlIterOrbE, varargoutOrbE, outStatusOrbE] = ...
            KeplerianOrbit(r0, theta0, v0, psi0, zone, T_simAhead);
        
        outStatus = outStatusOrbE;
        r_periapsisOrbE = varargoutOrbE(7);
        v_periapsisOrbE = varargoutOrbE(8);
        
        %------------------------------visualisation-----------------------
            if (outStatus == 'O')
                T_simAhead = varargoutOrbE(13);
            end
            
            clear projected;
            
            if (~retroOrbE)
                for i=1:1:iMaxOrbE
                    fi = 2*pi*i/1000 + theta0;
                    if (tFlIterOrbE(i) < T_simAhead)
                        trayectory(i + iVeh0, 1) = ltsOrbE*cos(fi)/(eOrbE*cos(fi - periapsisArgOrbE) + 1);
                        trayectory(i + iVeh0, 2) = ltsOrbE*sin(fi)/(eOrbE*cos(fi - periapsisArgOrbE) + 1);
                    else
                        projected(i, 1) = ltsOrbE*cos(fi)/(eOrbE*cos(fi - periapsisArgOrbE) + 1);
                        projected(i, 2) = ltsOrbE*sin(fi)/(eOrbE*cos(fi - periapsisArgOrbE) + 1);
                        pValid = 1;
                    end
                end
                iVeh0 = iMaxOrbE + iVeh0;
            else
                for i=iMaxOrbE:1:1000 %correct for retrograde trajectories (those followed clockwise)
                    fi = 2*pi*i/1000 + theta0;
                    if (tFlIterOrbE(1001 - i) < T_simAhead)
                        trayectory(1001 - i + iVeh0, 1) = ltsOrbE*cos(fi)/(eOrbE*cos(fi - periapsisArgOrbE) + 1);
                        trayectory(1001 - i + iVeh0, 2) = ltsOrbE*sin(fi)/(eOrbE*cos(fi - periapsisArgOrbE) + 1);
                    else
                        projected(1001 - i, 1) = ltsOrbE*cos(fi)/(eOrbE*cos(fi - periapsisArgOrbE) + 1);
                        projected(1001 - i, 2) = ltsOrbE*sin(fi)/(eOrbE*cos(fi - periapsisArgOrbE) + 1);
                        pValid = 1;
                    end
                end
                iVeh0 = 1000 - iMaxOrbE + iVeh0;
            end
        %------------------------------------------------------------------
        %------------------------Flight Data-------------------------------
        rOut = varargoutOrbE(1);
        thetaOut = varargoutOrbE(2);
        vOut = varargoutOrbE(3);
        psiOut = varargoutOrbE(4);
        
        %------------------------------------------------------------------
        
        if (outStatus == 'O')
            
            r_apoapsisOrbE = varargoutOrbE(9);
            v_apoapsisOrbE = varargoutOrbE(10);
            T_orbitOrbE = varargoutOrbE(11);
            k = varargoutOrbE(12);
            
            Tflight = T_sim_max;
            T_sim = T_sim_max;
            MnTpres = T_sim + MnTinit;
            
            history =  strcat(history,zone,'=> Orbiting Earth;');
            outStatus = 'Orbiting Earth';
            Orb = 1;
            break;
            
        else
            
            Tflight = varargoutOrbE(5);
            if (~pValid)
                T_sim = T_sim + Tflight;
            else
                T_sim = T_sim_max;
            end
            MnTpres = T_sim + MnTinit;
        
            if (outStatus == 'Z') | (outStatus == 'Zh')
                
                % entring the atmosphere
                history =  strcat(history,zone,'=>',outStatus,':Atmospheric Entry;');
                zone = outStatus;
                outStatus = 'Atmospheric Entry';
                break;
            
            else %'B' or 'Bh'
                
                % leaving the A zone
                history =  strcat(history,zone,'=>',outStatus,';');
                outStatus = 'B'; %changing from 'Bh' to allow for direct assignement zone = outStatus
                zone = '/';
            end
        end 
    end
        
    if (zone == 'B')
        
        % Zone B, below code corresponds to part of the flight conducted in 
        % "free space" where both Earth and Moon affect the flight path
        
        Bvisited = 1;
        T_simAhead = T_sim_max - T_sim;
        
        MnTpres= T_sim + MnTinit;
        [rVeh, rOut, thetaOut, vOut, psiOut, outStatus, Tflight, r_mn] = InBetween(r0, theta0, v0, psi0, dt_min, T_simAhead, MnTpres);
        
        %-----Updating the trajectory (for visualisation)---- 
        for i = 1:1:size(rVeh,1)
            trayectory(i + iVeh0, 1) = rVeh(i,1);
            trayectory(i + iVeh0, 2) = rVeh(i,2);
        end
        iVeh0 = size(rVeh,1) + iVeh0;
        %----------------------------------------------------
        
        T_sim = T_sim + Tflight;
        MnTpres= T_sim + MnTinit;
        
        history =  strcat(history,zone,'=>',outStatus,';');
        zone = '/';
        
    end

    if (zone == 'C')
        
        % Zone C, below code corresponds to part of the flight conducted in 
        % the neighborhood of the Moon
        
        Cvisited = 1;
        %---------Going to Moon referenced coordinate system--------------
        [rVehMn, thetaVehMn, vVehMn, psiVehMn] = AngCoordConv(MnTpres, 'EM', r0, theta0, v0, psi0);
        %-----------------------------------------------------------------
               
        T_simAhead = T_sim_max - T_sim;
        % Calculating flight trajectory characteristics according to
        % clasical celestial mechanics 
        % (algorith mcontained in KeplerianOrbit):
        [ltsOrbM, eOrbM, periapsisArgOrbM, retroOrbM, iMaxOrbM, tFlIterOrbM, varargoutOrbM, outStatusOrbM] = ...
            KeplerianOrbit(rVehMn, thetaVehMn, vVehMn, psiVehMn, zone, T_simAhead);
        
        outStatus = outStatusOrbM;
        
        if (outStatus == 'O')
            T_orbitMn = varargoutOrbM(11);
            T_simAhead = varargoutOrbM(13); %for visualisation in moon's coordinates
            k = varargoutOrbM(12);
            
            %-----Visualisation in Earth's coordinates (Orbits around the moon)----
            if (~retroOrbM)
                for j=1:1:(k+1)
                    for i=1:1:1000
                        %starting point in Moon's coordinates:
                        fiVehMn = 2*pi*i/1000 + thetaVehMn; 
                        %corresponding moon position:        
                        MeanAnomMn = 2*pi*(MnTpres + tFlIterOrbM(i) + (j-1)*T_orbitMn)/T_mn; 
                        TrueAnomMn = TrueAnom(MeanAnomMn,e_mn,4);
                        rMn = lts_mn/(1+e_mn*cos(TrueAnomMn));
                        %translation vector:
                        r_mn(1) = rMn*cos(TrueAnomMn);
                        r_mn(2) = rMn*sin(TrueAnomMn);
                        
                        if (j < k+1)
                            trayectory(i + iVeh0, 1) = ltsOrbM*cos(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1) + r_mn(1);
                            trayectory(i + iVeh0, 2) = ltsOrbM*sin(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1) + r_mn(2);
                        else
                            if (tFlIterOrbM(i) < T_simAhead)
                                trayectory(i + iVeh0, 1) = ltsOrbM*cos(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1) + r_mn(1);
                                trayectory(i + iVeh0, 2) = ltsOrbM*sin(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1) + r_mn(2);
                                projected(i, 1) = r_mn(1);
                                projected(i, 2) = r_mn(2);
                            else
                                projected(i, 1) = ltsOrbM*cos(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1) + r_mn(1);
                                projected(i, 2) = ltsOrbM*sin(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1) + r_mn(2);
                                pValid = 1;
                            end
                        end
                    
                    end
                    iVeh0 = 1000 + iVeh0;
                end
            else
                for j=1:1:(k+1)
                    for i=1:1:1000 %correct for retrograde trajectories (those followed clockwise)
                        %starting point in Moon's coordinates:
                        fiVehMn = 2*pi*i/1000 + thetaVehMn;
                        %corresponding moon position:
                        MeanAnomMn = 2*pi*(MnTpres + tFlIterOrbM(1001 - i) + (j-1)*T_orbitMn)/T_mn;
                        %this visualisation starts with points at the end of the 
                        %flight path and follows to those at the beginning
                        TrueAnomMn = TrueAnom(MeanAnomMn,e_mn,4);
                        rMn = lts_mn/(1+e_mn*cos(TrueAnomMn));
                        %translation vector:
                        r_mn(1) = rMn*cos(TrueAnomMn);
                        r_mn(2) = rMn*sin(TrueAnomMn);
                        
                        if (j < k+1)
                            trayectory(1001 - i + iVeh0, 1) = ltsOrbM*cos(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1) + r_mn(1);
                            trayectory(1001 - i + iVeh0, 2) = ltsOrbM*sin(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1) + r_mn(2);
                        else
                            if (tFlIterOrbM(1001 - i) < T_simAhead)
                                trayectory(1001 - i + iVeh0, 1) = ltsOrbM*cos(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1) + r_mn(1);
                                trayectory(1001 - i + iVeh0, 2) = ltsOrbM*sin(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1) + r_mn(2);
                                projected(1001 - i, 1) = r_mn(1);
                                projected(1001 - i, 2) = r_mn(2);
                            else
                                projected(1001 - i, 1) = ltsOrbM*cos(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1) + r_mn(1);
                                projected(1001 - i, 2) = ltsOrbM*sin(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1) + r_mn(2);
                                pValid = 1;
                            end
                        end
                    end
                    iVeh0 = 1000 + iVeh0;
                end
            end
            %-------------------------------------------------------------
        else
            %-------Visualisation in Earth's coordinates (others)---------
            clear projected;
            
            if (~retroOrbM)
                for i=1:1:iMaxOrbM+1
                    %starting point in Moon's coordinates:
                    fiVehMn = 2*pi*(i-1)/1000 + thetaVehMn; 
                    %corresponding moon position:        
                    MeanAnomMn = 2*pi*(MnTpres + tFlIterOrbM(i))/T_mn; 
                    TrueAnomMn = TrueAnom(MeanAnomMn,e_mn,4);
                    rMn = lts_mn/(1+e_mn*cos(TrueAnomMn));
                    %translation vector:
                    r_mn(1) = rMn*cos(TrueAnomMn);
                    r_mn(2) = rMn*sin(TrueAnomMn);
                
                    if (tFlIterOrbM(i) < T_simAhead)
                        trayectory(i + iVeh0, 1) = ltsOrbM*cos(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1) + r_mn(1);
                        trayectory(i + iVeh0, 2) = ltsOrbM*sin(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1) + r_mn(2);
                        projected(i, 1) = r_mn(1);
                        projected(i, 2) = r_mn(2);
                    else
                        projected(i, 1) = ltsOrbM*cos(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1) + r_mn(1);
                        projected(i, 2) = ltsOrbM*sin(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1) + r_mn(2);
                        pValid = 1;
                    end
                    
                end
                iVeh0 = iMaxOrbM + iVeh0;
            else
                for i=iMaxOrbM:1:1000 %correct for retrograde trajectories (those followed clockwise)
                    %starting point in Moon's coordinates:
                    fiVehMn = 2*pi*i/1000 + thetaVehMn;
                    %corresponding moon position:
                    MeanAnomMn = 2*pi*(MnTpres + tFlIterOrbM(1001 - i))/T_mn;
                    %this visualisation starts with points at the end of the 
                    %flight path and follows to those at the beginning
                    TrueAnomMn = TrueAnom(MeanAnomMn,e_mn,4);
                    rMn = lts_mn/(1+e_mn*cos(TrueAnomMn));
                    %translation vector:
                    r_mn(1) = rMn*cos(TrueAnomMn);
                    r_mn(2) = rMn*sin(TrueAnomMn);
                
                    if (tFlIterOrbM(1001 - i) < T_simAhead)
                        trayectory(1001 - i + iVeh0, 1) = ltsOrbM*cos(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1) + r_mn(1);
                        trayectory(1001 - i + iVeh0, 2) = ltsOrbM*sin(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1) + r_mn(2);
                        projected(1001 - i, 1) = r_mn(1);
                        projected(1001 - i, 2) = r_mn(2);
                    else
                        projected(1001 - i, 1) = ltsOrbM*cos(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1) + r_mn(1);
                        projected(1001 - i, 2) = ltsOrbM*sin(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1) + r_mn(2);
                        pValid = 1;
                    end
                end
                iVeh0 = 1000 - iMaxOrbM + iVeh0;
            end
            %-------------------------------------------------------------
        end
        
        
        %------------Visualisation in moon's coordinates------------------
        if (~retroOrbM)
            for i=1:1:iMaxOrbM
                fiVehMn = 2*pi*(i-1)/1000 + thetaVehMn; %starting point
                if (tFlIterOrbM(i) < T_simAhead)
                    trayectoryVehMn(i, 1) = ltsOrbM*cos(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1);
                    trayectoryVehMn(i, 2) = ltsOrbM*sin(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1);
                else
                    projectedVehMn(i, 1) = ltsOrbM*cos(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1);
                    projectedVehMn(i, 2) = ltsOrbM*sin(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1);
                    pValidMn = 1;
                end
            end
        else
            for i=iMaxOrbM:1:1000 %correct for retrograde trajectories (those followed clockwise)
                fiVehMn = 2*pi*i/1000 + thetaVehMn; %starting point
                if (tFlIterOrbM(1001 - i) < T_simAhead)
                    trayectoryVehMn(1001 - i, 1) = ltsOrbM*cos(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1);
                    trayectoryVehMn(1001 - i, 2) = ltsOrbM*sin(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1);
                else
                    projectedVehMn(1001 - i, 1) = ltsOrbM*cos(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1);
                    projectedVehMn(1001 - i, 2) = ltsOrbM*sin(fiVehMn)/(eOrbM*cos(fiVehMn - periapsisArgOrbM) + 1);
                    pValidMn = 1;
                end
            end
        end
        %-----------------------------------------------------------------
                        
        %------------------------------------------------------------------
        %------------------------Flight Data-------------------------------
        rVehMnOut = varargoutOrbM(1);
        thetaVehMnOut = varargoutOrbM(2);
        vVehMnOut = varargoutOrbM(3);
        psiVehMnOut = varargoutOrbM(4);
        r_periapsisOrbM = varargoutOrbM(7);
        v_periapsisOrbM = varargoutOrbM(8);
        
        %------------------------------------------------------------------
                
        if (outStatus == 'O')
            
            r_apoapsisOrbM = varargoutOrbM(9);
            v_apoapsisOrbM = varargoutOrbM(10);
                        
            Tflight = T_sim_max;
            T_sim = T_sim_max;
            MnTpres = T_sim + MnTinit;
                                
            history =  strcat(history,zone,'=>Orbiting Moon;');
            Orb = 1;
            
        else
            
            Tflight = varargoutOrbM(5);
            if (~pValid)
                T_sim = T_sim + Tflight;
            else
                T_sim = T_sim_max;
            end
            MnTpres = T_sim + MnTinit;
            
            if (outStatus == 'Z') | (outStatus == 'Zh')
                
                % crashing with the Lunar surface
                history =  strcat(history,zone,'=>',outStatus,':Collision with the surface;');
                zone = 'Z';
                
            else %'B' or 'Bh'
                
                % leaving C zone
                history =  strcat(history,zone,'=>',outStatus,';');
                outStatus = 'B'; %changing from 'Bh' to allow for direct assignement zone = outStatus
                zone = '/';
            end
        end
        %----------------Back to Earth coordinate system-------------------
        [rOut, thetaOut, vOut, psiOut] = AngCoordConv(MnTpres, 'ME', rVehMnOut, thetaVehMnOut, vVehMnOut, psiVehMnOut);
        %------------------------------------------------------------------
        if (outStatus ~= 'B')
            break;
        end
    end

end

%------------------Standarization of variables---------------------
if (~Orb)
    while (abs(thetaOut) - pi > pi),
        if (thetaOut > 2*pi) thetaOut = thetaOut - 2*pi; end
        if (thetaOut < 0) thetaOut = thetaOut + 2*pi; end
    end
        
    while (abs(psiOut) > (pi/2)),
        if (psiOut > pi) 
            psiOut = psiOut - 2*pi; 
        elseif (psiOut > pi/2)
            psiOut = pi - psiOut;
            vOut = -vOut;
        end
        if (psiOut < -pi) 
            psiOut = psiOut + 2*pi; 
        elseif (psiOut < -pi/2)
            psiOut = -(psiOut + pi);
            vOut = -vOut;
        end
    end
end
%------------------------------------------------------------------

%---------------------Lunar orbit visualisation--------------------
for i = 1:1:1000
    tMn = i/1000*(MnTpres-MnTinit);
    MeanAnom = 2*pi*(MnTinit+tMn)/T_mn;
    %true anomaly as calculated by Kepler:
    fi = TrueAnom(MeanAnom,e_mn,4);
    rMn_polar = lts_mn/(1+e_mn*cos(fi));
    rMnVis(i,1) = rMn_polar*cos(fi);
    rMnVis(i,2) = rMn_polar*sin(fi);
end
%------------------------------------------------------------------

for i=0:1:200 %silhouettes of earth's and moon's surfaces
    E(i+1,1)= R_E*cos(2*pi*i/200);
    E(i+1,2)= R_E*sin(2*pi*i/200);
    M(i+1,1)= R_mn*cos(2*pi*i/200);
    M(i+1,2)= R_mn*sin(2*pi*i/200);    
end

%----------------------------PLOTTING------------------------------

figure(1);
title('Trajectory in Terrestrial coordinates');
if (pValid)
    if (Bvisited+Cvisited)
        plot(real(trayectory(:,1)), real(trayectory(:,2)), 'b', projected(:,1), projected(:,2), 'b:',...
            E(:,1), E(:,2), 'g', rMnVis(:,1), rMnVis(:,2), 'y');
        legend('Covered path(Object)','Path to go','Earth','Covered path(Moon)',2);
    else
        plot(real(trayectory(:,1)), real(trayectory(:,2)), 'b', projected(:,1), projected(:,2), 'b:',...
            E(:,1), E(:,2), 'g');
        legend('Covered path(Object)','Path to go','Earth',2);
        if (Orb)
            if (k==1)
                title('Trajectory in Terrestrial coordinates (after: 1 orbit)');
            elseif (k>1);
                outTitle = strcat('Trajectory in Terrestrial coordinates (after:', num2str(k),' orbits)');
                title(outTitle);
            end
        end
    end
else
    if (Bvisited+Cvisited)
        plot(real(trayectory(:,1)), real(trayectory(:,2)), 'b', E(:,1), E(:,2), 'g', rMnVis(:,1), rMnVis(:,2), 'y');
        legend('Covered path(Object)','Earth','Covered path(Moon)',2);
        else
        plot(real(trayectory(:,1)), real(trayectory(:,2)), 'b', E(:,1), E(:,2), 'g');
        legend('Covered path(Object)','Earth',2);
    end
        
end
axis equal;
xlabel('metres');
ylabel('metres');

if (Cvisited)
    figure(2);
    if (pValidMn)
        plot(real(trayectoryVehMn(:,1)), real(trayectoryVehMn(:,2)), 'b', projectedVehMn(:,1), projectedVehMn(:,2), 'b:',...
            M(:,1), M(:,2), 'g');
        legend('Covered path(Object)','Path to go','Moon',2);
    else
        plot(real(trayectoryVehMn(:,1)), real(trayectoryVehMn(:,2)), 'b', M(:,1), M(:,2), 'g');
        legend('Covered path(Object)','Moon',2);
    end
    axis equal;
    xlabel('metres');
    ylabel('metres');
    if (Orb)
        if (k==1)
            title('Trajectory in Lunar coordinates (after: 1 orbit)');
        elseif (k>1);
            outTitle = strcat('Trajectory in Lunar coordinates (after:', num2str(k),' orbits)');
            title(outTitle);
        end
    else
        title('Trajectory in Lunar coordinates');
    end
end

%------------------------Clearing workspace-----------------------------
%Erasing internal data only: 

clear Bvisited Cvisited E G M M_M MeanAnom Orb fi fiVehMn i j pValid;
clear pValidMn rMn_polar tMn u outStatus outStatusOrbE outStatusOrbM;
clear iVeh0 iMaxOrbE iMaxOrbM T_simAhead zone tFlIterOrbM  tFlIterOrbE;
clear rMn MeanAnomMn0 EccAnomMn0 Tflight;

%-----------------------------------------------------------------------


