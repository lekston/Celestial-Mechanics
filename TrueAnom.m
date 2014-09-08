function Theta = TrueAnom(MeanAnom, eps, I)

%First below algorithm calculates numerically the eccentric anomaly using
%Newton-Raphson's Method
%from Kepler's procedure: MeanAnom(t) = EccAnom - eps*sin(EccAnom)
%input variables: 
%MeanAnom - mean anomaly, 0 <= MeanAnom <= 2*pi 
%eps - orbit eccentricity
%I - number of iterations (usually 3 - 4 will suffice)


EccAnom(1) = MeanAnom;

for i=2:1:(I+1)
    EccAnom(i) = EccAnom(i-1) - (MeanAnom - EccAnom(i-1) + eps*sin(EccAnom(i-1)))/(eps*cos(EccAnom(i-1))-1);
end

EccAnomOut = EccAnom(I+1);

%and now for the True anomaly:
Theta = 2*atan2((sqrt(1+eps)*sin(EccAnomOut/2)), (sqrt(1-eps)*cos(EccAnomOut/2))); % Theta = 2*arg(Y,X);
%atan2(Y,X) finds the angular coordinate of the point in a range from -pi to pi.