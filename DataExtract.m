clear all; close all; clc;

L   = 1;                  % Length of beam in meters
Ply = .71e-3/3;           % Ply thickness
N_P = [2;3;4;5];          % Number of Plies
F_W = [6;12;18;24]*10^-3; % Flange Width
N_B = [1;2;3;4];          % Number of batteries
B_W = 8e-3;               % Battery width
H   = 203.19997e-3;       % Height
t   = N_P*Ply;            % Laminate thickness
d   = (N_B*B_W + t)/2;    % Mid-Plane height
E_f = 65000e6;            % Modulus of fibers
E_c = 125e6;              % Modulus of core
E_C = 37;                 % Watt-Hours (Energy Capacity of a Cell) 

%% Importing Data
data = readtable('data.xlsx');

for i = 1:4
    for j = 1:4
        for k = 1:4
           
            %Get the position within the vector of recorded values
            index = 16*(i-1)+4*(j-1)+k;

            %Get the force of the roller in the center of the beam
            fRoller = table2array(data(index,3));

            %Get the vertical displacement in the center of the beam
            disp = table2array(data(index,2))*10^-3;

            %Get the mass of each model and then the energy density
            mass = table2array(data(index,4));
            eDensity(i,j,k) = 3*i*E_C/mass;
            
            d = (N_B(i)*B_W + t(j))/2; % Mid-Plane height
            Area = 4*F_W(k)*N_P(j)*Ply + N_B(i)*B_W*H + 2*N_P(j)*Ply*H; % Cross Sectional Area
            D = E_f*H*t(j)^3/6 + E_f*H*t(j)*d^2/2 + E_c*H*N_B(i)*B_W;

            %Calcualte the shear stiffness of the beam
            G(i,j,k) = ((disp - fRoller.*L.^3./(48.*D)).^-1)*fRoller.*L./(4.*Area);

        end
    end
end

%% Optimize Pack Level Energy Density

[maxDensity, optLoc] = optimize(-eDensity, G);

function [maxDensity, minPos] = optimize(battery, stiffness)

buffer    = ones(6,6,6)*1e9; 
debug     = ones(6,6,6);

%Apply penalty method to the battery data
for i = 1:4
    for j = 1:4
        for k = 1:4

            buffer(i+1, j+1, k+1) = battery(i,j,k);
            debug(i+1, j+1, k+1) = battery(i,j,k);

            if stiffness(i,j,k) < 2e8

                buffer(i+1,j+1,k+1) = battery(i,j,k) + 1e17/(stiffness(i,j,k));
                debug(i+1,j+1,k+1) = battery(i,j,k) + 1e11/(stiffness(i,j,k));

            end

       
        end
    end
end

batteryNew = buffer; 

%function [optVal] = hookJeeves(battery)

x0 = [5 5 5];
xCurrent = x0;

xPlus  = [1 0 0; 
          0 1 0; 
          0 0 1]; 

xMinus = -xPlus; 

for i = 1:128

    fCurrent = batteryNew(xCurrent(1), xCurrent(2), xCurrent(3)); 

    for j = 1:3

        %Build the pattern of search
        xPosP(j,:) = xCurrent + xPlus(j,:);
        fPlus(j) = batteryNew(xPosP(j,1), xPosP(j,2), xPosP(j,3));

        xPosM(j,:) = xCurrent + xMinus(j,:);
        fMinus(j) = batteryNew(xPosM(j,1), xPosM(j,2), xPosM(j,3));
           
    end

    xNew = vertcat(xCurrent, xPosP, xPosM);
    fNew = [fCurrent fPlus fMinus];

    [minVal, minPos] = min(fNew);

    xCurrent = xNew(minPos,:);


end


maxDensity = -minVal;
minPos = xCurrent;


end 
