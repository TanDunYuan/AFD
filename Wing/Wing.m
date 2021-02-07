clear 
clc

%% Default Settings for  Figures
set(0,'defaulttextInterpreter','Latex','DefaultAxesFontSize', 15);

%% Load variables
var=load_mainWing();

%% Find station perimeters
% Discretise wing
station.SpanMesh = 0:var.spandx:var.SemiSpanActual;
Nstations=length(station.SpanMesh);
stations_FA = station.SpanMesh ./ cosd(var.FAAngle); %get FA length in m at each station (this will be the moment arm for each station)

% Import airfoil coordinates
fid = fopen(var.datfileLoc, 'r');
x_airfoil = fscanf(fid,'%f %f',[2 inf]);
y_airfoil = x_airfoil(2,[1:end]);
x_airfoil = x_airfoil(1,[1:end]);
fclose(fid);

% Get cross-sectional area of each station
station.Chord = var.GeoRootChord - (var.GeoRootChord-var.TipChord) * (var.FuseRad+station.SpanMesh) / var.SemiSpan; %get station chords
for i = 1:length(station.Chord) %scale airfoil coordinates for each station and put it to global airplane xyz coord
    % where y is direction of span, z is direction of gravity
    station.x_airfoil(i,:) = station.Chord(i) .* x_airfoil;
    station.z_airfoil(i,:) = station.Chord(i) .* y_airfoil;
    station.y_airfoil(i,:)=ones(1,length(station.x_airfoil(i,:))) *station.SpanMesh(i);

    %apply the function polygeom
    [geom(i,:),iner,cpmo] = polygeom(station.x_airfoil(i,:),station.z_airfoil(i,:)); %geom(:,1) is the cross-sectional area of each station in m^2
end

% Wing mass for each station based on volume ratio and Wing Mass
station.Vol = geom(:,1)' * var.spandx / cosd(var.FAAngle); %get station volume in m^3
station.totalVol = sum(station.Vol);%get total discretised wing volume
station.StructMass = station.Vol / station.totalVol * var.SemiWingMass;%station estimated mass based on Wing Mass in var

% Fuel mass for each station based on perc Span the fuel takes up
wing_FuelVol = sum(station.Vol(1:round(Nstations*var.FuelpercSpan))); %find wing volume used for fuel
station.FuelMass = zeros(1,length(station.StructMass)); %ensure same array size as station.SpanMesh
station.FuelMass(1:round(Nstations*var.FuelpercSpan) ) = station.Vol(1:round(Nstations*var.FuelpercSpan)) / wing_FuelVol * var.FuelMass; 

% find the station which holds the uc
ucmass_y_wing=var.UCyPos-var.FuseRad;
[~,station_idx_UC]=min(abs(ucmass_y_wing-station.SpanMesh)); 

% find the station which holds the engine
enginemass_y_wing=var.EngineYPos-var.FuseRad;
[~,station_idx_Engine]=min(abs(enginemass_y_wing-station.SpanMesh)); 

% total mass contain the struct, fuel and single loads such as Engine and
% UC
station.TotalMass = station.StructMass + station.FuelMass;
station.TotalMass(station_idx_Engine)=station.TotalMass(station_idx_Engine)+var.EngineMass;
station.TotalMass(station_idx_UC)=station.TotalMass(station_idx_UC)+var.UCMass;

%% Plot wing shape
figure

global_wing.y=station.y_airfoil;
global_wing.z=station.z_airfoil;
for i =1:Nstations
global_wing.x(i,:)=station.x_airfoil(i,:)+station.SpanMesh(i)*tand(var.Sweep);
end
for i=1:100:length(station.SpanMesh)
    plot3(global_wing.y(i,:),global_wing.x(i,:),station.z_airfoil(i,:),"k")
    hold on
end
plot3(global_wing.y(1:200:end,:),global_wing.x(1:200:end,:),global_wing.z(1:200:end,:),"k")

axis equal
grid on 

%% Find Aerodynamic Loads
aero_Lift = var.n * var.MZFW * 9.81 / 2; %lift required per wing in N
% get lift for each station
aero_L0 = 1000; %initial guess for L0
dL = aero_L0 * sqrt(1-(station.SpanMesh./var.SemiSpanActual).^2); 
aero_LiftError = abs(sum(dL) - aero_Lift);

% Simple gradient descent to find dL&L0
while aero_LiftError > 0.1 
    if sum(dL) > aero_Lift
        aero_L0 = aero_L0 - aero_LiftError/100/100;
    else
        aero_L0 = aero_L0 + aero_LiftError/100/100;
    end
    dL = aero_L0 * sqrt(1-(station.SpanMesh./var.SemiSpanActual).^2); 
    
    aero_LiftError = abs(sum(dL) - aero_Lift);
end

%% inertialF Loads
[SF.inertial,BM.inertial]=getSFnBM(-station.TotalMass.*9.81 ,station.SpanMesh,var.FAAngle);

[SF.aero,BM.aero]=getSFnBM(dL ,station.SpanMesh,var.FAAngle);

%% Overall loading and Torque Distribution
SF.tot = SF.aero + SF.inertial;
BM.tot = BM.aero + BM.inertial;

% get aerodynamic pitching moment for each station
T_M0 = 0.5 * var.rho * var.Vmo^2 * station.Chord.^2 * var.Cm0; %formula from notes p54

% get the moment arms for each stations
station.a = (var.FALocation - 0.25).* station.Chord; %lift moment arm 
station.b = (var.CGLocation - var.FALocation).* station.Chord; %weight moment arm (CG - flexural axis)

% get total torque
dT = station.a.*dL + station.b.*station.StructMass * 9.81 + T_M0; %N (+T_M0 because Cm0 is -ve )
T=zeros(1,Nstations);
for i = 1:Nstations-1
    T(length(dT)-i) = sum(dT(end-i:end));
end
%% Plots

% plot force distribution
figure_Force=figure; 
hold on
area1 = area(station.SpanMesh,dL);
area1.FaceColor = [0 0 1];
area1.FaceAlpha = 0.3;
area2 = area(station.SpanMesh,-station.StructMass*9.81);
area2.FaceColor = [1 0 0];
area2.FaceAlpha = 0.3;
area3= area(station.SpanMesh,-station.FuelMass*9.81);
area3.FaceColor = [ 0 1 0];
area3.FaceAlpha = 0.3;

plot([station.SpanMesh(station_idx_UC),station.SpanMesh(station_idx_UC)],[0,-var.UCMass*9.81],'k',"LineWidth",2)
plot([station.SpanMesh(station_idx_Engine),station.SpanMesh(station_idx_UC)],[0,-var.EngineMass*9.81],'r',"LineWidth",2)
ylim([-200 1000]) 
ylabel('Force Distribution')
xlabel('Span station Y (m)')
grid on
legend({'Lift','Wing Self Weight','FuelMass','MainUC Weight','Engine Weight'},'location','Best') 

annotation(figure_Force,'textarrow',[0.376071428571429 0.211785714285714],...
    [0.397619047619048 0.186190476190476],'String',num2str(round(var.UCMass*-9.81,0))+" N","FontSize",15);

annotation(figure_Force,'textarrow',[0.431785714285714 0.313928571428572],...
    [0.325238095238095 0.175714285714286],'Color',[1 0 0],'String',num2str(round(var.EngineMass*-9.81,0))+" N","FontSize",15);

%% Plot BM SF T distribution
figure
subplot(3,1,1)
hold on
grid on
ylabel('SF (N)')
plot(station.SpanMesh,SF.tot,'b')
% plot BM
subplot(3,1,2)
hold on
grid on
ylabel('BM (Nm)')
plot(station.SpanMesh,BM.tot,'b')
% plot T
subplot(3,1,3)
hold on
grid on
ylabel('Torque (Nm)')
xlabel('Span stations (m)')
plot(station.SpanMesh,T,'b')

%%  find spar coordinates in %chord 
 
% find the FS and RS locations (in terms of idx of airfoil)
[~,yRS_idx]=mink(abs(x_airfoil-var.RSLoc),2);
[~,yFS_idx]=mink(abs(x_airfoil-var.FSLoc),2);

% h and l are the wing box coordinates in terms of perc chord
h = min([abs(diff(y_airfoil(yRS_idx))) abs(diff(y_airfoil(yFS_idx))) ]); % use lower of the two as the height of rectangular wing box
l = var.RSLoc - var.FSLoc; % distance between spar as perc of chord


figure
plot(x_airfoil ,y_airfoil )
line([var.FSLoc,var.RSLoc,var.RSLoc,var.FSLoc,var.FSLoc],[ y_airfoil(yFS_idx(2))  flip(y_airfoil(yFS_idx)) y_airfoil(yFS_idx)],'LineWidth',0.5,'Color','k');hold on
line([var.FSLoc,var.RSLoc,var.RSLoc,var.FSLoc,var.FSLoc],[ y_airfoil(yRS_idx(2))  flip(y_airfoil(yRS_idx)) y_airfoil(yRS_idx)],'LineWidth',0.5,'Color','r');
axis equal
title("wing box based on FS/RS")


% wing box dimensions
wingbox.c=l*station.Chord';
wingbox.b2=h*station.Chord';

