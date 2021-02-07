function var=load_mainWing()
% define the Main Wing relevant variables here 

% Weights and or position
var.SemiWingMass = 7743.93 / 2; % single wing
var.FuelMass = 0 ; % (worst case assumes no fuel)
var.FuelpercSpan = 0.7;   % Percentage of Span which fuel is occupying

var.spandx=0.01;  %dx for station spacing

var.MZFW = 106125-37588; %

var.UCMass = 3086.29/2 ; % 
var.UCyPos = 4.24; %  m from fuse centre line

var.EngineMass = 2952.9; % 
var.EngineYPos = 7; % m from fuse centre line


% Fuselage Params
var.FuseRad = 5.47/2; %m

% Wing Params
var.datfileLoc="n63412.dat";

var.WingSpan = 41.73;
var.SemiSpan = var.WingSpan/2; 
var.TaperRatio = 0.3;
var.Sweep = 22;
var.GeoRootChord = 6.978; % geometric root chord
var.TipChord = 2.0935; 
var.FSLoc = 0.1; % in perc chord
var.RSLoc = 0.7; %in perc chord

var.SemiSpanActual = var.SemiSpan - var.FuseRad; 
var.RootChordActual = var.GeoRootChord - (var.GeoRootChord-var.TipChord) * var.FuseRad / var.SemiSpan; 

% Load Case Considered
var.n = 2.5 * 1.5; % CFR 25.333  * 1.5 for ultimate load 

% Flight Conditions
var.rho = 0.3796; % cruise altitude
var.Vmo = 237.28 ; % cruise condition velocity 
var.Cm0 = -0.08; % moment coefficient of airfoil at cruise AoA 
var.CGLocation = 0.4; % initial assumption of wing CG %chord 

% Material properties (Spar: Al-Li 8090, Nose: Al-Zn 7055)
% Spar
var.Spar.E = 77e9; 
var.Spar.rho = 2540; 
var.Spar.nu = 0.34;
var.Spar.sigmay = 370e6; %Pa


% Skin Stringer
var.SS.E = 71.7e9; %Pa
var.SS.rho = 2860; %/m^3

var.Rib.rho = 2540; %/m^3


var.FALocation = (var.FSLoc + var.RSLoc)/2; %chord
var.FAAngle = 20.5; % deg

end