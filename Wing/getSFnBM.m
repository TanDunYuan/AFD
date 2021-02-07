function [SF,BM]=getSFnBM(stationForceDistribution,stationSpanmesh,FAAngle)
% obtain the Shear Force and BM of each station 
%outpu 

Nstations=length(stationSpanmesh);
SF=zeros(1,Nstations);
BM=zeros(1,Nstations);
dM_inertial=zeros(1,Nstations);

stationsFA= stationSpanmesh ./ cosd(FAAngle); %get FA length in m at each station (this will be the moment arm for each station)

% get spanwise shear force distribution
for i = 1:Nstations-1
    SF(Nstations-i) = sum(stationForceDistribution(end-i:end));
end

% get bending moment contribution for each station
for i = 1:Nstations-1
    dM_inertial(i) = (SF(i)+SF(i+1)) * (stationsFA(i+1)-stationsFA(i)) / 2;
end

% get spanwise bending moment distribution
for i = 1:Nstations-1
    BM(length(dM_inertial)-i) = sum(dM_inertial(end-i:end));
end

end