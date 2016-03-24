function [ r_XYZ, v_XYZ ] = OrbitalElements2rvGeo( h, muo, e, i, Omega, w, theta )
% this function is used to get state vectors r and v in geocentric equatorial frame of reference
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
%% INPUTS:
% h              : specific angular momentum vector in km^2/s
% muo :  Gravitational Parameter
% e              : eccentricity
% i               : inclination angle in degree
% Omega    : right ascension of the ascending node in degree
% w              : argument of perigee in degree
% theta        : true anomaly in degree
%% OUTPUTS:
% r    : geocentric state position in km
% v    : geocentric state velocity in km/s
% ---------------------------------------------------------------------------------------------------------------------------------------------------------
% r @ perifocal coordinates
r_xyz_bar=h^2/muo/(1+e*cosd(theta))*[cosd(theta);sind(theta);0];
% v @ perifocal coordinates
v_xyz_bar=muo/h*[-sind(theta);e+cosd(theta);0];
% transformation matrix from perifocal to geocentric equatorial coordinates
QxX=[-sind(Omega)*cosd(i)*sind(w)+cosd(Omega)*cosd(w),-sind(Omega)*cosd(i)*cosd(w)-cosd(Omega)*sind(w),sind(Omega)*sind(i);...
           cosd(Omega)*cosd(i)*sind(w)+sind(Omega)*cosd(w),cosd(Omega)*cosd(i)*cosd(w)-sind(Omega)*sind(w),-cosd(Omega)*sind(i);...
           sind(i)*sind(w),sind(i)*cosd(w),cosd(i)];
% geocentric  r
r_XYZ=QxX*r_xyz_bar;
% geocentric  v
v_XYZ=QxX*v_xyz_bar;
end