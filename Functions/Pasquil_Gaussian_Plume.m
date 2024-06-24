function [ C ] = Pasquil_Gaussian_Plume( s, p )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

s.D=s.ci;
s.t = s.cii;
lamda = sqrt((s.D.*s.t)./(1+ (s.u.^2.*s.t)./(4*s.D)));

module_dist = sqrt(((s.x-p.x_matrix)).^2 + ((s.y-p.y_matrix)).^2 + ((s.z-p.z_matrix)).^2);
angazimuth= atan2((p.x_matrix-s.x),(p.y_matrix-s.y));
angazimuth(angazimuth<0)=angazimuth(angazimuth<0)+2*pi;
ang_diff_azimuth = angazimuth+s.phi-pi/2;
ang_diff_azimuth(ang_diff_azimuth<0)=ang_diff_azimuth(ang_diff_azimuth<0)+2*pi;
xr = module_dist .* cos(ang_diff_azimuth);

% C = s.Q./(4*pi.*s.D.*module_dist).*exp(1*(xr).*s.u./(2.*s.D)).*exp(-1*module_dist./lamda);
C = s.Q./(4*pi.*s.D.*module_dist).*exp((1*(xr).*s.u./(2.*s.D))+(-1*module_dist./lamda));


end

