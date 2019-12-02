function albedo = albedo_function(theta, custom_parameters)
% ALBEDO_FUNCTION calculates the albedo at solar incidence angle theta.
%
% albedo = ALBEDO_FUNCTION(theta) calculates the albedo for array of theta
% values, theta in degrees.
%
% albedo = ALBEDO_FUNCTION(theta, custom_parameters) calculates the albedo
% for array of theta values with custom albedo parameters (e.g.
% custom_parameters = struct('a', 0.06, 'b', 0.25)

if nargin == 1
    custom_parameters = struct;
end

par = define_parameters(custom_parameters);
a = par.a;
b = par.b;
A0 = par.A0;
a_OVER_pi_OVER_4_POWER_3 = a*(1/(pi/4))^3;
b_OVER_pi_OVER_2_POWER_8 = b*(1/(pi/2))^8;

theta = deg2rad(theta);

albedo = A0 + a_OVER_pi_OVER_4_POWER_3*theta.*theta.*theta + b_OVER_pi_OVER_2_POWER_8*theta.*theta.*theta.*theta.*theta.*theta.*theta.*theta;
end