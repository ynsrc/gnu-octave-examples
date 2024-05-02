% clear screen and variables
clear, clc

% import YNSRC EFT Library
em_ynsrc;

coordinate_system = 'Cylindrical';

I = 10;  % current
R = 3;   % radius
h = 4;   % find H at this point

L = [ 0 sym('2*pi') ]; % full circle

k = I/(4*sym(pi)*(R^2+h^2)^(3/2));

H = double([ 0 0 int(k*R^2, phi, L(1), L(2)) ]);

fprintf("In r=%g radius circle at (z=0) flowing %g A current ", R, I);
fprintf("at h=%g point\ncauses H = %s A/m\n", h,...
  vec2str(H, coordinate_system));

