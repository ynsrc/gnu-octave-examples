% clear screen and variables
clear, clc

% import YNSRC EFT Library
em_ynsrc;

I = 2;        % current (A)
A = [0 0 0];  % start point
B = [0 0 10]; % final point

C = [5 0 0];  % target

fprintf("\n%g A current flows in path %s->%s H at %s\n\tH = %s\n", I,...
  vec2strd(A,'Cartesian'), vec2strd(B,'Cartesian'), vec2strd(C,'Cartesian'),...
  num2engvec(em_i2h(I, A, B, C),'Cartesian', 'A/m')
);

C = [5 5 0];  % target

fprintf("\n%g A current flows in path %s->%s H at %s\n\tH = %s\n", I,...
  vec2strd(A,'Cartesian'), vec2strd(B,'Cartesian'), vec2strd(C,'Cartesian'),...
  num2engvec(em_i2h(I, A, B, C),'Cartesian', 'A/m')
);

C = [5 15 0];  % target

fprintf("\n%g A current flows in path %s->%s H at %s\n\tH = %s\n", I,...
  vec2strd(A,'Cartesian'), vec2strd(B,'Cartesian'), vec2strd(C,'Cartesian'),...
  num2engvec(em_i2h(I, A, B, C),'Cartesian', 'A/m')
);

C = [5 -15 0];  % target

fprintf("\n%g A current flows in path %s->%s H at %s\n\tH = %s\n", I,...
  vec2strd(A,'Cartesian'), vec2strd(B,'Cartesian'), vec2strd(C,'Cartesian'),...
  num2engvec(em_i2h(I, A, B, C),'Cartesian', 'A/m')
);

