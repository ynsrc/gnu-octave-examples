% clear screen and variables
clear, clc

% import YNSRC EFT Library
em_ynsrc;

% half infinity

I = 3;
A = [ 0 0 Inf ];
B = [ 0 0 0 ];
C = [-3 4 0 ];

H1 = em_i2h(3, A, B, C);

fprintf("\n%g A current flows in path %s->%s causes at %s\n\tH1 = %s\n", I,...
  vec2strd(A,'Cartesian'), vec2strd(B,'Cartesian'), vec2strd(C,'Cartesian'),...
  num2engvec(H1, 'Cartesian', 'A/m')
);

I = 3;
A = [ 0 0 0 ];
B = [ Inf 0 0 ];
C = [-3 4 0 ];

H2 = em_i2h(3, A, B, C);

fprintf("\n%g A current flows in path %s->%s causes at %s\n\tH2 = %s\n", I,...
  vec2strd(A,'Cartesian'), vec2strd(B,'Cartesian'), vec2strd(C,'Cartesian'),...
  num2engvec(H2, 'Cartesian', 'A/m')
);

HT = H1 + H2;

fprintf("\nTotal (HT) = H1+H2 = %s\n", num2engvec(H1+H2, 'Cartesian', 'A/m'));

