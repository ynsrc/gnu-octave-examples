% clear screen and variables
clear, clc

% import YNSRC EFT Library
em_ynsrc;

%%%%%%%%%%%% EXAMPLE 1 %%%%%%%%%%%%%

I = 6;
A = [1 0 0];
B = [0 1 0];
C = [0 0 0];

csys = 'Cartesian';

H = em_i2h(I, A, B, C);

fprintf("%g A current flows in path %s->%s\nH at %s point = %s\n\n", I,...
  vec2strd(A,csys), vec2strd(B,csys), vec2strd(C,csys),...
  num2engvec(H, csys, 'A/m'));

%%%%%%%%%%%%% EXAMPLE 2 %%%%%%%%%%%%

I = 10;
P = [0 0 0; 2 0 0; 1 1 0; 0 0 0];
C = [0 0 5];

em_calc_h(I,P,C);

