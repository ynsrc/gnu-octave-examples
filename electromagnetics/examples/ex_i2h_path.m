% clear screen and variables
clear, clc

% import YNSRC EFT Library
em_ynsrc;

I = 10;
P = [0 0 0; 8 0 0; 8 4 0; 0 4 0; 0 0 0];

C = [2 2 0];
em_calc_h(I, P, C);

C = [4 2 0];
em_calc_h(I, P, C);

C = [4 8 0];
em_calc_h(I, P, C);

C = [0 0 2];
em_calc_h(I, P, C);

