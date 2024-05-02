# YNSRC - Electromagnetic Field Theory (with GNU Octave)

This project contains electormagnetic field theory calculations.

# How To Run
You can use [Octave Online](https://octave-online.net) to run this examples.
Or you can download and install GNU Octave application into your computer.

# Usage
Type as command or place after typical "clear, clc;" line at the beginning of script file.
```
em_ynsrc;
```

This will load all symbolic variables and functions into workspace.

# Examples

You can look at example usages in folder located `examples` next to em_ynsrc.m file.

Or type a function name without parameters, that will cause to show you function's
help text from `em_ynsrc_help.txt` if located next to `em_ynsrc.m` file.

# Using Help
If you copy `em_ynsrc_help.txt` file next to `em_ynsrc.m` file, when you type a
function name without parameters, that causes to show you function help text
from help file.

If you type `int2` as command help will show as;
```
R = int2(v)

Computes 2-fold (double) integral.

E.g.: R = int2([ x*y, x 1 3, y 2 5 ])

The above code integrates the function x*y with respect to x and y.
The x boundaries are [1 3] and the y boundaries are [2 5].
```

or type `em_calc_gauss` as command and help will show

```
[fS,psi] = em_calc_gauss(D, [N, L, L2, coord_sys])

1. Finds the volume charge density (pV) from the electric flux density (D).
2. Calculates the pV(N) density value at the specified point N.
3. Calculates the number of surfaces of the object whose boundaries are given by L.
4. Finds the fluxes (fS) coming from these surfaces and their sum.
5. Verifies the divergence theorem.
6. Calculates the flux coming out of the surface limited to L2.

Example:
  D = [ 3*r 2*cos(phi) r*sin(theta) ];
  L = [ 0 3 0 sym(pi) 0 2*sym(pi) ];
  L2 = [ 3 3 0 sym(pi)/2 0 sym(pi)/2 ];
  N = [ 1 sym(pi)/6 sym(pi)/3 ];
  [fS,psi] = em_calc_gauss(D,L,N,L2);

Result:
Electric flux density (D):
  [3*r  2*cos(phi)  r*sin(theta)]

a) Volume charge density pV (∇.D):
      2*cos(phi)*cos(theta)
  9 + ---------------------
           r*sin(theta)

b) Volume charge density at point (1,30,60) (pV) = 10.7321 C/m³

c) 0<=r<=3, 0<=θ<=180, 0<=φ<=360 bounded sphere
   There are 1 surface(s) in total.

 âr directed electric flux ∫(D)dS,dS=r²sin(θ)dθdφ(âr) = 1.01788 kC
   Ψ1 = 324*pi [C]


Left side of the equation ∮(D)dS = Ψ = Q = 1.01788 kC
  324*pi

Right side of the equation ∫(pV)dV = Ψ = Q = 1.01788 kC
  324*pi

Since left is equal to right, Divergence theorem is satisfied.

e) Flux leaving the bounded surface 3<=r<=3, 0<=θ<=90, 0<=φ<=90
   Ψ = Q = 127.235 C
```

# License
CC BY 4.0 LEGAL CODE
Attribution 4.0 International

Feel free to use this code in your personal, open-source or even commercial projects. Only attribution needed.
