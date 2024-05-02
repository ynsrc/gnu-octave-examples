# YNSRC - Electromagnetic Field Theory (with GNU Octave)

This folder contains example usage of YNSRC EFT library.

# How To Run
You can use [Octave Online](https://octave-online.net) to run this examples.
Or you can download and install GNU Octave application into your computer.

# Usage
Type as command or place after typical "clear, clc;" line at the beginning of script file.
```
em_ynsrc;
```

This will load all symbolic variables and functions into workspace.

# Outputs of Examples (from Octave Online)

## ex_i2h.m
```
6 A current flows in path (1,0,0)->(0,1,0)
H at (0,0,0) point = +954.9297(âz)  mA/m

10 A current flowing in the conductor;

 (0,0,0)->(2,0,0) path
 H(0,0,5) = -0.059109(ây)  A/m

 (2,0,0)->(1,1,0) path
 H(0,0,5) = +27.3651(âx) +27.3651(ây) +10.946(âz)  mA/m

 (1,1,0)->(0,0,0) path
 H(0,0,5) = -30.6294(âx) +30.6294(ây)  mA/m

Total H(0,0,5) = -3.2643(âx) -1.1142(ây) +10.946(âz)  mA/m
```

## ex_i2h_2.m
```
2 A current flows in path (0,0,0)->(0,0,10) H at (5,0,0)
	H = +28.4705(ây)  mA/m

2 A current flows in path (0,0,0)->(0,0,10) H at (5,5,0)
	H = -12.9949(âx) +12.9949(ây)  mA/m

2 A current flows in path (0,0,0)->(0,0,10) H at (5,15,0)
	H = -5.1043(âx) +1.7014(ây)  mA/m

2 A current flows in path (0,0,0)->(0,0,10) H at (5,-15,0)
	H = +5.1043(âx) +1.7014(ây)  mA/m
```

## ex_i2h_circle.m
```
In r=3 radius circle at (z=0) flowing 10 A current at h=4 point
causes H = +0.36(âz)  A/m
```

## ex_i2h_halfinf.m
```
3 A current flows in path (0,0,Inf)->(0,0,0) causes at (-3,4,0)
	H1 = -0.038197(âx) -0.028648(ây)  A/m

3 A current flows in path (0,0,0)->(Inf,0,0) causes at (-3,4,0)
	H2 = +23.8732(âz)  mA/m

Total (HT) = H1+H2 = -38.1972(âx) -28.6479(ây) +23.8732(âz)  mA/m
```

## ex_i2h_path.m
```
10 A current flowing in the conductor;

 (0,0,0)->(8,0,0) path
 H(2,2,0) = +658.8179(âz)  mA/m

 (8,0,0)->(8,4,0) path
 H(2,2,0) = +83.882(âz)  mA/m

 (8,4,0)->(0,4,0) path
 H(2,2,0) = +658.8179(âz)  mA/m

 (0,4,0)->(0,0,0) path
 H(2,2,0) = +562.6977(âz)  mA/m

Total H(2,2,0) = +1.9642(âz)  A/m

10 A current flowing in the conductor;

 (0,0,0)->(8,0,0) path
 H(4,2,0) = +711.7625(âz)  mA/m

 (8,0,0)->(8,4,0) path
 H(4,2,0) = +177.9406(âz)  mA/m

 (8,4,0)->(0,4,0) path
 H(4,2,0) = +711.7625(âz)  mA/m

 (0,4,0)->(0,0,0) path
 H(4,2,0) = +177.9406(âz)  mA/m

Total H(4,2,0) = +1.7794(âz)  A/m

10 A current flowing in the conductor;

 (0,0,0)->(8,0,0) path
 H(4,8,0) = +88.9703(âz)  mA/m

 (8,0,0)->(8,4,0) path
 H(4,8,0) = +37.2662(âz)  mA/m

 (8,4,0)->(0,4,0) path
 H(4,8,0) = -0.28135(âz)  A/m

 (0,4,0)->(0,0,0) path
 H(4,8,0) = +37.2662(âz)  mA/m

Total H(4,8,0) = -0.11785(âz)  A/m

10 A current flowing in the conductor;

 (0,0,0)->(8,0,0) path
 H(0,0,2) = -0.38601(ây)  A/m

 (8,0,0)->(8,4,0) path
 H(0,0,2) = +10.2148(âx) +40.8594(âz)  mA/m

 (8,4,0)->(0,4,0) path
 H(0,0,2) = +69.4609(ây) +138.9218(âz)  mA/m

 (0,4,0)->(0,0,0) path
 H(0,0,2) = -0.35588(âx)  A/m

Total H(0,0,2) = -345.6664(âx) -316.5465(ây) +179.7812(âz)  mA/m
```