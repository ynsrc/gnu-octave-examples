% YNSRC - ELECTROMAGNETIC FIELD THEORY LIBRARY

% import required symbols
pkg load symbolic;
syms x y z rho phi r theta epsilon_0;

'%%%%%%%%%%%%%%%% REQUIRED HELPER FUNCTIONS %%%%%%%%%%%%%%';

% show help for caller function
function show_ynsrc_help()
  if ~exist("em_ynsrc_help.txt")
    printf("Help file em_ynsrc_help.txt not found!");
    return;
  end

  fnc = dbstack(1).name;
  file = fopen("em_ynsrc_help.txt", "r");

  line = fgetl(file);
  found = strcmp(line, ["#begin " fnc]);

  while ~found && ~feof(file)
    line = fgetl(file);
    found = strcmp(line, ["#begin " fnc]);
  end

  if found
    line = fgetl(file);
    while ~startsWith(line, "#end") && ~feof(file)
      disp(line);
      line = fgetl(file);
    end
  end

  fclose(file);
end

% 2-fold (double) integral
function R = int2(v)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  R = int(int(v(1),v(2),v(3),v(4)),v(5),v(6),v(7));
end

% Triple  integral
function R = int3(v)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  R = int(int2(v(1:7)),v(8),v(9),v(10));
end

% Integral and substitute
function R = int1subs(v)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  R = em_subs(int(v(1),v(2),v(3),v(4)), [v(5) v(7)], [v(6) v(8)]);
end

% 2-fold integral and substitute
function R = int2subs(v)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  R = em_subs(int2(v(1:7)), v(8), v(9));
end

% Convert number to engineering notation string (e.g. 3 uC, 5 mF)
function [engstr] = num2eng(value, unit="")
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  value = double(value);
  if value == 0
      engstr = sprintf("%d %s", value, unit);
      return
  end

  ranges = "fpnum kMGTPE";
  range = 7;

  k1 = floor(log10(value));
  k2 = floor(k1 / 3);

  dispval = value / 10^(k2*3);

  if k2 >= -5 && k2 <= 6
      range_t = strtrim(ranges(k2+6));
      engstr = sprintf("%d %s%s", dispval, range_t, unit);
  else
      engstr = sprintf("%dx10^%d %s", dispval, k1, unit);
  end
end

% If V is a vector, it gives the coordinate system,
% and if it is text, it gives its symbols.
function cs = vdetcoordsys(V)
if nargin == 0; show_ynsrc_help(); return; end; % show help
  if ischar(V)
    switch V
      case { 'cartesian', 'Cartesian' }
        syms x y z;
        cs = [ x y z ];
      case { 'cylindrical', 'Cylindrical' }
        syms rho phi z;
        cs = [ rho phi z ];
      case { 'spherical', 'Spherical' }
        syms r theta phi;
        cs = [ r theta phi ];
    end
  else
    syms x y z rho phi r theta;
    vars = symvar(V);
    cs = 'Cartesian';
    if any([ has(vars, x) has(vars, y) ])
      cs = 'Cartesian';
      return;
    elseif any([ has(vars, rho) ])
      cs = 'Cylindrical';
      return;
    elseif any([ has(vars, r) has(vars, theta) ])
      cs = 'Spherical';
      return;
    elseif any([ has(vars, phi) ])
      if any(has(vars, z))
        cs = 'Cylindrical';
      elseif any(not([ has(vars, r) has(vars, theta) ]))
          cs = 'Cylindrical';
          warning("Assuming coordinate system is Cylindrical!");
      else
        warning("Assuming coordinate system is Cylindrical!");
        cs = 'Spherical';
      end
    else
      warning("Assuming coordinate system is Cartesian!");
      cs = 'Cartesian';
    end
  endif
end

% Divergence
function Div = divergence_sym(V,X,coordinate_system)
if nargin == 0; show_ynsrc_help(); return; end; % show help
if nargin < 3; coordinate_system=vdetcoordsys(V); end;
switch coordinate_system
    case {'cartesian','Cartesian'}
        Div=divergence(V,X);
    case {'cylindrical','Cylindrical'}
        Div=1/X(1)*diff(X(1)*V(1),X(1))+1/X(1)*diff(V(2),X(2))+diff(V(3),X(3));
    case {'spherical','Spherical'}
        Div=1/X(1)^2*diff(X(1)^2*V(1),X(1))+1/(X(1)*sin(X(2)))*diff(sin(X(2))*V(2),X(2))+1/(X(1)*sin(X(2)))*diff(V(3),X(3));
end
end

% Gradient
function gradientSym = gradient_sym(V,X,coordinate_system)
if nargin == 0; show_ynsrc_help(); return; end; % show help
if nargin < 3; coordinate_system=vdetcoordsys(V); end;
switch coordinate_system
    case {'cartesian','Cartesian'}
        gradientSym=gradient(V,X);
    case {'cylindrical','Cylindrical'}
        gradientSym=[diff(V,X(1)),1/X(1)*diff(V,X(2)),diff(V,X(3))];
    case {'spherical','Spherical'}
        gradientSym=[diff(V,X(1)),1/X(1)*diff(V,X(2)),1/(X(1)*sin(X(2)))*diff(V,X(3))];
end
end

% Rotational (Curl)
function CurlSym = curl_sym(V,X,coordinate_system)
if nargin == 0; show_ynsrc_help(); return; end; % show help
if nargin < 3; coordinate_system=vdetcoordsys(V); end;
switch coordinate_system
    case {'cartesian','Cartesian'}
        CurlSym=curl(V,X);
    case {'cylindrical','Cylindrical'}
        CurlSym=[1/X(1)*diff(V(3),X(2))-diff(V(2),X(3)),...
            diff(V(1),X(3))-diff(V(3),X(1)),...
            1/X(1)*diff(X(1)*V(2),X(1))-1/X(1)*diff(V(1),X(2))];
    case {'spherical','Spherical'}
        CurlSym=[1/(X(1)*sin(X(2)))*(diff(V(3)*sin(X(2)),X(2))-diff(V(2),X(3))),...
            1/X(1)*(1/sin(X(2))*diff(V(1),X(3))-diff(X(1)*V(3),X(1))),...
            1/X(1)*(diff(X(1)*V(2),X(1))-diff(V(1),X(2)))];
end
end

% rot_sym is just a reference for curl_sym
rot_sym = @curl_sym;

'%%%%%%%%%%%%%%%%%%%%%%%%%%%% VECTORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%';

% Calculates by substituting variables such as epsilon_0 into the equation
function R = em_subs(F, X=[], Y=[])
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  syms epsilon_0;
  R = subs(F, [X epsilon_0], [Y sym('10^-9/(36*pi)')]);
end

% Places the values of symbols such as [rho phi z] in the expression F
function R = em_subs_vec(F, Y=[], coordinate_system)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  if nargin < 3; coordinate_system=vdetcoordsys(F); end;
  R = em_subs(F, vdetcoordsys(coordinate_system), Y);
end

% Like em_subs_vec but angles must be entered in degrees
function R = em_subs_vecd(F, Y, coordinate_system)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  if nargin < 3; coordinate_system=vdetcoordsys(F); end;
  X = vdetcoordsys(coordinate_system);
  switch coordinate_system
    case { 'cylindrical', 'Cylindrical' }
      Y = vcyl(Y);
    case { 'spherical', 'Spherical' }
      Y = vsph(Y);
  end
  R = em_subs(F, X, Y);
end

% vector in cartesian coordinate system
function R = vcart(varargin)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  if nargin == 1
    V = varargin{1};
    R = [ sym(num2str(V(1))) sym(num2str(V(2))) sym(num2str(V(3))) ];
  elseif nargin == 3
    R = [ sym(num2str(varargin{1})) sym(num2str(varargin{2})) ...
          sym(num2str(varargin{3})) ];
  else
    error("vcart was used incorrectly!");
  end
end

% vector in cylindrical coordinate system (degree angle)
function R = vcyl(varargin)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  if nargin == 1
    V = varargin{1};
    for i = 1:rows(V)
      R(i,1:3) = [ sym(num2str(V(i,1))), sym(num2str(V(i,2)))*sym(pi)/180, ...
                   sym(num2str(V(i,3))) ];
    end
  elseif nargin == 3
    k = sym('pi/180');
    R = [ sym(num2str(varargin{1})) k*sym(num2str(varargin{2})) ...
          sym(num2str(varargin{3})) ];
  else
    error("vcyl was used incorrectly!");
  end
end

% vector in global coordinate system (degree angle)
function R = vsph(varargin)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  if nargin == 1
    V = varargin{1};
    for i = 1:rows(V)
      R(i,1:3) = [ sym(num2str(V(i,1))), sym(num2str(V(i,2)))*sym(pi)/180, ...
                   sym(num2str(V(i,3)))*sym(pi)/180 ];
    end
  elseif nargin == 3
    k = sym('pi/180');
    R = [ sym(num2str(varargin{1})) k*sym(num2str(varargin{2})) ...
          k*sym(num2str(varargin{3})) ];
  else
    error("vsph was used incorrectly!");
  end
end

% Converts vector to sym expression
function R = vec2sym(V, coordinate_system)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  switch coordinate_system
    case {'cartesian', 'Cartesian'}
      R = vcart(V);
    case {'cylindrical', 'Cylindrical'}
      R = vcyl(V);
    case {'spherical', 'Spherical'}
      R = vsph(V);
    otherwise
      disp("Invalid coord. sys. (cartesian, cylindrical or spherical)!");
  end
end

% Converts the vector typed sym to a number
function R = vec2dbl(V, coordinate_system)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  switch coordinate_system
    case {'cartesian', 'Cartesian'}
      R = double(V);
    case {'cylindrical', 'Cylindrical'}
      R = [ double(V(1)) rad2deg(double(V(2))) double(V(3))];
    case {'spherical', 'Spherical'}
      R = [ double(V(1)) rad2deg(double(V(2:3))) ];
    otherwise
      disp("Invalid coord. sys. (cartesian, cylindrical or spherical)!");
  end
end

% Converts sym vector boundaries to numbers
function R = veclims2dbl(L, coordinate_system)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  switch coordinate_system
    case {'cartesian', 'Cartesian'}
      R = double(L);
    case {'cylindrical', 'Cylindrical'}
      R = [ double(L(1:2)) rad2deg(double(L(3:4))) double(L(5:6))];
    case {'spherical', 'Spherical'}
      R = [ double(L(1:2)) rad2deg(double(L(3:6))) ];
    otherwise
      disp("Invalid coord. sys. (cartesian, cylindrical or spherical)!");
  end
end

% Converts values entered in degrees to radians
function R = em_vdeg2rad(V, coordinate_system)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  if all(size(V) == [1 3]) || all(rows(V) > 1 && length(V) == 3)
    switch coordinate_system
      case {'cartesian', 'Cartesian'}
        R = vcart(V);
      case {'cylindrical', 'Cylindrical'}
        R = vcyl(V);
      case {'spherical', 'Spherical'}
        R = vsph(V);
    end
  elseif all(size(V) == [1 6])
    switch coordinate_system
      case {'cartesian', 'Cartesian' }
        R(1:2:6) = vcart(V(1:2:end));
        R(2:2:6) = vcart(V(2:2:end));
      case { 'cylindrical', 'Cylindrical' }
        R(1:2:6) = vcyl(V(1:2:end));
        R(2:2:6) = vcyl(V(2:2:end));
      case { 'spherical', 'Spherical' }
        R(1:2:6) = vsph(V(1:2:end));
        R(2:2:6) = vsph(V(2:2:end));
    end
  else
    error("An invalid vector was entered!");
  end
end

% Returns the vector as (a,b,c) (with degree symbol)
function R = vec2strd(V, coordinate_system)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  switch coordinate_system
    case {'cartesian', 'Cartesian'}
      ft = "(%g,%g,%g)";
    case {'cylindrical', 'Cylindrical'}
      ft = "(%g,%g°,%g)";
    case {'spherical', 'Spherical'}
      ft = "(%g,%g°,%g°)";
  end
  if ~isnumeric(V)
    V = vec2dbl(V, coordinate_system);
  end
  R = sprintf(ft, V);
end

% Returns the vector as +a(âx) +b(ây) +c(âz)
function R = vec2str(V, coordinate_system)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  D = double(V);
  switch coordinate_system
      case { 'cartesian', 'Cartesian' }
        A = { 'âx'; 'ây'; 'âz' };
      case { 'cylindrical', 'Cylindrical' }
        A = { 'âρ'; 'âφ'; 'âz' };
      case { 'spherical', 'Spherical' }
        A = { 'âr'; 'âθ'; 'âφ' };
      otherwise
        error("Coordinate system is not detected!");
    end
    R = "";
    for i = 1:3
      if sign(D(i)) == -1; S = ""; else; S = "+"; end;
      if D(i) ~= 0
        R = [ R S num2str(D(i)) "(" A(i){:} ") " ];
      end
    end
end

% displays the vector appropriately
function R = num2engvec(V, coordinate_system, unit)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  if nargin < 3; unit=""; end
  V = double(V);
  k_min = Inf;
  k_max = -Inf;
  for value = V
    if value ~= 0
      range = floor(log10(value));
      if range < k_min; k_min = range; end
      if range > k_max; k_max = range; end
    end
  end

  if isfinite(k_min) && isfinite(k_max)
    ranges = "fpnum kMGTPE";
    range = 7;
    k2 = floor(k_max/3);
    if k2 >= -5 && k2 <= 6
        dispval = V ./ 10^(k2*3);
        range_t = strtrim(ranges(k2+6));
        R = [vec2str(dispval, coordinate_system), " ", range_t, unit ];
    else
        R = [vec2str(V, coordinate_system), " ", unit];
    end
  else
    R = [vec2str(V, coordinate_system), " ", unit];
  end
end

'%%%%%%%%% WORK AND ENERGY IN THE FIELD OF STATIC ELECTRICITY  %%%%%%%%%';

% Calculates the work done against the electric field 
% if W<0, it is the work done by the electric field
function em_calc_work(E, Q, P, coordinate_system)
if nargin == 0; show_ynsrc_help(); return; end; % show help
if nargin < 4; coordinate_system=vdetcoordsys(E); end;
X = vdetcoordsys(coordinate_system);

if isnumeric(Q)
  Q = sym(num2str(Q));
end

if rows(P) == 2
  A = P(1,:);
  B = P(end,:);
  C = [ A ];

  while ~isequal(A,B)
    if A(1) ~= B(1)
      A(1) = B(1);
      C = [ C; A ];
    elseif A(2) ~= B(2)
      A(2) = B(2);
      C = [ C; A ];
    elseif A(3) ~= B(3)
      A(3) = B(3);
      C = [ C; A ];
    end
  end

  if ~isequal(P,C)
    P = C;
    printf("Attention: Missing steps on the path are completed!\n");
    printf("This is the new path:\n");
    for i = 1:rows(P)-1
      printf("%s -> ", vec2strd(P(i,:), coordinate_system));
    end
    printf("%s\n\n", vec2strd(P(end,:), coordinate_system));
  end
end

W_total = 0;

printf("Q charge = %s\nFrom the formula W = -Q.∫(E)dL;\n\n", num2eng(Q, "C"));

for i = 1 : rows(P)-1
  A = P(i,:);
  B = P(i+1,:);
  dL = [ 1, 1, 1 ];

  switch coordinate_system
    case {'cartesian', 'Cartesian'}
      dL = [ 1, 1, 1 ];
    case {'cylindrical', 'Cylindrical'}
      dL = [ 1, X(1), 1 ];
    case {'spherical', 'Spherical'}
      dL = [ 1, X(1), X(1)*sin(X(2)) ];
  end

  W(1) = -Q * subs( int(E(1)*dL(1),X(1),A(1),B(1)), [X(2) X(3)], [A(2) A(3)]);
  W(2) = -Q * subs( int(E(2)*dL(2),X(2),A(2),B(2)), [X(1) X(3)], [A(1) A(3)]);
  W(3) = -Q * subs( int(E(3)*dL(3),X(3),A(3),B(3)), [X(1) X(2)], [A(1) A(2)]);

  W_sum = double(sum(W));
  W_total = W_total + W_sum;

  printf(["(%d-) By moving charge Q %s -> %s\n     work done (W) = %s\n"],
    i, vec2strd(A,coordinate_system), vec2strd(B, coordinate_system),
    num2eng(W_sum,'J'));
end
printf("\nTotal work done = %s\n\n", num2eng(W_total,'J'));
end

% Same as em_calc_work but in path P angles are in degrees
function em_calc_workd(E, Q, P, coordinate_system)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  if nargin < 4; coordinate_system=vdetcoordsys(E); end;
  switch coordinate_system
    case {'cylindrical', 'Cylindrical'}
      P = vcyl(P);
    case {'spherical', 'Spherical'}
      P = vsph(P);
  end
  em_calc_work(E, Q, P, coordinate_system);
end

% Calculating the energy of the system from the electrical charges given their positions
function em_calc_energy(Q, P)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  printf("Electric charges in system:\n");

  for i = 1:rows(P)
    printf(" %8s charge located at %8s,\n", num2eng(Q(i), "C"),
      vec2strd(P(i,:),'Cartesian'));
  end

  printf("then;\n\n");

  if length(Q) == rows(P)
    W = 0;
    for i = 1 : length(Q)
      S = 0;
      for j = 1 : length(Q)
        if i ~= j
          S = S + 9e9 * Q(j) / norm(P(i,:) - P(j,:));
        end
      end
      W = W + Q(i)/2 * S;
    end
    printf("Total energy of the system = %s\n", num2eng(W, "J"));
  else
    error("The location of each load must be specified!");
  end
end

'%%%%%%%%%%%%%%%%%%%%%%%%% USEFUL TRANSFORMATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%';

% calculate electric field from potential difference
function E = em_v2e(V, coordinate_system)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  if nargin < 2; coordinate_system=vdetcoordsys(V); end;
  X = vdetcoordsys(coordinate_system);
  E = (-1)*gradient_sym(V, X, coordinate_system);
end

% calculate electric flux density from potential difference
function D = em_v2d(V, coordinate_system, epsilon_r=1)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  if nargin < 2; coordinate_system=vdetcoordsys(V); end;
  D = em_e2d(em_v2e(V, coordinate_system), epsilon_r);
end

% finds volume charge density (pV) from electric flux density
function pV = em_d2pv(D, coordinate_system)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  if nargin < 2; coordinate_system=vdetcoordsys(D); end;
  X = vdetcoordsys(coordinate_system);
  pV = divergence_sym(D, X, coordinate_system);
end

% finds the volume charge density (pV) from the electric field
function pV = em_e2pv(E, coordinate_system, epsilon_r=1)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  if nargin < 2; coordinate_system=vdetcoordsys(E); end;
  pV = em_d2pv(em_e2d(E, epsilon_r), coordinate_system);
end

% calculates the electric flux density from the electric field expression
function D = em_e2d(E, epsilon_r = 1)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  syms epsilon_0;
  D = epsilon_0 * epsilon_r * E;
end

% calculate electric field from electric flux density expression
function E = em_d2e(D, epsilon_r = 1)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  syms epsilon_0;
  E = D / (epsilon_0 * epsilon_r);
end

'%%%%%%%%%%%%%%%%%%%% FORCE APPLIED BY THE LOAD %%%%%%%%%%%%%%%%%%%%%%';

% Calculates the force exerted by charge Q1 on charge Q2
function F = em_f12(Q1, A, Q2, B, epsilon_r=1, coordinate_system='Cartesian')
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  if ~strcmp(class(Q1),'sym'); Q1 = sym(num2str(Q1)); end;
  if ~strcmp(class(Q2),'sym'); Q2 = sym(num2str(Q2)); end;
  syms epsilon_0;
  R = vec2sym(B - A, coordinate_system);
  k = 1/(4 * sym('pi') * epsilon_r * epsilon_0);
  F = em_subs(k*Q1*Q2*R/norm(R)^3);
  if length(symvar(F)) == 0; F = double(F); end;
end

% The total force of the charges Qn at points B on the charge Q at point A
function F = em_fnet(Q1, A, Qn, Bn, epsilon_r=1, coordinate_system='Cartesian')
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  F = 0;
  for i = 1 : length(Qn)
    Q2 = Qn(i);
    B = Bn(i,:);
    F = F + em_f12(Q2, B, Q1, A, epsilon_r, coordinate_system);
  end
end

'%%%%% POINT, INFINITE LINEAR AND SURFACE LOAD DISTRIBUTIONS %%%%%%%%%';

% electric field created by a point charge at a point
function E = em_e_point(Q, A, B, epsilon_r=1, coordinate_system='Cartesian')
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  syms epsilon_0;
  if ~strcmp(class(Q),'sym'); Q = sym(num2str(Q)); end;
  R = vec2sym(B - A, coordinate_system);
  k = 1/(4 * sym('pi') * epsilon_r * epsilon_0);
  E = em_subs(k*Q*R/(norm(R)^3));
  if length(symvar(E)) == 0; E = double(E); end;
end

% electric field of infinite linear charge density
% The component with Inf is assigned to the component at the position of the point.
function E = em_e_linear(pL, A, B, epsilon_r=1, coordinate_system='Cartesian')
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  syms epsilon_0;
  S = find(~isfinite(A));
  A(S) = B(S);
  if ~strcmp(class(pL),'sym'); pL = sym(num2str(pL)); end;
  R = vec2sym(B - A, coordinate_system);
  k = 1/(2 * sym('pi') * epsilon_r * epsilon_0);
  E = em_subs(k * pL * R / (norm(R)^2));
  if length(symvar(E)) == 0; E = double(E); end;
end

% electric field of infinite planar charge density
% Components with Inf are assigned to the component at the position of the dot.
function E = em_e_planar(pS,A,B, epsilon_r=1, coordinate_system='Cartesian')
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  if ~strcmp(class(pS),'sym'); pS = sym(num2str(pS)); end;
  syms epsilon_0;
  S = find(~isfinite(A));
  A(S) = B(S);
  R = vec2sym(B-A, coordinate_system);
  k = 1 / (2 * epsilon_r * epsilon_0);
  E = em_subs(k * pS * sign(R));
  if length(symvar(E)) == 0; E = double(E); end;
end

'%%%%%%%%%%%%%%%%%%%%%%%%%%% SURFACE INTEGRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%';

% Calculates the line integral of vector V
function [fS,fT,S] = em_intline(V, L, coord_sys)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  if nargin < 3; coord_sys=vdetcoordsys(V); end;
  X = vdetcoordsys(coord_sys);
  S = [ 0 0 0 0 0 0 ];

  for i = 1:2:5
    S(i) = all(L(i+1) > L(i));
    S(i+1) = all(L(i+1) < L(i));
  end

  switch coord_sys
    case { 'cartesian', 'Cartesian' }
      kS = { 1, 1, 1 };
    case { 'cylindrical', 'Cylindrical' }
      kS = { 1, X(1), 1 };
    case { 'spherical', 'Spherical' }
      kS = { 1, X(1), X(1)*sin(X(2)) };
  end

  clear fS;
  fS = sym([]);
  fS(6) = 0;

  K = [ 1 1 2 2 3 3 ];

  dA = [ X(1) L(1:2) ];
  dB = [ X(2) L(3:4) ];
  dC = [ X(3) L(5:6) ];

  dX = [ dA; dB; dC ];

  sX = [ X(2) L(3) X(3) L(5); X(1) L(1) X(3) L(5); X(1) L(1) X(2) L(3) ];

  for i = 1:6
    k = K(i);
    if mod(i,2) == 1; sN = 1; else; sN = -1; end;
    if S(i)
      fS(i) = int1subs([ sN*V(k)*kS{k}, dX(k,:), sX(k,:) ]);
    end
  end
  fT = sum(fS(find(~isnan(fS))));
end

% Same as "em_intline" but limit values are in degrees
function [fS, fT, S] = em_intlined(V, L, coord_sys)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  if nargin < 3; coord_sys=vdetcoordsys(V); end;
  [fS, fT, S] = em_intline(V, em_vdeg2rad(L, coord_sys), coord_sys);
end


% Calculates the surface integral of the vector V
function [fS,fT,S] = em_intsurf(V, L, coord_sys)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  if nargin < 3; coord_sys=vdetcoordsys(V); end;
  X = vdetcoordsys(coord_sys);
  S = [ 0 0 0 0 0 0 ];
  switch coord_sys
    case { 'cartesian', 'Cartesian' }
      S(1) = all(L(3)~=L(4) && L(5) ~= L(6));
      S(2) = all(S(1) && L(1)~=L(2));
      S(3) = all(L(1)~=L(2) && L(5) ~= L(6));
      S(4) = all(S(3) && L(3)~=L(4));
      S(5) = all(L(1)~=L(2) && L(3) ~= L(4));
      S(6) = all(S(5) && L(5)~=L(6));
      kS = { 1, 1, 1 };
    case { 'cylindrical', 'Cylindrical' }
      S(1) = all(L(3)~=L(4) && L(5)~=L(6));
      S(2) = all(S(1) && L(1)~=L(2) && L(1)>0);
      S(3) = all(L(1)~=L(2) && L(5)~=L(6) && abs(L(4)-L(3))~=sym('2*pi'));
      S(4) = all(S(3) && L(1)~=L(2) && L(3)~=L(4));
      S(5) = all(L(1)~=L(2) && L(3)~=L(4));
      S(6) = all(S(5) && L(5)~=L(6));
      kS = { X(1), 1, X(2) };
    case { 'spherical', 'Spherical' }
      S(1) = all(L(3)~=L(4) && L(5)~=L(6));
      S(2) = all(S(1) && L(1)~=L(2) && L(1) > 0);
      S(3) = all(L(1)~=L(2) && L(5)~=L(6) && L(4) < sym('pi'));
      S(4) = all(L(1)~=L(2) && L(5)~=L(6) && L(3) > 0);
      S(5) = all(L(1)~=L(2) && L(3)~=L(4) && abs(L(6)-L(5))~=sym('2*pi'));
      S(6) = all(S(5) && L(5)~=L(6));
      kS = { X(1)^2*sin(X(2)), X(1)*sin(X(2)), X(1) };
  end

  clear fS;
  fS = sym([]);
  fS(6) = 0;

  K = [ 1 1 2 2 3 3 ];

  dA = [ X(1) L(1:2) ];
  dB = [ X(2) L(3:4) ];
  dC = [ X(3) L(5:6) ];

  dX = [ dB dC; dA dC; dA dB ];
  dC = [ 2 1 4 3 6 5 ];

  for i = 1:6
    k = K(i);
    if mod(i,2) == 1; sN = 1; else; sN = -1; end;
    if S(i)
      fS(i) = int2subs([ sN*V(k)*kS{k},  dX(k,:), X(k), L(dC(i)) ]);
    end
  end
  fT = sum(fS(find(~isnan(fS))));
end

% Same as "em_intsurf" but with limit values in degrees
function [fS, fT, S] = em_intsurfd(V, L, coord_sys)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  if nargin < 3; coord_sys=vdetcoordsys(V); end;
  [fS, fT, S] = em_intsurf(V, em_vdeg2rad(L, coord_sys), coord_sys);
end

'%%%%%%%%%%%%%%%%% GAUSS LAW AND DIVERGENCE THEOREM %%%%%%%%%%%%%%%%%%%%';

% angles and boundaries must be specified in radians (see em_calc_gaussd)
function [fS,psi] = em_calc_gauss(D, N=[], L=[], L2=[], coord_sys)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  if nargin < 5; coord_sys=vdetcoordsys(D); end;
  X = vdetcoordsys(coord_sys);

  disp("Electric flux density (D): "), disp(D), disp("");

  pV = em_d2pv(D, coord_sys);
  disp("a) Volume charge density pV (∇.D): "), disp(pV);

  if ~isempty(N)
    N_pV = double(em_subs(pV, vdetcoordsys(coord_sys), N));
    printf("\nb) Volume charge density at (%g,%g,%g) (pV) = %s\n\n",
        vec2dbl(N, coord_sys), num2eng(N_pV, "C/m³"));
  end

  switch coord_sys
    case { 'cartesian', 'Cartesian' }
      T = 'kutu';
      V = { 'x', 'y', 'z' };
      dS = { 'dydz', 'dxdx', 'dxdy' };
      dV = pV;

    case { 'cylindrical', 'Cylindrical' }
      T = 'silindir';
      V = { 'ρ', 'φ', 'z' };
      dS = { 'ρdφdz', 'dρdz', 'ρdρdφ' };
      dV = pV * X(1);

    case { 'spherical', 'Spherical' }
      T = 'küre';
      V = { 'r', 'θ', 'φ' };
      dS = { 'r²sin(θ)dθdφ', 'rsin(θ)drdφ', 'rdrdθ' };
      dV = pV * X(1)^2 * sin(X(2));
  end

  if ~isempty(L)
    printf(["c) For %g<=" V{1} "<=%g, %g<=" V{2} "<=%g, %g<=" V{3} ...
      "<=%g bounded %s\n"], veclims2dbl(L, coord_sys), T);

    [fS, fT, S] = em_intsurf(D, L);

    printf("   There are %d surface(s) in total.\n\n", sum(S));

    psi_c = 1;

    K = [ 1 1 2 2 3 3 ];

    for i = 1:6
      k = K(i);
      if S(i)
        if mod(i,2) == 1; sP = ""; sS=" "; sN = 1;
        else; sP = "-"; sS="-"; sN = -1; end;
        printf([ sS "â" V(k){:} " directed electric flux ∫(D)dS,dS="  dS(k){:} ...
          "(%sâ%s) = %s\n   Ψ%d = %s [C]\n\n"], sP, V(k){:},
           num2eng(fS(i),"C"), psi_c, char(fS(i)));
        psi_c = psi_c + 1;
      end
    end

    printf("\nLeft side of the equation ∮(D)dS = Ψ = Q = %s\n", num2eng(fT, "C"));
    disp(simplify(fT));

    psi = int3([dV, X(1) L(1:2), X(2) L(3:4), X(3) L(5:6)]);

    printf("\nRight side of the equation ∫(pV)dV = Ψ = Q = %s\n", num2eng(psi, "C"));
    disp(simplify(psi));

    disp("\nSince left is equal to right, Divergence theorem is satisfied.\n");
  else
    fS = 0;
    psi = 0;
  end

  if ~isempty(L2)
    if L2(1)==L2(2) || L2(3)==L2(4) || L2(5)==L2(6)
      [fS2, fT2, S2] = em_intsurf(D, L2);

      printf(["e) Flux leaving the bounded surface %g<=" V{1} "<=%g, %g<=" V{2} "<=%g, %g<=" V{3} ...
      "<=%g\n"], veclims2dbl(L2, coord_sys));

      fL2 = fS2(find(S2==1));

      printf("   Ψ = Q = %s\n\n", num2eng(fL2, "C"));
    else
      disp("Caution: L2 does not specify a valid surface boundary!");
      disp("Therefore the flux emerging from the confined surface was not calculated.");
    end
  end
  disp("");
end

% Like "em_calc_gauss" but angles are entered in degrees
function [fS, psi] = em_calc_gaussd(D, Nd=[], Ld=[], L2d=[], coord_sys)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  if nargin < 5; coord_sys=vdetcoordsys(D); end;
  if ~isempty(Nd); N = em_vdeg2rad(Nd, coord_sys); end;
  if ~isempty(Ld); L = em_vdeg2rad(Ld, coord_sys); end;
  if ~isempty(L2d); L2 = em_vdeg2rad(L2d, coord_sys); end;

  if ~exist('N', 'var'); N = []; end;
  if ~exist('L', 'var'); L = []; end;
  if ~exist('L2', 'var'); L2 = []; end;

  [fS, psi] = em_calc_gauss(D, N, L, L2, coord_sys);
end

'%%%%%%%%%%%%%%%%%%% CURRENT DENSITY, MATERIAL ETC. %%%%%%%%%%%%%%%%%%%%';

% Finds the index in the direction vector S of the specified direction
function R = find_surface_index(direction, coordinate_system)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  dir = char(direction);
  neg = startsWith(dir,'-');
  if neg; dir = dir(2:end); end;
  switch dir
    case { 'x', 'rho', 'r' }
      R = 1;
    case { 'y', 'theta' }
      R = 3;
    case { 'z' }
      R = 5;
    case { 'phi' }
      switch coordinate_system
        case { 'cylindrical', 'Cylindrical' }
          R = 3;
        case { 'spherical', 'Spherical' }
          R = 5;
      end
  end
  if neg; R = R + 1; end;
end

function I = em_i_qt(Q,t)
  I = Q/t;
end

% Calculates the current coming out of the surface from the current density
function I = em_j2i(J,S,surface_direction=[], coordinate_system)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  if nargin < 4; coordinate_system = vdetcoordsys(J); end;
  if ~isempty(surface_direction)
    surface_index = find_surface_index(surface_direction, coordinate_system);
  else
    surface_index = -1;
  end
  I = em_intsurf(J, S);
  if surface_index > 0; I = I(surface_index); end;
  if length(symvar(I)) == 0; I = double(I); end;
end

% Like em_j2i but the boundaries of the shape are in degrees
function I = em_j2id(J,S,surface_direction=[], coordinate_system)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
   if nargin < 4; coordinate_system = vdetcoordsys(J); end;
   I = em_j2id(E, em_vdeg2rad(S), surface_direction, coordinate_system);
end

% Calculates electric field and current from sigma
function I = em_e2i(E, S, sigma, surface_direction=[], coordinate_system)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  if nargin < 5; coordinate_system = vdetcoordsys(E); end;
  if ~isempty(surface_direction)
    surface_index = find_surface_index(surface_direction, coordinate_system);
  else
    surface_index = -1;
  end
  I = em_intsurf(sigma*E, S);
  if surface_index > 0; I = I(surface_index); end;
  if length(symvar(I)) == 0; I = double(I); end;
end

% Like "em_e2i" but the boundaries of the shape will be degrees
function I = em_e2id(E, S, sigma, surface_direction=[], coordinate_system)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  if nargin < 5; coordinate_system = vdetcoordsys(E); end;
  I = em_e2i(E, em_vdeg2rad(S), sigma, surface_direction, coordinate_system);
end

'%%%%%%%%%%%%%%%%%%%%%%%%%    STOKES THEOREM    %%%%%%%%%%%%%%%%%%%%%%%%';

% H (cartesian) formed at point C by the A->B path passing I current
function H = em_i2h(I, A, B, C)
  if nargin == 0; show_ynsrc_help(); return; end; % show help

  V1 = C - A;
  V2 = B - A;
  V3 = C - B;

  if all(isfinite(A)) && all(isfinite(B))

    % finite (for limited conductor)

    cos_beta = dot(V1,V2) / ( norm(V1) * norm(V2) );
    cos_alpha_1 = -cos_beta;
    cos_alpha_2 = -dot(V2,V3) / ( norm(V2) * norm(V3) );

    _rho = norm(V3) * sin(acos(cos_alpha_2));
    _al = V2/norm(V2);
    K = B - (norm(V3) * cos_alpha_2) * (-_al);
    _arho = (C-K)/_rho;
  else
    % for infinite or semi-infinite conductor

    K = V1;
    K(find(~isfinite(K))) = 0;
    _rho = norm(K);
    _arho = K/norm(K);

    _al = double(sign(V2) .* ~isfinite(V2));

    if ~all(isfinite(A)) && ~all(isfinite(B)) % full infinity
      A_inf = A(~isfinite(A));
      B_inf = B(~isfinite(B));
      if A < B
        cos_alpha_1 = -1;
        cos_alpha_2 = 1;
      else
        cos_alpha_1 = 1;
        cos_alpha_2 = -1;
      end

    elseif ~all(isfinite(A)) % half infinite Inf(A) -> B
      index = find(~isfinite(A));
      if B(index) == C(index)
        cos_alpha_1 = 1;
        cos_alpha_2 = 0;
      else
        disp("CAUTION: Not checked, this result may be inaccurate (!)");
        R2 = V2;
        R2(~isfinite(R2)) = 1;
        R2 = R2/norm(R2);
        R1 = V1/norm(V1);
        cos_beta = dot(R1,R2) / ( norm(R1) * norm(R2) );
        cos_alpha_1 = 1;
        cos_alpha_2 = -cos_beta;
        _rho = norm(V1) * sin(acos(cos_alpha_1));
        K = A - (norm(V1) * cos_beta) * (-_al);
        _arho = (C-K)/_rho;
      end

    elseif ~all(isfinite(B)) % half infinite A -> Inf(B)
      index = find(~isfinite(B));
      if A(index) == C(index)
        cos_alpha_1 = 0;
        cos_alpha_2 = 1;
      else
        R2 = V2;
        R2(~isfinite(R2)) = 1;
        R2 = R2/norm(R2);
        R1 = V1/norm(V1);
        cos_beta = dot(R1,R2) / ( norm(R1) * norm(R2) );
        cos_alpha_1 = -cos_beta;
        cos_alpha_2 = 1;
        _rho = norm(V1) * sin(acos(cos_alpha_1));
        K = A - (norm(V1) * cos_beta) * (-_al);
        _arho = (C-K)/_rho;
      end
    end
  end

  _aphi = cross( _al, _arho );

  H = (I / (4 * pi * _rho) * (cos_alpha_2 - cos_alpha_1)) * _aphi;

  if 0 %debug
  fprintf("\n\tcos_a1=%g, cos_a2=%g", cos_alpha_1, cos_alpha_2);
  fprintf("\n\ta1=%g, a2=%g, rho=%g\n\tâp=%s, âl=%s\n\tâphi=%s\n", ...
    acosd(cos_alpha_1), acosd(cos_alpha_2), _rho, ...
    vec2strd(_arho,'Cartesian'), vec2strd(_al,'Cartesian'), ...
    vec2strd(_aphi,'Cartesian'));
  end
end

% H formed at point C by the current I flowing along the path P
function H = em_calc_h(I, P, C)
  if nargin == 0; show_ynsrc_help(); return; end; % show help
  
  if size(P,1) < 2
    disp("Invalid path, number of lines must be greater than 2!");
  else
    H = 0;
    fprintf("%g A current flowing in the conductor;\n\n", I);
    for i = 1 : size(P,1) - 1
      A = P(i,:);
      B = P(i+1,:);
      H = H + em_i2h(I, A, B, C);
      fprintf(" %s->%s path\n H%s = %s\n\n",...
        vec2strd(A,'Cartesian'), vec2strd(B,'Cartesian'),...
        vec2strd(C,'Cartesian'),...
        num2engvec(em_i2h(I, A, B, C), 'Cartesian', 'A/m'));
    end
    fprintf("Total H%s = %s\n\n", vec2strd(C,'Cartesian'),...
      num2engvec(H,'Cartesian', 'A/m'));
  end
end
