% DC and AC Analysis of Common Emitter
% Emitter Degenerated BJT - NPN
% Amplifier Circuit 
% And Calculating Voltage Gain

pkg load matgeom;
clear; clc;

%---------- Transistor Parameters -----------%

hFE = 200;   % Current Gain (β)
Vt = 26e-3;  % 26 mV
VBE = 0.7; % V

%----------- Circuit Parameters ----------------%

VCC = 12;    % V

% Resistors

R1 = 330e3;  % bias resistor 1
R2 = 33e3;   % bias resistor 2
RC = 4.7e3;  % collector resistor
RE = 470;    % emitter resistor
RL = 100e3;  % load resistor

% Capacitors

CB = 100e-9;  % 100 nF
CL = 100e-9;  % 100 nF

%------------------ DC Analysis -------------------%

parallel_res = @(z) 1 / sum(1./z); % parallel resistors

Vth = (R2/(R1+R2)) * VCC;      % Thevenin Eq. Voltage
Rth = parallel_res([R1, R2]);  % Thevenin Eq. Resistance (R1//R2)

IB = (Vth-VBE) / (Rth + (hFE+1)*RE); % Base Current
IC = hFE * IB;                       % Collector Current
IE = IC + IB;                        % Emitter Current

VC = VCC - IC * RC;  % Collector Voltage
VE = RE * IE;        % Emitter Voltage
VB = VBE + VE;       % Base Voltage

VCE = VC - VE;       % Collector-Emitter Voltage

%-------------------- AC Analysis -------------------%

gm = IC/Vt;     % Transconductance
rPi = Vt/IB;    % rπ Resistance

% Voltage Gain
Av = -parallel_res([RC RL]) / ((1/gm) + RE);

% Input Impedance (rπ+(β+1)RE)//R1//R2
Rin = parallel_res([(rPi+(hFE+1)*RE), R1, R2]);

% Output Impedance (RC//RL)
Rout = parallel_res([RC RL]);

%%%%%%%%%%%%%%%  END OF CALCULATIONS  %%%%%%%%%%%%%%%%

%------------ Graph Drawing Functions ------------%

function draw_power(x, y, angle=90, label="")
    phi = sin(angle*pi/180);
    theta = cos(angle*pi/180);
    x2 = x + 5 * theta;
    y2 = y + 5 * phi;
    arrow = drawArrow([ x y x2 y2], 1, 1, 1);
    set(arrow.body, 'color', 'black');
    set(arrow.wing, 'color', 'black');
    if length(label) > 0
        tx = x + 7.5 * (theta) - 1;
        ty = y + 7 * phi;
        text(tx, ty, label);
    end
end

function draw_component(coordinates, x, y, angle=0)
    R = [
        cosd(angle) -sind(angle);
        sind(angle) cosd(angle)
    ];
    
    rot_points = coordinates * R;
    
    [rows, columns] = size(rot_points);
    
    rot_points(:,1) = rot_points(:,1) + x;
    rot_points(:,2) = rot_points(:,2) + y;
    
    points = [
        rot_points(1:2:rows,:) rot_points(2:2:rows,:)
    ];
    
    drawEdge(points, 'color', 'black');
end

function draw_npn(x, y, angle=0)
    coordinates = [
        [-6 0]; [-1.5 0];
        [-1.5 1.6]; [-1.5 -1.6];
        [-1.5 1.5]; [0 3];
        [0 3]; [0 5];
        [-1.5 -1.5]; [0 -3];
        [-0.5 -2.5]; [0 -1.5];
        [-0.5 -2.5]; [-1.8 -2];
        [0 -3]; [0 -5];
    ];
    
    drawCircle(x, y, 3, 'color', 'black');
    draw_component(coordinates, x, y, angle);
    
    text(x-4.5, y+4.5, "NPN");
    text(x, y+2, "C");
    text(x-4.5, y+1, "B");
    text(x+0.5, y-1.5, "E");
end

function [engstr] = num2eng(value, unit="")        
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

function draw_resistor(x, y, angle=0, name="", value=0)
    coordinates = [
        [-5 0]; [-3 0];
        [-3 0]; [-2 1];
        [-2 1]; [-1 -1];
        [-1 -1]; [0 1];
        [0 1]; [1 -1];
        [1 -1]; [2 0]; 
        [2 0]; [5 0]
    ];
    
    draw_component(coordinates, x, y, angle);
    
    label = num2eng(value, "Ω");
    
    if angle == 0 || angle == 180
        if isempty(name)
            text(x-1, y-3.5, label);
        else
            text(x-1, y-2.5, name);
            text(x-1, y-4.5, label);
        end
    else
        if isempty(name)
            text(x+2, y, label);
        else
            text(x+2, y+1, name);
            text(x+2, y-1, label);
        end
    end
end

function draw_capacitor(x, y, angle=0, name="", value=0)
    coordinates = [
        [-0.5 -2]; [-0.5 2];
        [0.5 -2]; [0.5 2];
        [-0.5 0]; [-5 0];
        [0.5 0]; [5 0]
    ];
    
    draw_component(coordinates, x, y, angle);
    
    label = num2eng(value, "F");
    
    if length(label) > 0
        if angle == 0 || angle == 180
            if isempty(name)
                text(x-1, y-4.5, label)
            else
                text(x-1, y-3.5, name)
                text(x-1, y-5.5, label)
            end
        else
            if isempty(name)
                text(x+3, y, label)
            else
                text(x+3, y-1, label)
                text(x+3, y+1, label)
            end
        end
    end
end

function draw_wire(points)
    drawEdge(points, 'color', 'black');
end

function draw_info(x, y, name, value, unit, color='black')
    label = sprintf("%s = %s", name, num2eng(value, unit));
    text(x, y, label, 'color', color);
end

%------------- Drawing Circuit ---------------%

figure, hold on, axis equal;
axis([5 70 0 50]);
axis off;

title("Common Emitter (Degenerated) Calculation");

x_left = 20;     % left edge of circuit
y_top = 40;      % top edge of circuit

y_center = y_top - 15;
y_bottom = y_center - 15;

y_res_top = y_top - 5;  
y_res_center = y_res_top - 10;
y_res_bottom = y_res_center - 10;

x_center = x_left + 10;
x_right = x_center + 15;

draw_power(x_right, y_top, 90, sprintf("VCC = %d V", VCC));
draw_power(x_right, y_bottom, 270, "GND");
draw_power(x_left, y_center, 90, "Vin");
draw_power(x_right+20, y_res_top-5, 90, "Vout");

draw_npn(x_right, y_center);

draw_resistor(x_right, y_res_top, 90, "RC", RC);
draw_resistor(x_right, y_res_bottom, 90, "RE", RE);
draw_resistor(x_center, y_res_top, 90, "R1", R1);
draw_resistor(x_center, y_res_bottom, 90, "R2", R2);
draw_resistor(x_right + 20, y_bottom + 5, 90, "RL", RL);

draw_capacitor(x_center - 5, y_center, 0, "CB", CB);
draw_capacitor(x_right + 10, y_res_top - 5, 0, "CL", CL);

draw_wire([
    [x_center y_res_top-5 x_center y_res_bottom+5];
    [x_center y_top x_right y_top];
    [x_center y_bottom x_right y_bottom];
    [x_center y_center x_right-6 y_center ];
    [x_right y_bottom x_right+20 y_bottom];
    [x_right y_res_top-5 x_right+5 y_res_top-5];
    [x_right+15 y_res_top-5 x_right+20 y_res_top-5];
    [x_right+20 y_res_top-5 x_right+20 y_res_center-5]
]);

draw_info(5, 48, 'hFE (β)', hFE, '', 'blue');
draw_info(5, 46, 'Vt', Vt, 'V', 'blue');
draw_info(5, 44, 'VBE', VBE, 'V', 'blue');

draw_info(5, 40, 'IB', IB, 'A', [0 .5 .2]);
draw_info(5, 38, 'IC', IC, 'A', [0 .5 .2]);
draw_info(5, 36, 'IE', IE, 'A', [0 .5 .2]);

draw_info(5, 32, 'VC', VC, 'V');
draw_info(5, 30, 'VB', VB, 'V');
draw_info(5, 28, 'VE', VE, 'V');

draw_info(5, 25, 'VCE', VCE, 'V', [0 .4 0]);

draw_info(5, 22, 'Rin', Rin, 'Ω');
draw_info(5, 20, 'Rout', Rout, 'Ω');
draw_info(5, 18, 'Av', Av, 'V/V', [.7 .2 .7]);

draw_info(5, 14, 'Vth', Vth, 'V');
draw_info(5, 12, 'Rth', Rth, 'Ω');

draw_info(5, 8, 'gm', gm, '℧', [.8 0 0]);
draw_info(5, 6, 'rπ', rPi, 'Ω', [.8 0 0]);
