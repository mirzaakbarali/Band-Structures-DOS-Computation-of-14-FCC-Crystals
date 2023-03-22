%% START DEFAULT COMMANDS OF *.m SCRIPT
clear all;       % Clear all variables/functions in memory
close all force; % Close all figures already opened
clc;             % Clear screen in the command window

iota = sqrt(-1);

%% THE SEVEN EXACT DEFINING CONSTANTS OF THE SI UNIT SYSTEM (2019 UPDATE)
hyperfine = 9192631770;       % Hyperfine transition frequency of Cs-133 [Hz]
celeritas = 299792458;        % Speed of light in vacuum [m/s]
Planck    = 6.62607015E-34;   % Planck's constant [Js]
qel       = 1.602176634E-19;  % Elementary charge (Absolute value of electron charge) [C]
kB        = 1.380640E-23;     % Boltzmann's constant [J/K]
Avogadro  = 6.02214076E-23;   % Avogadro's constant [1/mole]
kcd       = 683;              % Luminous efficacy of 540 THz radiation [candela=lumen/Watt]
                              % Green light at 555,016 nm = maximum possible luminous efficacy.
                              % Originally = peak sensitivity of "average" human eye.

%% PHYSICAL CONSTANTS OF ELECTROMAGNETISM
epsilon0  = 8.8541878128E-12; % Electrical constant (vacuum dielectric permittivity) [F/m]
mu0       = 1.25663706212E-6; % Magnetic constant (vacuum magnetic  permeability) [N/A^2]
                              % close to 4*pi*1.E-7 =1.25663706143E-6
Klitzing  = 25812.80745;      % von Klitzing's constant = Planck/qel^2
                              % respecting significant digits [Ohm]

%% UNITS USED IN ATOMIC, MOLECULAR AND SOLID STATE PHYSICS
hbar      = 1.054571817E-34;                    % Reduced Planck's constant = Planck/(2*pi)
                                                % respecting significant digits [Js]
Angstroem = 1.0E-10;                            % Angström [m]
amu       = 1.66053906660E-27;                  % Atomic mass unit [kg] = 1 Dalton
                                                % = mass of Carbon_12 atom / 12
elm       = 9.1093837015E-31;                   % Electron mass [kg]
nem       = 1.67492749804E-27;                  % Neutron mass [kg]
prm       = 1.67262192369E-27;                  % Proton mass [kg]
elecint   = qel^2/(4*pi*epsilon0);              % Scale of electron-electron interaction
Bohr      = hbar^2/(elm*elecint)/Angstroem;     % Bohr radius [Angström]
Rydberg   = (elm/(2*hbar^2)) * (elecint)^2/qel; % Rydberg [eV]
Hartree   = 2*Rydberg;                          % Hartree [eV]

%% Units used in calculations of Band Structures

recipunit = 1.0e+10; % Inverse Angström
ekin_fact = ((hbar*recipunit)^2/(2 * elm))/qel;

%% WARNING
% The rest of this file can only be executed when sourced in a *.m file
% that includes plotting commands.

%% DEFAULT PLOT CONFIGURATION
% get(groot,'factory')           % List factory-defined plot configurations
% get(groot,'factoryObjectType') % List factory-defined properties for a specific object.
                                 % Examples of 'ObjectType' are: 'Axes', 'Figure', 'Image',
                                 % 'Line','Surface', 'Text', 'ui' (user interface), etc.

% groot is the "handle index" (identifier) of the object
% that is "parent" of all plot objects of the session.

%% THICKER LINES & CHARACTERS FOR VIDEOPROJECTION IN CLASSROOM
set(groot,'defaultLineLineWidth',1);  % Function "set" is not case sensitive.
                                      % Upper cases for
                                      % readability.
set(groot,'defaultAxesFontSize',12);
set(groot,'defaultAxesFontWeight','bold');
set(groot,'defaultAxesLineWidth',2);
set(groot,'defaultAxesXaxisLocation','bottom'); % Location of abcissae ticks & labels
                                                % Possible values:'top','bottom',
                                                % left', 'right', 'origin'.
                                                % 'origin' puts abcissae labels
                                                % on y=0 line

%% DEFAULT COLOR ORDER
% default_colors=get(gca,'colororder'); % gca = "get current axes" = axe system identifier
% default_colors = default order for coloring curves in last versions of Matlab/Octave:
% (1) Tropical blue, (2) Deep orange, (3) Deep yellow,
% (4) Violet,        (5) Grass green, (6) Azure blue,   (7) Sienna.

%% DEFINING YOUR OWN COLOR ORDER

mycolors=[ % Color order of successive curves must be defined by a 7x3 matrix
0    0    0      % Black
1    0    0      % Red
0    0    1      % Blue
0    0.5  0      % Dark Green
0.9  0.5  0.1    % Orange
0    0.75 0.75   % Turquoise
0.5  0.5  0.5    % Grey
];

%set(groot,'defaultAxesColorOrder',mycolors);

%% TRACING HORIZONTAL AND VERTICAL LINES ACROSS THE PLOT WINDOW

% yline(val) traces horizontal line y(x)=val using current xlim
yline = @(yval, varargin) line(xlim, [yval yval], varargin{:});

% xline(val) traces vertical   line at x=val using current ylim
xline = @(xval, varargin) line([xval xval], ylim, varargin{:});

%% TRACING X=0 AND Y=0 AXES
function z = plotzeros() % Warning: Octave implementation only!
                         % Matlab requires to define this
                         % function in an independent file plotzeros.m
                         % or at the end of the main script
    xline = @(xval, varargin) line([xval xval], ylim, varargin{:});
    yline = @(yval, varargin) line(xlim, [yval yval], varargin{:});
    xline(0,'color',[0 0 0],'linewidth',0.5,'linestyle','-');
    yline(0,'color',[0 0 0],'linewidth',0.5,'linestyle','-');
end

%% END DEFAULT COMMANDS OF *.m SCRIPT

