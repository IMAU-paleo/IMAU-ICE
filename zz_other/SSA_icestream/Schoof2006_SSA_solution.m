function [u,tau_c] = Schoof2006_SSA_solution( tantheta, h0, B, L, m, y)

% parameters_module
ice_density           =  910.0;     % Ice density [kg m^-3]
grav                  = 9.81;       % Acceleration of gravity [m s^-2]

% The experimental configuraton (Bueler & Brown 2009)
% theta                 = atan(0.001); % Bedrock slope in the x-direction (1 in 1000)
% h0                    = 2000;        % Ice thickness
% B                     = 3.7e8;       % Ice hardness
% L                     = 40000;       % Half-width of weak till
% m                     = 10;          % Ice stream width scaling factor (Schoof 2006 shows solutions for m=1, m=10, m=20)

% Calculate the gravitational driving stress f
f = ice_density * grav * h0 * tantheta;

% Calculate the "ice stream half-width" W
W = L * (m+1)^(1/m);

% Calculate the till yield stress across the stream
% y = linspace(-3*L,3*L,10000);
tau_c = f * abs(y/L).^m;

% Calculate the analytical solution for u
ua = -2 * f^3 * L^4 / (B^3 * h0^3);
ub = ( 1 / 4                  ) * (   (y/L).^     4  - (m+1)^(   4/m) );
uc = (-3 / ((m+1)   * (  m+4))) * (abs(y/L).^(  m+4) - (m+1)^(1+(4/m)));
ud = ( 3 / ((m+1)^2 * (2*m+4))) * (abs(y/L).^(2*m+4) - (m+1)^(2+(4/m)));
ue = (-1 / ((m+1)^3 * (3*m+4))) * (abs(y/L).^(3*m+4) - (m+1)^(3+(4/m)));

u = ua * (ub + uc + ud + ue);
u( abs(y)>W) = 0;

end