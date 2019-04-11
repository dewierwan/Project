function [y,S_0,S_E,E_price,G_c,Load,T_in,T_out,v] = calcs()

Dt = 1;               % Size of time step [hours]
n = 10;                 % Number of EVs
J = 24/Dt;              % Number of time steps
etaP = 0.99;            % Charging efficiency
Q_s = 70;               % Battery capacity (70 kWh)
S_max = 0.85;           % Maximum Battery SoC
S_min = 0.15;           % Minimum Battery SoC
%w1 = 2e0;          	% Cost weight
%w2 = 1.5e0;             % S - 0.5 weight
%w3 = 3e-2;              % P^2 weight
Pmax = 7;              % Max charging/discharging power
%alpha_s = 1e0;         % Scaling Factor
Z = etaP * Dt / Q_s;    % Constant 

G_c = zeros(J,1);       % Working point irradiance 
Load = zeros(J,1);      % Basic Load at k

L = 30;               % Base Load

%% y 

y(1) = Z;  
y(2) = S_max; 
y(3) = S_min; 
y(4) = n;
y(5) = J;
y(6) = Dt;
%y(7) = w1;     % Cost
%y(8) = w2;     % S - 0.5
%y(9) = w3;     % P^2
y(10) = Pmax;
%y(10) = alpha_s;

%% Wind Speed

v(1:1/Dt) = 21.7*0.514;     %met office wow valley 18 march 2018
v(1+1/Dt:2/Dt) = 20*0.514;
v(1+2/Dt:3/Dt) = 22.6*0.514;
v(1+3/Dt:4/Dt) = 21.7*0.514;
v(1+4/Dt:5/Dt) = 20.9*0.514;
v(1+5/Dt:6/Dt) = 18.3*0.514;
v(1+6/Dt:7/Dt) = 22.6*0.514;
v(1+7/Dt:8/Dt) = 20.0*0.514;
v(1+8/Dt:9/Dt) = 20*0.514;
v(1+9/Dt:10/Dt) = 19.1*0.514;
v(1+10/Dt:11/Dt) = 19.1*0.514;
v(1+11/Dt:12/Dt) = 14.8*0.514;
v(1+12/Dt:13/Dt) = 13.9*0.514;
v(1+13/Dt:14/Dt) = 11.3*0.514;
v(1+14/Dt:15/Dt) = 12.5*0.514;
v(1+15/Dt:16/Dt) = 14.8*0.514;
v(1+16/Dt:17/Dt) = 13*0.514;
v(1+17/Dt:18/Dt) = 8.7*0.514;
v(1+18/Dt:19/Dt) = 7.5*0.514;
v(1+19/Dt:20/Dt) = 6.1*0.514;
v(1+20/Dt:21/Dt) = 8.7*0.514;
v(1+21/Dt:22/Dt) = 7.8*0.514;
v(1+22/Dt:23/Dt) = 12.2*0.514;
v(1+23/Dt:24/Dt) = 7*0.514;

%% Electricity Price 

E_price = zeros(J,1);

%E_price = 0.13;        % Electricity price £ per kWh TIDE
E_price(1:7/Dt) = 0.0791;
E_price(8/Dt:16/Dt) = 0.1627;
E_price(17/Dt:20/Dt) = 0.3255;
E_price(21/Dt:24/Dt) = 0.1627;

%% Load

Load(1:1/Dt) = L * 0.329;  %Mendeley hourly energy consumption data for commerical type consumer
Load(1+1/Dt:2/Dt) = L * 0.336;
Load(1+2/Dt:3/Dt) = L * 0.326;
Load(1+3/Dt:4/Dt) = L * 0.327;
Load(1+4/Dt:5/Dt) = L * 0.311;
Load(1+5/Dt:6/Dt) = L * 0.386;
Load(1+6/Dt:7/Dt) = L * 0.573;
Load(1+7/Dt:8/Dt) = L * 0.675;
Load(1+8/Dt:9/Dt) = L * 0.732;
Load(1+9/Dt:10/Dt) = L * 0.781;
Load(1+10/Dt:11/Dt) = L * 0.807;
Load(1+11/Dt:12/Dt) = L * 0.810;
Load(1+12/Dt:13/Dt) = L * 0.821;
Load(1+13/Dt:14/Dt) = L * 0.870;
Load(1+14/Dt:15/Dt) = L * 0.831;
Load(1+15/Dt:16/Dt) = L * 0.822;
Load(1+16/Dt:17/Dt) = L * 0.792;
Load(1+17/Dt:18/Dt) = L * 0.799;
Load(1+18/Dt:19/Dt) = L * 0.802;
Load(1+19/Dt:20/Dt) = L * 0.783;
Load(1+20/Dt:21/Dt) = L * 0.727;
Load(1+21/Dt:22/Dt) = L * 0.608;
Load(1+22/Dt:23/Dt) = L * 0.384;
Load(1+23/Dt:24/Dt) = L * 0.323;

%% Solar Output

G_c(1:6/Dt) = 0;    % https://geoflow.com.au/wp-content/uploads/2018/04/word-image-28.png
G_c(1+7/Dt:8/Dt) = 20;
G_c(1+8/Dt:9/Dt) = 192;
G_c(1+9/Dt:10/Dt) = 450;
G_c(1+10/Dt:11/Dt) = 692;
G_c(1+11/Dt:12/Dt) = 850;
G_c(1+12/Dt:13/Dt) = 954;
G_c(1+13/Dt:14/Dt) = 948;
G_c(1+14/Dt:15/Dt) = 856;
G_c(1+15/Dt:16/Dt) = 684;
G_c(1+16/Dt:17/Dt) = 458;
G_c(1+17/Dt:18/Dt) = 200;
G_c(1+18/Dt:19/Dt) = 26;
G_c(1+19/Dt:20/Dt) = 5;
G_c(1+20/Dt:24/Dt) = 0;

%% Initial and Expected State of Charge

S_0m = 0.25;
S_0 = zeros(n,1);
S_Em = 0.75;
S_E = ones(n,1);

for e = 1:n
    S_0(e) = (S_0m + 0.1*rand(1)-0.1*rand(1));
    S_E(e) = (S_Em + 0.1*rand(1)-0.1*rand(1));    
end    

%% Time in and time out 
%T_in = ones(n,1)+1;
%T_out = J*ones(n,1)-1;

T_in = zeros(n,1);
T_out = zeros(n,1);

for e = 1:n
    T_in(e,1) = 7*J/24 + round((1.3/Dt)*randn(1),0);        % Normal distribution around 7am
    T_out(e,1) = 18*J/24 + round((1.3/Dt)*randn(1),0);      % Normal distribution around 6pm 
end

end