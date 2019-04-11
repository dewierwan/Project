%profile on
clc;clear;close all;

%% Problems 

%% System definitions

[y,S_0,S_E,E_price,G_c,Load,T_in,T_out,v] = calcs();

Z = y(1); 
n = y(4); 
J = y(5);
Dt = y(6);
%w1 = y(7);     % Cost
%w2 = y(8);     % S + Zx - 0.5
%w3 = y(9);     % P^2
Pmax = y(10);
%alpha_s = y(10);

P_st = 30;           %Rated PV power
G_st = 1000;            %Solar Irradiance
nu = 0.002;             %Temp coefficient
T_st = 278.15;          %Standard module temp

P_r = 30;             %Rated WT power - Evance Iskra R-9000 5kW 5.5m 
v_in = 3;               %Cut in speed
v_r = 12;                %Rated speed
v_out = 30;              %Cut out speed

T_c = zeros(J,1);       %Module temp
P_pv = zeros(J,1);      %Solar Power
P_w = zeros(J,1);       %Wind Power
P_gen = zeros(J,1);
%v = zeros(J,1);         %Wind speed
v=v.';

for k = 1:J
      
    T_c(k,1)=270;
    
%PV Output. G_c(1,k) and T_c(1,k) are the variables.     
    
    P_pv(k,1) = P_st * (G_c(k,1)/G_st) * (1 + nu * (T_c(k,1)-T_st));
    
%Wind Turbine Output

    if (v_in <= v(k,1)) &&  (v(k,1)< v_r)
        P_w(k,1) = (P_r * (v_in .^ 3)) / ((v_in .^ 3) - (v_r .^ 3)) + ...
            (P_r * (v(k,1) .^ 3)) / ((v_r .^ 3) - (v_in .^ 3));
    elseif (v_r <= v(k,1)) && (v(k,1) < v_out)
        P_w(k,1) = P_r;
    else
        P_w(k,1) = 0;
    end
    
%Create random Load demand 

    P_w(k) = round(P_w(k),0);
    P_pv(k) = round(P_pv(k),0);
    Load(k,1) = round(Load(k,1),0);
    P_gen(k) = P_pv(k) + P_w(k);
   
end

%% Optimization
disp('Optimization in Progress');
fprintf('Total number of EVs: %0.0f\n',n);
fprintf('Total number of time steps: %0.0f\n',J);

tic

%% FMINCON %%%%
x0 = 1.5*ones(n*J,1);

for k = 1:J
    for e = 1:n
        if k < T_in(e) || k > T_out(e)
            x0(n*(k-1)+e) = 0;
        end
    end
end

%x0_s = alpha_s*x0;
lb = - Pmax * ones(n*J,1);
ub = Pmax * ones(n*J,1);
options = optimoptions(@fmincon,'SpecifyConstraintGradient',true,'SpecifyObjectiveGradient',true);
xopt = zeros(110,n*J);
j = 1;
for w1 = 0:0.1:1
    for w2 = 1 - w1
            w3 = 1e-2;
            Fun = @(x) fminconFun(x,y,P_gen,Load,S_0,S_E,E_price,T_in,T_out,w1,w2,w3);
            NLcon = @(x) fminconCon(x,y,P_gen,Load,S_0,S_E,E_price,T_in,T_out,w1,w2,w3);
            [xopt(j,:), fopt] = fmincon(Fun,x0,[],[],[],[],lb,ub,NLcon,options); 
            disp(j);
            j = j+1;
    end
end

%% Fmincon w/o gradient

%{
options = optimoptions(@fmincon);
Fun = @(x) fminconFun(x,y,P_gen,Load,S_0,S_E,E_price,T_in,T_out,w1,w2,w3);
            NLcon = @(x) fminconCon(x,y,P_gen,Load,S_0,S_E,E_price,T_in,T_out,w1,w2,w3);
            tic
            [xoptFw, foptFw] = fmincon(Fun,x0,[],[],[],[],lb,ub,NLcon,options); 
            toc

%% GA   

Pfun = @(x) PSOobj(x,y,P_gen,Load,S_0,S_E,E_price,T_in,T_out,w1,w2,w3);
PNLcon = @(x) PSOcon(x,y,P_gen,Load,S_0,S_E,E_price,T_in,T_out,w1,w2,w3);

options = optimoptions(@ga);
tic
[xoptGA, foptGA] = ga(Pfun,n*J,[],[],[],[],lb,ub,PNLcon,options);
toc
xoptGA = xoptGA.';
   
%% PSO %%%%
tic
disp('waiting for pso');
[xoptP, foptP] = pso(Pfun,n*J,[],[],[],[],lb,ub,PNLcon);
disp('still waiting?');
toc
xoptP = xoptP.';

%}

%% Tables, graphs, SOC, etc. 



%% State of Charge

I = 1;      % Iteration number

S = zeros(n*J,1);

for e = 1:n
    S(e) = S_0(e);
end
    
for k = 1:J
    for e = 1:n
        S(n*k+e) = S(n*(k-1)+e) + xopt(I,n*(k-1)+e)*Z; % S_n+1 = S_n + P_n
    end
end

SOC  = reshape(S,[n,J+1]);

%% Table

totalP_ev = zeros(J,1);
E_cost = zeros(J,1);

for k = 1:J
    totalP_ev(k) = sum(xopt(I,n*(k-1)+1:n*k));
    E_cost(k) = 1e3*(E_price(k)/(1000/Dt))*(Load(k)+totalP_ev(k)-P_gen(k));
end

Time = transpose(1:J);
Generated_Power_kW = P_gen;
Total_EV_Power_kW = totalP_ev;
Electricity_Cost = round(E_cost,3);
Load_kW = Load;

T = table(Time,Generated_Power_kW,Load_kW,Total_EV_Power_kW,Electricity_Cost);
disp(T);
fprintf('State of Charge of each EV per time step:\n\n');
disp(SOC.');

TotalCost = sum(E_cost);
if TotalCost > 0
    fprintf('Total Electricity Cost = £%0.2f\n',TotalCost);
else
    TotalCost = abs(TotalCost);
    fprintf('Total Electricity Income = £%0.2f\n',TotalCost);
end

%% Graphs 

figure
yyaxis left
plot(1:J,P_gen,'r--o',1:J,Load,'b-',1:J,totalP_ev,'k:')
xlabel('Hour');
ylabel('Power [kW]');
yyaxis right
plot(1:J,E_price)
ylabel('TIDE Electricity Price [£/kWh]');
ylim([0 0.40]);
legend({'Generated Power','Load Power','Total EV Power','TIDE Price'},'Location','best','NumColumns',2);