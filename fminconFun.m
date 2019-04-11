function [f,g] = fminconFun(x,y,P_gen,Load,S_0,~,E_price,T_in,T_out,w1,w2,w3)
%% Initialise 
Z = y(1); 
n = y(4); 
J = y(5);
Dt = y(6);
%w1 = y(7);     % Cost
%w2 = y(8);     % S + Zx - 0.5
%w3 = y(9);     % P^2
%alpha_s = y(10);
E_cost = zeros(J,1);
P_sqr = zeros(J*n,1);
L_res = zeros(J*n,1);
P_ev = zeros(J*n,1);
SumP_ev = zeros(J,1);
%x = x_s / alpha_s;              %Unscaling

%% SoC
S = zeros(n*J+n);
for k = 1:J
    for e = 1:n        
%% Initial SOC        
        if k == T_in(e)
            S(n*(k-1)+e) = S_0(e);
        end
        if (k >= T_in(e)) && (k <= T_out(e))
%% Reduce Power per EV
            P_sqr(n*(k-1)+e) = x(n*(k-1)+e)^2;
%% Move SoC towards point of least resistance
            L_res(n*(k-1)+e) = abs(S(n*(k-1)+ e) + Z * x(n*(k-1)+e) - 0.5);
            S(n*k+e) = S(n*(k-1)+e) + x(n*(k-1)+e) * Z;
%% Sum total EV power 
            P_ev(n*(k-1)+e) = x(n*(k-1)+e);
        end
    end
%% Reduce Total Electricity Cost
    SumP_ev(k) = sum(P_ev(n*(k-1)+1:n*k));
    E_cost(k) = 1e3*(E_price(k)/(1000/Dt))*(Load(k)+SumP_ev(k)-P_gen(k));        
end
%% Gradient and Hessian
g = zeros(n*J,1);
if nargout > 1 % gradient required
    for k = 1:J
        for e = 1:n
            if k >= T_in(e) && k <= T_out(e)
                R = zeros(J,1);
                for u = k:J
                    y = zeros(u,1);
                    for v = 1:u
                        y(v) = x(n*(v-1)+e);
                    end
                    B = S_0(e)+Z*sum(y)-0.5;
                    R(u) = abs(B)/B;
                end
                   g(n*(k-1)+e) = w1*1e3*E_price(k)/(1000/Dt)...
                    + w2*Z*sum(R)...
                    + w3*2*x(n*(k-1)+e);
            end
        end
    end 
end
%% Objective Function
f =  w1*sum(E_cost) + w2*sum(L_res) + w3*sum(P_sqr);    %Unscaled objective function
%f_s = f * alpha_s;                                  %Scaled variable
%g_s = g * alpha_s;
%disp('ObjFunEnd');
end