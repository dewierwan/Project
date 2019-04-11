function [c,ceq] = PSOcon(x,y,~,~,S_0,S_E,~,T_in,T_out,~,~,~)

Z = y(1); 
S_max = y(2); 
S_min = y(3);
n = y(4); 
J = y(5);
%alpha_s = y(10);
c = zeros(1,2*n*J);
ceq = zeros(1,n*J);
S = zeros(n+J*n,1);
%x = x_s / alpha_s;      %Unscaling

%% SOC, opertionality, charge/discharge power

for k = 1:J
    for e = 1:n
        if k == T_in(e)
            S(n*(k-1)+e) = S_0(e);
        end
                 
        if k >= T_in(e) && k < T_out(e)
            c(2*(n*(k-1)+e)-1) = (S(n*(k-1)+e) + x(n*(k-1)+e) * Z) - S_max;         % SoC <= SoC_Max
            c(2*(n*(k-1)+e)) = - (S(n*(k-1)+e) + x(n*(k-1)+e) * Z) + S_min;         % SoC >= SoC_Min
            S(n*k+e) = S(n*(k-1)+e) + x(n*(k-1)+e) * Z;                             % Update SoC
        end
        
        if k == T_out(e)
            c(2*(n*(k-1)+e)-1) = (S(n*(k-1)+e) + x(n*(k-1)+e) * Z) - S_E(e) - 0.05; % Final SoC ~= Expected SoC
            c(2*(n*(k-1)+e)) = - (S(n*(k-1)+e) + x(n*(k-1)+e) * Z) + S_E(e) - 0.05; 
        end
        
        if k < T_in(e) || k > T_out(e)
            ceq(n*(k-1)+e) = x(n*(k-1)+e);
        end
    end
end

%% Scaling
%c_s = c * alpha_s;
%ceq_s = ceq * alpha_s;
end