T_values = [1, 0.98, 0.96, 0.94, 0.92, 0.9, 0.85, 0.8, 0.7, 0.6];
itmax = 100;
tol = 1e-10;

V1_results = zeros(length(T_values), 1);
V2_results = zeros(length(T_values), 1);
P_results = zeros(length(T_values), 1);

Oa = 1/(9*(2^(1/3)-1));
Ob = (2^(1/3)-1)/3;

for i = 1:length(T_values)
    T = T_values(i);
    assignin('base', 'T', T);
    x0 = [0.4; 5.0];
    sol = vnewton_silent(x0, itmax, tol, T, Oa, Ob);
    V1_results(i) = sol(1);
    V2_results(i) = sol(2);
    P_results(i) = 3*T/(sol(1)-3*Ob) - 9*Oa/(T^0.5*sol(1)*(sol(1)+3*Ob));
end

fprintf('%-6s  %-14s  %-14s  %-14s\n', 'T', 'V1', 'V2', 'P');
for i = 1:length(T_values)
    fprintf('%-6.2f  %-14.8f  %-14.8f  %-14.8f\n', T_values(i), V1_results(i), V2_results(i), P_results(i));
end

save('phase_data.mat', 'T_values', 'V1_results', 'V2_results', 'P_results');

function sol = vnewton_silent(x0, itmax, tol, T, Oa, Ob)
    x = x0;
    for k = 1:itmax
        fx = f_rk(x, T, Oa, Ob);
        Jx = jac_rk(x, T, Oa, Ob);
        x = x + Jx\(-fx);
        if norm(f_rk(x, T, Oa, Ob)) <= tol
            break
        end
    end
    sol = x(1:2);
end

function y = f_rk(x, T, Oa, Ob)
    V1 = x(1); V2 = x(2);
    P1  = 3*T/(V1-3*Ob) - 9*Oa/(T^0.5*V1*(V1+3*Ob));
    P2  = 3*T/(V2-3*Ob) - 9*Oa/(T^0.5*V2*(V2+3*Ob));
    mu1 = 3*V1*T/(V1-3*Ob) - 9*Oa/(T^0.5*(V1+3*Ob)) - 3*T*log(V1-3*Ob) - (3*Oa/(T^0.5*Ob))*log((V1+3*Ob)/V1);
    mu2 = 3*V2*T/(V2-3*Ob) - 9*Oa/(T^0.5*(V2+3*Ob)) - 3*T*log(V2-3*Ob) - (3*Oa/(T^0.5*Ob))*log((V2+3*Ob)/V2);
    y(1,1) = P1 - P2;
    y(2,1) = mu1 - mu2;
end

function J = jac_rk(x, T, Oa, Ob)
    V1 = x(1); V2 = x(2);
    dPdV1 = -3*T/(V1-3*Ob)^2 + 9*Oa*(2*V1+3*Ob)/(T^0.5*V1^2*(V1+3*Ob)^2);
    dPdV2 = -3*T/(V2-3*Ob)^2 + 9*Oa*(2*V2+3*Ob)/(T^0.5*V2^2*(V2+3*Ob)^2);
    J = [dPdV1, -dPdV2; V1*dPdV1, -V2*dPdV2];
end
