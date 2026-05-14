function [rx,ry] = vnewton(x0,itmax,tol)
r=zeros(itmax+1,1);
x=x0;r(1,1)=norm(f(x0));
fprintf('\n');
disp(' iter  L2 norm of f(x)')
fprintf('%4.0f',0)
fprintf('%19.14f',r(1,1))
fprintf('\n');
for k=1:itmax
    J=jacobian(x);rhs=-f(x);
    delta=J\rhs;
    x=x+delta;
    r(k+1,1)=norm(f(x));
    fprintf('%4.0f',k)
    fprintf('%19.14f',r(k+1,1))
    fprintf('\n');
    if abs(r(k+1,1))<=tol
        fprintf('\n');
        disp(['Method has converged after ',num2str(k),' iterations.'])
        break
    end
end
if abs(r(k+1,1))>tol
    fprintf('\n');
    disp(['Method has NOT converged after ',num2str(k),' iterations.'])
end
fprintf('\n');
disp('Solution vector:')
x = x(1:2);
x
rx=zeros(k,1);ry=zeros(k,1);
for i=1:k
    rx(i,1)=r(i,1);
    ry(i,1)=r(i+1,1);
end

function y=f(x)
T  = evalin('base','T');
Oa = 1/(9*(2^(1/3)-1));
Ob = (2^(1/3)-1)/3;
V1 = x(1);
V2 = x(2);
P1  = 3*T/(V1-3*Ob) - 9*Oa/(T^0.5*V1*(V1+3*Ob));
P2  = 3*T/(V2-3*Ob) - 9*Oa/(T^0.5*V2*(V2+3*Ob));
mu1 = 3*V1*T/(V1-3*Ob) - 9*Oa/(T^0.5*(V1+3*Ob)) ...
    - 3*T*log(V1-3*Ob) - (3*Oa/(T^0.5*Ob))*log((V1+3*Ob)/V1);
mu2 = 3*V2*T/(V2-3*Ob) - 9*Oa/(T^0.5*(V2+3*Ob)) ...
    - 3*T*log(V2-3*Ob) - (3*Oa/(T^0.5*Ob))*log((V2+3*Ob)/V2);
y(1,1) = P1 - P2;
y(2,1) = mu1 - mu2;

function J=jacobian(x)
T  = evalin('base','T');
Oa = 1/(9*(2^(1/3)-1));
Ob = (2^(1/3)-1)/3;
J=zeros(2,2);
V1 = x(1);
V2 = x(2);
dPdV1 = -3*T/(V1-3*Ob)^2 + 9*Oa*(2*V1+3*Ob)/(T^0.5*V1^2*(V1+3*Ob)^2);
dPdV2 = -3*T/(V2-3*Ob)^2 + 9*Oa*(2*V2+3*Ob)/(T^0.5*V2^2*(V2+3*Ob)^2);
J(1,1) =      dPdV1;
J(1,2) =     -dPdV2;
J(2,1) = V1 * dPdV1;
J(2,2) = -V2 * dPdV2;
