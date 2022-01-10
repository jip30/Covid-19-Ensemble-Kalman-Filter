function [t,st_up] = SIRDsolve(stk, tspan, param)

fun = @(t,st)odeSIRD(t,st,param); 

[t, st_up] = ode15s(fun, tspan, stk);


    function dstdt =  odeSIRD(t ,st, param)
        %State Variables
        S = st(1); I= st(2); R = st(3); D = st(4);
        
        %Parameters
        beta = param(1); gamma = param(2); delta = param(3); N= param(4);
        
        %ODE system
        dSdt =  -((beta)/N)*I*S; %add back in -
        dIdt = ((beta)/N)*I*S - gamma*I - delta*I;
        dRdt = gamma*I;
        dDdt = delta*I;
        
        dstdt = [dSdt; dIdt; dRdt; dDdt];
        
        
    end

end