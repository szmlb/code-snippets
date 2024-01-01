% update (runge=kutta 4th itnegration)
function ret = rk4_update(type, xvec, time, dt, param)
        
        %{
        % parameter map
        Jm = param(1);
        Jl = param(2);
        Bm = param(3);
        Bl = param(4);
        Ks = param(5);
        Ds = param(6);
        n = param(7);
        taum = param(8);
        taud = param(9);
        taul = param(10);
        %}

        k1 = calc_deri(type, xvec, param(1), param(2), param(3), param(4), param(5), param(6), param(7), param(8), param(9), param(10));
        k2 = calc_deri(type, xvec + dt / 2.0 * k1, param(1), param(2), param(3), param(4), param(5), param(6), param(7), param(8), param(9), param(10));
        k3 = calc_deri(type, xvec + dt / 2.0 * k2, param(1), param(2), param(3), param(4), param(5), param(6), param(7), param(8), param(9), param(10));
        k4 = calc_deri(type, xvec + dt * k3, param(1), param(2), param(3), param(4), param(5), param(6), param(7), param(8), param(9), param(10));
        
        xvec = xvec + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);   
        ret = xvec;