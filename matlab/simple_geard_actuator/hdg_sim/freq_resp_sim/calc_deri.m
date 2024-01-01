 % calculate derivative of state value
 function ret = calc_deri(type, xvec, Jm, Jl, Bm, Bl, Ks, Ds, n, taum, taud, taul)
    if type == 0 %type1
        tmp(1, 1) = xvec(2, 1);
        tmp(2, 1) = -Ks / Jm * xvec(1, 1) -Bm / Jm * xvec(2, 1) - Ds / Jm * xvec(2, 1) + n * Ks / Jm * xvec(3, 1) + n * Ds / Jm * xvec(4, 1) + taum / Jm - taud / Jm; 
        tmp(3, 1) = xvec(4, 1);
        tmp(4, 1) = n * Ks / Jl * xvec(1, 1) + n * Ds / Jl * xvec(2, 1) - n^2 * Ks / Jl * xvec(3, 1) - Bl / Jl * xvec(4, 1) - n^2 * Ds / Jl * xvec(4, 1) - taul / Jl; 
    elseif type == 1 %type2
        tmp(1, 1) = xvec(2, 1);
        tmp(2, 1) = -Ks / (Jm * n^2) * xvec(1, 1) - Ds / (Jm * n^2) * xvec(2, 1) + Ks / (Jm * n) * xvec(3, 1) + Ds / (Jm * n)  * xvec(4, 1) - Bm / Jm * xvec(2, 1)  + taum / Jm - taud / Jm; 
        tmp(3, 1) = xvec(4, 1);
        tmp(4, 1) = Ks / (Jl * n) * xvec(1, 1) + Ds / (Jl * n) * xvec(2, 1) - Ks / Jl * xvec(3, 1) - Ds / Jl * xvec(4, 1) - Bl / Jl * xvec(4, 1) - taul / Jl; 
    end
    
    ret = tmp;