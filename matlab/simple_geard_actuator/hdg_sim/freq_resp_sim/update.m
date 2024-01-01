% update (euler itnegration)
function ret = update(xvec, dxvec, sampling_time)
        tmp(1, 1) = xvec(1, 1) + dxvec(1, 1) * sampling_time;
        tmp(2, 1) = xvec(2, 1) + dxvec(2, 1) * sampling_time;
        tmp(3, 1) = xvec(3, 1) + dxvec(3, 1) * sampling_time;
        tmp(4, 1) = xvec(4, 1) + dxvec(4, 1) * sampling_time;
        ret = tmp;
