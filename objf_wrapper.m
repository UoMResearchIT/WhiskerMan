function [mode, objf, objgrd, user] = objf_wrapper(mode, n, x, nstate, user)
    % objf_wrapper   Wrapper function used to pass objective function
    %               in to nag e04dg routine
    % mode, n, nstate, and user are set in the NAG function e04dg
    % x is the starting guess at the minimum x-value of the objective
    % function
    z = x;
    r = user.r;
    np0 = user.np0;
    np1 = user.np1;
    imh = user.imh;
    imv = user.imv;
    dt = user.dt;
    sigma1 = user.sigma1;
    sigma2 = user.sigma2;
    calib = user.calib;
    
    [objf, objgrd] = Bimfun_wm4(z, r, np0, np1, imh, imv, dt, sigma1, sigma2, calib);
end