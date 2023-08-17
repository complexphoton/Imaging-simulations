function win_prof = planckwin(N, p)
%% Generate an N-point Planck-taper window, where p controls the window width.
    win_prof = zeros(N, 1);
    pN = round(p*N);
    x1 = 2:(pN-1);
    trans = @(x) 1./(1+exp(pN./x-pN./(pN-x)));

    win_prof(x1) = trans(x1);

    x2 = pN:(N-pN-1);
    win_prof(x2) = 1;

    x3 = (N-pN):N;
    win_prof(x3) = trans(N-x3);
end
