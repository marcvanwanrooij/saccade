function id = pa_instdir(R, Phi, A)

% FUNCTION id = instdir(R, Phi, A)
%
%  This function gives the value of Phi at R = 0.3*A or at the point
%  closest to that. R, Phi are vectors; A = scalar

dr = abs(R-0.3*A);
i  = find(dr==min(dr));
id = mean(Phi(i));
