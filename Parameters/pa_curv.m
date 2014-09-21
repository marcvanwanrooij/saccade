%   calculate the Arend Smit curvature index for saccades.
% c=curv(h,v)
%   calculate the Arend Smit curvature index for saccades.
%  
%   Jeroen Goossens  

function c = pa_curv(h,v);

h = h(isnan(h)==0) - h(1);    %  saccade starts in (0,0)
v = v(isnan(v)==0) - v(1);

[p,r] = cart2pol(h,v);    % rotate to positive horizontal axis
p = p - p(length(p));
[h,v] = pol2cart(p,r);


c  = max(abs(v))/max(h);

if max(abs(v)) ~= max(v), 
  c = -1 * c; 
end;
