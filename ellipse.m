% ellipse(Xo,Yo,L,S,Phi,Sty)
%  draw an ellipse with long and short axes L and S
%  with orientation Phi (in deg) at point Xo,Yo.
%
%  Marc van Wanrooij

function h = ellipse(Xo,Yo,L,S,Phi,Sty)
if nargin<6
	Sty = 'r';
end
DTR = pi/180;
Phi = Phi*DTR;
wt  = (0:360).*DTR;
X   = Xo + L*cos(Phi)*cos(wt) - S*sin(Phi)*sin(wt);
Y   = Yo + L*sin(Phi)*cos(wt) + S*cos(Phi)*sin(wt);
% h	= plot(X,Y,Sty,'LineWidth',2);


h = patch(X,Y,Sty);
hold on
alpha(h,.2);
set(h,'EdgeColor',Sty,'LineWidth',2);

% wt = [0 180]*DTR;
% X   = Xo + L*cos(Phi)*cos(wt) - S*sin(Phi)*sin(wt);
% Y   = Yo + L*sin(Phi)*cos(wt) + S*cos(Phi)*sin(wt);
% plot(X,Y,'-','Color',Sty);
% 
% wt = [90 270]*DTR;
% X   = Xo + L*cos(Phi)*cos(wt) - S*sin(Phi)*sin(wt);
% Y   = Yo + L*sin(Phi)*cos(wt) + S*cos(Phi)*sin(wt);
% plot(X,Y,'-','Color',Sty);
