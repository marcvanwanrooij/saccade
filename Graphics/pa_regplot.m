function [h,b,r]=pa_regplot(X,Y)

%            H = PA_REGPLOT(X,Y)
%
% PA_REGPLOT plots X vs Y, and performs linear regression on X and Y.
%
% See also PA_LOC, PA_PLOTLOC

% (c) 2011 Marc van Wanrooij
% E-mail: marcvanwanrooij@gmail.com

%% Initialization
X       = X(:)';
Y       = Y(:)';
mrkr    = 'o';

%% Regression

b = regstats(Y,X,'linear','beta');
b = b.beta;
gain                        = b(2);
bias                        = b(1);
r                           = corrcoef(X,Y);
r                           = r(2);
% corrpear                    = pa_pearson(X,Y);

%% Text
if bias>0
    linstr                  = ['Y = ' num2str(gain,2) 'X + ' num2str(bias,2) ];
elseif bias<=0
    linstr                  = ['Y = ' num2str(gain,2) 'X - ' num2str(abs(bias),2) ];
end
corrstr                     = ['r^2 = ' num2str(r^2,2)];

%% Graphics
h                           = plot(X, Y, ['k' mrkr]); set(h,'MarkerSize', 5, 'LineWidth',2,'MarkerFaceColor','w');
hold on
lsline;
axis square;
box off;
title(linstr)
text(range(1)+10,range(4)-10,corrstr,'HorizontalAlignment','left')
ax = axis;
% mm = minmax(ax);
% axis([mm mm]);
plot(ax([1 2]),gain*ax([1 2])+bias,'k-','LineWidth',2);
xlabel('X');
ylabel('Y');
