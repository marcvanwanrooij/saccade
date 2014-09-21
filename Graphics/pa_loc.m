function [corr,gain,gainerr,bias]=pa_loc(Xpoints,Ypoints,which,range,inc)

% [CORR,CORERR,GAIN,GAINERR,FITERR,LOCERR,BIAS,INSIGNERR] =
%            LOC(XPOINTS,YPOINTS,WHICH,RANGE,INC)
%
% LOC has a two-fold function:
%
% a) Plot X vs Y, usually Stimulus-azimuth vs Response-azimuth
% b) Performs linear regression on X and Y:
%        Y=aX+b
%    and returns Correlation, Gain (a), bias (b), their respective errors,
%    and the number of data points.
%
% RANGE states the range in which the Data-points can vary; a data-point
% with a larger range is discarded for further analysis in this m-function.
% This is necessary, especially when responses fall outside the unambiguous
% range created by the magnetic fields.
% RANGE is optional, as is INC, which is used for determining XTicks.
% Default value is [-60 60 -60 60].
% WHICH consists of 'az' or 'el', and is also optional.
%
%

% (c) 2011 Marc van Wanrooij
% E-mail: marcvanwanrooij@gmail.com

%% Initialization
if (nargin<3)
    which   = 'az';
end
if (nargin<4)
    range   = [-90 90 -90 90];
end
if (nargin<5)
    inc     = 30;
end
numpoints           = sum(Ypoints<range(1) | Ypoints>range(2));
if numpoints
    disp(['Warning: ' int2str(numpoints) ' points fall outside the specified range.']);
    disp('These points are still evaluated in the regression.');
end

%% Regression
Xpoint                      = [Xpoints ones(size(Xpoints))];
b                           = regress(Ypoints,Xpoint);
gain                        = b(1);
bias                        = b(2);
b2                          = bootstrp(100,@regress,Ypoints,Xpoint);
gainerr                     = std(b2(1,:));
% biaserr                     = std(b2(2,:));
corr                        = corrcoef(Xpoints,Ypoints);
corr                        = corr(2);
corrpear                    = pa_pearson(Xpoints,Ypoints);

%% Text
switch which
    case 'az'
        mrkr = '.';
        if bias>0
            linstr              = ['\alpha_R = ' num2str(gain,2) '\alpha_T + ' num2str(bias,2) ];
        elseif bias<=0
            linstr              = ['\alpha_R = ' num2str(gain,2) '\alpha_T - ' num2str(abs(bias),2) ];
        end
    case 'el'
        mrkr = 'o';
        if bias>0
            linstr              = ['\epsilon_R = ' num2str(gain,2) '\epsilon_T + ' num2str(bias,2)];
        elseif bias<0
            linstr              = ['\epsilon_R = ' num2str(gain,2) '\epsilon_T - ' num2str(abs(bias),2) ];
        end
end
corrstr                         = ['r^2 = ' num2str(corrpear^2,2)];

%% Graphics
h                           = plot(range([1 2]), range([3 4]), 'k-'); set(h,'LineWidth',1,'Color',[0.7 0.7 0.7]);
hold on;
h                           = plot(range([1 2]), [0 0], 'k-'); set(h,'LineWidth',1,'Color',[0.7 0.7 0.7]);
h                           = plot([0 0], range([3 4]), 'k-'); set(h,'LineWidth',1,'Color',[0.7 0.7 0.7]);
h                           = plot(range([1 2]), gain*range([1 2])+bias,'k--'); set(h,'LineWidth',2,'Color',[0.5 0.5 0.5]);
h                           = plot(Xpoints, Ypoints, ['k' mrkr]); set(h,'MarkerSize', 5, 'LineWidth',2,'MarkerFaceColor','w');
axis(range); 
axis square;
box on;
% grid on;
title(linstr)
% text(0,range(3)+10,linstr,'HorizontalAlignment','center')
text(range(1)+10,range(4)-10,corrstr,'HorizontalAlignment','left')

%--------------- Graphic settings ----------------------%
set(gca, 'YTick', range(1):inc:range(2));
set(gca, 'XTick', range(3):inc:range(4));
