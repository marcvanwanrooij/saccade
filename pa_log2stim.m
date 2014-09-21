function Stim = pa_log2stim(Log)
% Convert log-matrix to Stim-matrix
%
% Stim = LOG2Stim(LOG)
%
%   converts Log matrix to Stim matrix
%
%       Stim(:,1) = Trial number
%       Stim(:,2) = Stimulus number
%       Stim(:,3) = Type of Stimulus
%                       0 = LED
%                       1 = SKY
%                       2 = SND1
%                       3 = SND2
%                       4 = Acquisition (of data)
%                       5 = Trigger
%       Stim(:,4) = Stimulus location AZIMUTH
%       Stim(:,5) = Stimulus location ELEVATION
%       Stim(:,6) = Stimulus location R
%       Stim(:,7) = Stimulus location PHI
%       Stim(:,8) = Stimulus ONSET
%       Stim(:,9) = Stimulus OFFSET
%       Stim(:,10)= Stimulus INTENSITY
%                      for LED:
%                      0 = lowest ... 7 = highest
%                      for SND:
%                      0 = lowest ... 100 = higest
%       Stim(:,11)= Stimulus ATTRIBUTE
%                       For LED:
%                       ? = red
%                       ? = green
%                       For SND:
%                       frequency??
%
%       N.B. when value = NaN then this is not used
%
%   See also: INDEX, PA_READCSV, PA_READDAT, PA_AZELRPHI, PA_FART2AZEL

% (c) 21 Feb 2008
% TomG Oct 2006
% MarcW March 2007
% e-mail: marcvanwanrooij@neural-code.com

%% Get speaker location
X                   = Log(:,6);
Y                   = Log(:,7);
% Convert to Double-Polar
[AZ,EL]             = pa_fart2azel(X,Y);

sel		= Y==31;
AZ(sel) = 30;
EL(sel) = 45;
sel		= Y==30;
AZ(sel) = 30;
EL(sel) = -30;

% Exception: Sky-LED
sel                 = Log(:,5)     == 1; % SkyLED
[AZ(sel),EL(sel)]   = pa_sky2azel(Y(sel),X(sel));

% Exception: TDT2.0 set-up
if Log(1,12)==1000; % Bit is set to 1000 for TDT2 set-up in PA_LOGCAL2CSVCAL
    AZ = X;
    EL = Y;
end

% Convert to Single-Polar
[R,PHI]             = pa_azel2pol(AZ,EL);

%% Create actual Stimulus-matrix
Stim                = [Log(:,1), ...
    Log(:,2), ...
    Log(:,5), ...
    AZ, ...
    EL, ...
    R, ...
    PHI, ...
    Log(:,8), ...
    Log(:,9), ...
    Log(:,10), ...
    Log(:,11)];
sel = ~isnan(Log(:,12));
Stim(sel,8) = Stim(sel,8)+Log(sel,12);
Stim(sel,9) = Stim(sel,9)+Log(sel,12);

%     Log(:,8)+Log(:,12), ...
%     Log(:,9)+Log(:,12), ...
