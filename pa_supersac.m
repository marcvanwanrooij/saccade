function SupSac = pa_supersac(Sac,Stim,type,stimnr,sfreq)

% SUPSAC = SUPERSAC(SAC,STIM,TYPE,STIMNR,SFREQ)
%
% Add stimulus info of a single stimulus during the saccadic trial
% to the sac-matrix for a more convenient analysis.
% Input are the SAC and STIM matrices. Optional input is the modality
% of the particular stimulus (default: 2 = SND1, see also INDEX), and 
% stimulus nr in the trial (default: 1 = first stimulus of the correct 
% type in that trial).
%
% SUPERSAC also recalculates Saccade Latency (column 5) to accomodate the
% stimulus type and number, if you supply the sample frequency.
%
%
% See also SUPINDEX, MATCHSAC, SPLIT
%
%  02-05-01 Marcus
% a minor change

%% Initialization
if nargin<3
    type        = 2; % SND
end
if nargin<4
    stimnr      = 1; % first stimulus in a trial of the correct type
end
if nargin<5
    sfreq       = 1000;
end
if nargin<2
    Sac     = [];
    Stim    = [];
    fname   =[];
    fname   = fcheckexist(fname,'.mat');
    load(fname);
    SupSac  = pa_supersac(Sac,Stim); %#ok<NASGU>
end


%% Convert from intensities
Stim            = pa_int2int(Stim);

%% Let's make a SuperSaccade matrix

STM             = pa_matchsac(Sac,Stim,type,stimnr);
Sac             = [Sac STM(:,2:11)];

%% Remove NaNs
sel = ~isnan(STM(:,1));
Sac = Sac(sel,:);

%% Correct Saccade Latency
if ~isempty(sfreq)
    Sac(:,5)    = 1000*Sac(:,3)/sfreq-Sac(:,27);
end

%% And call it SupSac
SupSac          = Sac;