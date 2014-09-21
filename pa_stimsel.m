function [Stim,indx] = pa_stimsel(Stim,type,stimnr)

% SUPSAC = SUPERSAC(STIM,TYPE,STIMNR,SFREQ)
%
% Select stimulus info of a single stimulus during the trial
%
% Input is the STIM matrix. Optional input is the modality
% of the particular stimulus (default: 2 = SND1, see also INDEX), and
% stimulus nr in the trial (default: 1 = first stimulus of the correct
% type in that trial).
%
% See also PA_SUPERSAC, PA_SUPINDEX, PA_MATCHSAC, SPLIT

% (c) 2012 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

%% Initialization
if nargin<3
	type        = 2; % SND
end
if nargin<4
	stimnr      = 1; % first stimulus in a trial of the correct type
end
if nargin<2
	fname   = [];
	fname   = fcheckexist(fname,'.mat');
	S		= load(fname);
	Stim	= S.Stim;
	Stim	= pa_stimsel(Stim);
end


%% Convert from intensities
Stim            = pa_int2int(Stim);

%% Let's make a selected Stim matrix
sel  = Stim(:,3) == type;
Stim = Stim(sel,:);

uS = unique(Stim(:,2));
nS = numel(uS);

if nS>1
	dim		= cat(1,true,diff(Stim(:,2))<0);
	indx	= find(dim)+stimnr-1;
	Stim = Stim(indx,:);
else
	indx = [];
end
