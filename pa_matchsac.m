function STM = pa_matchsac (Sac,Stim,type,stimnr)
% STIM = MATCHSAC(SAC,STIM,TYPE,STMNR)
%
% Insert or deletes lines in the STIM matrix such that it is matched, with 
% respect to the number of rows, to the SAC matrix.
% Inputs are the SAC and STIM matrices. Optional input is the modality
% of the particular stimulus (default: 1 = SND, see also INDEX), and 
% stimulus nr in the trial (default: 1 = first stimulus of the correct 
% type in that trial).
%
% See also INDEX, SUPERSAC

% 2012 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

if nargin<3
    type = 2;
end
if nargin<4
    stimnr = 1;
end

STM                     = NaN*ones(size(Sac,1),11);
for i                   = 1:size(Sac,1)
    trialNr             = Sac(i,1);
    indx                = find(Stim(:,1)  == trialNr & Stim(:,3) == type);
    if ~isempty(indx)
        indx            = indx(stimnr);
        STM(i,:)        = Stim(indx,1:11);
    end
end
