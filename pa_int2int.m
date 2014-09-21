function intStim = int2int(attStim)
% STIM = INT2INT(STIM)
%
% Convert TDT intensities (0-100%) to sound level (dBA)
%
% Author: Marcus

intStim             = attStim;
% As measured by Marc for BB, HP and LP wav-stimuli
% 48% = 60 dB
% 38% = 50 dB
% 28% = 40 dB
sel                   = intStim(:,3) == 1;
intStim(sel,10)       = attStim(sel,10)+12;
