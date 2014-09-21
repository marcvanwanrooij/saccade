function pa_logcal2csvcal(fname)
%          PA_LOGCAL2CSVCAL(FNIN,FNOUT)
%                or
%          PA_LOGCAL2CSVCAL(FNIN)
%               or
%          PA_LOGCAL2CSVCAL
%
% Converts LOG-files (FNIN, optional) to CSV-files (FNOUT, optional). Data
% and file formats are used in the Biophysics department of the Donders
% Institute. This is only used for backwards-compatibility.
%
% See also PA_TDTCAL2DATCAL, PA_TDT2DAT, PA_LOG2CSV

% 2011 Marc van Wanrooij
% e-mail:marcvanwanrooij@neural-code.com

%% Initialization
if nargin==0
    fname = pa_fcheckexist([],'*.log');
elseif nargin==1
    fname = pa_fcheckext(fname,'log');
end
fout =  pa_fcheckext(fname,'csv');

%% Read
[stim,log] = pa_readstim(fname);
nchan       = log.NChan;
maxtrials   = log.NTest;
nrep        = 1;
ntrials     = log.NTest;
% 1) 0
% 2) Maximum number of trials per repetition
% 3) Number of repetitions
% 4) Number of trials (2)*3))
% 5) Inter Trial Interval start
% 6) Inter Trial Interval end
% 7) Randomization Type
% 8) Number of (Acq) channels
fid                     = fopen(fout,'w');

% firstLine = sprintf('0;%d;%d;%d;0;0;1;%d;',maxtrials,nrep,ntrials,nchan)
fprintf(fid,'0;%d;%d;%d;0;0;1;%d;\n',maxtrials,nrep,ntrials,nchan);

% 1) 0
% 2) Channel Number
% 3) Channel NumberFart Manual for Audio Dummies
% 4) Low-pass cut-off frequency
% 5) Desired Sampling rate
% 6) Number of samples
% chan1 = sprintf('0;1;1;150;%d;%d;',log.FSample,log.NSample)
% chan2 = sprintf('0;2;2;150;%d;%d;',log.FSample,log.NSample)
% chan3 = sprintf('0;3;3;150;%d;%d;',log.FSample,log.NSample)
fprintf(fid,'0;1;1;150;%d;%d;\n',log.FSample,100);
fprintf(fid,'0;2;2;150;%d;%d;\n',log.FSample,100);
fprintf(fid,'0;3;3;150;%d;%d;\n',log.FSample,100);

% Should become:
% 1) Trial number
% 2) Stimulus number in trial
% 3) Desired Inter Trial Interval
% 4) Actual Inter Trial Interval (can deviate from 3) when e.g. disk is too busy)
% 5) Modality of stimulus
% 6) Hoop location (degrees, -180..180) / Spoke (1-12)
% 7) Speaker/Led location (number)
% 8) Stimulus onset relative to trial onset
% 9) Stimulus offset relative to trial onset
% 10) Stimulus intensity (for LED: 0 = lowest ... 7 = highest, for SND: 0 = lowest ... 100 = higest)
% 11) Stimulus attribute (For LED: 0 = red, 1 = green; for SND: XXX (fragment of wav-name ->
% sndXXX.wav))
% 12) Bit (Micro-controller Trg0: 5; For SND: XXX (fragment of wav-name-> sndXXX.wav))
% 13) line number in EXP-file (which deviates from 1) when trials are randomized according to
% the EXP-file)
% But is
% Stimuli :
%    1 Trial Number
%    2 Number of Targets
%  Fixation Spot
%    3 Modality (1=Vis,2=Aud,3=Bi)
%    4 Amplitude                [deg]
%    5 Direction                [deg]
%    6 Azimuth                  [deg]
%    7 Elevation                [deg]
%    8 Onset  Time              [ms]
%    9 Offset Time              [ms]
%   10 Intensity                [ms]
%   11 Attribute (color/freq)   [ms]
%  Target 1
%   12 Modality (1=Vis,2=Aud,3=Bi)
%   13 Amplitude                [deg]
%   14 Direction                [deg]
%   15 Azimuth                  [deg]
%   16 Elevation                [deg]
%   17 Onset  Time              [ms]
%   18 Offset Time              [ms]
%   19 Intensity                [ms]
%   20 Attribute (color/freq)   [ms]
%  Target 2
%   Etc.

for ii  = 1:ntrials
    ntargets    = stim(ii,2);
    for jj  = 1:ntargets
        Mod     = stim(ii,3+(jj-1)*9);
        %         sprintf('%d;\t%d;\t0;\t0;\t%d;\t%f;\t%f;\t%d;\t%d;\t%d;\t%d;\tNaN;\t1;',stim(ii,1),jj,getmod(Mod),...
        %             stim(ii,(6)+(jj-1)*9),...
        %             stim(ii,(7)+(jj-1)*9),...
        %             stim(ii,(8)+(jj-1)*9),...
        %             stim(ii,(9)+(jj-1)*9),...
        %             stim(ii,(10)+(jj-1)*9),...
        %             stim(ii,(11)+(jj-1)*9))
        if Mod>0
            Modstr = getmod(Mod);
            
            fprintf(fid,'%d;\t%d;\t0;\t0;\t%s;\t%f;\t%f;\t%d;\t%d;\t%d;\t%d;\t1000;\t1\r',stim(ii,1),jj,Modstr,...
                stim(ii,(6)+(jj-1)*9),...
                stim(ii,(7)+(jj-1)*9),...
                stim(ii,(8)+(jj-1)*9),...
                stim(ii,(9)+(jj-1)*9),...
                stim(ii,(10)+(jj-1)*9),...
                stim(ii,(11)+(jj-1)*9));
        end
    end
end
%% write
% fprintf(fid, '%s', a_str);
fclose(fid);

function Modstr = getmod(Mod)
switch Mod
    case 0
        Modstr = 'Fix';
        %     Modstr = 0;
    case 1
        Modstr = 'Led';
    case 2
        Modstr = 'Snd';
end
