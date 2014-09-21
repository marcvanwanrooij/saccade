function [expinfo,chaninfo,cLog] = pa_readevent_sac(logfile)
% Extract experimental parameters from a csv-file
%
% [EXPINFO,CHANINFO,LOG] = READCSV(CSVFILE)
%
% The EXPINFO-vector contains experimental information,
% extracted directly from the EXP-file, with column:
%               1 = 0
%               2 = Maximum Number of trials (as set in exp-file)
%               3 = Number of repeats
%               4 = Number of '=>' (trial indicator in exp-file)
%               5 = ITI start (ITI: Intertrial Interval)
%               6 = ITI stop
%               7 = Random type (0 = no, 1 = per set, 2 = all trials)
%               8 = Number of channels
%
% The CHANINFO-matrix contains data acquistion channel information,
% extracted directly from the CFG-file, with each row corresponding to a
% single channel and each column corresponding to the following parameter:
%               1 = 0
%               2 = Channel number
%               3 = Channel name
%               4 = Low-Pass frequency
%               5 = Sample rate
%               6 = Number of samples
%
% For example, to obrain the sample rate of channel 2:
%  >> chaninfo(2,5) 
%
%
% The LOG-matrix contains trial information
%
%   1) =  Trial number (1 = trial no. 1)
%   2) =  Stimulus number (e.g. Acq-Led-Led-Snd -> 1 2 3 4, but trialnr stays the same)
%   3) =  random ITI calculated
%   4) =  random ITI actual (is larger when for instance disk is too busy)
%   5) =  Modality of stimulus (0 = LED, 1 = SKY, 2 = SND1, 3 = SND2, 4 = Acquisition, 5 = Trg0, 6 = Input 1, 7 = Input 2)
%   6) =  Stimulus location x BOOG (degrees, -180..180)
%   7) =  Stimulus location y speaker (number)
%   8) =  Stimulus onset
%   9) =  Stimulus offset
%   10) = Stimulus intensity (for LED: 0 = lowest ... 7 = highest, for SND: 0 = lowest ... 100 = higest)
%   11) = Stimulus attribute (For LED: 0 = red, 1 = green; for SND: XXX (fragment of wav-name -> sndXXX.wav))
%   12) = delay of second sound
%   13) = line number in EXP-file
%
%                       
%   Type of trial effects which columns are not zero. The columns used
%   are presented below:
%
%                      ITI ITI
%              Trl Stm cal act Typ  x   y  On  Off Int Atr Bit Edg
%               1   2   3   4   5   6   7   8   9   10  11  12  13
%       Led:    X   X   X   X   X   X   X   X   X   X   X       X   
%       Snd:    X   X   X   X   X   X   X   X       X   X       X
%       Acq:    X   X   X   X   X           X                   X
%       Trg0:   X   X   X   X   X           X   X           X   X
% 
%   'Cells' that are empty in the logfile are denoted with 'NaN' in the
%   cLog-matrix
%
%  See also LOADDAT

% Copyright 2006
% Author: Tomg Oct 2006
% Modified by: Marcw 2007

%% Initialization
if nargin<1
    logfile             = '';
end
chaninfopar             = 6;
logpar                  = 13;

%% check INPUT %%%
logfile                 = fcheckext(logfile,'.csv');
logfile                 = fcheckexist(logfile);
if isempty(logfile)
    disp('No Log File has been chosen (readcsv)');
    expinfo             = [];
    chaninfo            = [];
    cLog                = [];
    return
end

%% Open csv-file
fid                     = fopen(logfile);
firstLine               = fgetl(fid);

%% Get Experimental Info from first line
expinfo(1,:)            = sscanf(firstLine,'%f;')';
Nchan                   = expinfo(1,8);
Ntrials                 = expinfo(1,4);


%% Get Channel Information from next Nchannels lines
chaninfo                = NaN*ones(Nchan,chaninfopar);
ident                   = 0;
k                       = 0;
while ~ident
    A                   = fscanf(fid,'%g;',chaninfopar)';
    ident               = A(1);
    if ~ident
        k               = k+1;
        chaninfo(k,:)   = A;
    end
end
chaninfo                = chaninfo(1:k,:);
if Nchan~=k
%     disp('Uh-oh - Number of relevant channels does not correspond to number of channel configs');
end

%% load cLog-file body
frewind(fid);
count                   = 0;
cLog                    = NaN*ones(Ntrials,logpar);
while ~feof(fid)
    curLine             = fgetl(fid);
    firstCell           = sscanf(curLine,'%g;',1);
    seps                = [0 findstr(';',curLine) length(curLine)+1];
    if firstCell
        count           = count+1;
        for i           = 2:length(seps)
            icol        = i-1;
            curcell     = curLine(seps(i-1)+1:seps(i)-1);
            if sum(isletter(curcell))
                if strcmpi(strtrim(curcell),'led')
                    cLog(count,icol)    = 0;
                elseif strcmpi(strtrim(curcell),'sky')
                    cLog(count,icol)    = 1;
                elseif strcmpi(strtrim(curcell),'snd') || strcmpi(strtrim(curcell),'snd1') 
                    cLog(count,icol)    = 2;
                elseif strcmpi(strtrim(curcell),'snd2')
                    cLog(count,icol)    = 3;
                elseif strcmpi(strtrim(curcell),'acq')
                    cLog(count,icol)    = 4;
                elseif strcmpi(strtrim(curcell),'trg0')
                    cLog(count,icol)    = 5;
                elseif strcmpi(strtrim(curcell),'inp1')
                    cLog(count,icol)    = 6;
                elseif strcmpi(strtrim(curcell),'inp2')
                    cLog(count,icol)    = 7;
                elseif strcmpi(strtrim(curcell),'nan')
                    cLog(count,icol)    = NaN;
                end
            elseif isempty(str2double(curcell))
                %Skip empty cell at end of file
                warning('READCSV:emptycell','from readcsv')
                disp(['   Cell (' num2str(count) ',' num2str(icol) ') is empty ... skip'])
                disp(['   logfile= ''' logfile ''')'])
            else
                cLog(count,icol)        = str2double(curcell);
            end
        end
    end
end


%% close Log-file %%
fclose(fid);