function [Stim, Log] = pa_readstim (fname)

% [STIM, LOG] = PA_READSTIM(FNAME)
%
% Extract the stimulus matrix from a log file FNAME.log.F ile formats
% are used at the Biophysics department of the Donders Institute. This is
% used for backwards-compatibility.
%
% See also PA_TDT2DAT, PA_READLOG

% 2011 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

%% Initialization
if nargin<1
    fname = [];
end
LogFile    = pa_fcheckexist(fname);
LogFile    = pa_fcheckext(LogFile,'.log');

%% get settings (#trials, #targets, etc..)
Log         = pa_readlog (LogFile);

%% Open file and read stuff
FpLog       = fopen(LogFile,'rt','l');
if FpLog==-1,
    disp(['   Error opening file ' LogFile]);
    return;
end;

Stim        = fgetstim(FpLog,1:Log.NTrial,Log.NTarget);
Nstim       = size(Stim,1);
fclose(FpLog);

%% Complain
if Nstim==0,
    disp(['   No stimuli found in LogFile ' LogFile]);
    return;
end;
if Nstim<Log.NTrial,
    disp([num2str(Log.NTrial-Nstim) ' trials missing in log file ']);
end;


function stim=fgetstim(fpLOG,TrialNr,Ntarget)
% stim=fgetstim(fpLOG,TrialNr,Ntarget)
%   read stimuli from log file. TrialNr is a vector of 
%   trial numbers to be read. Ntarget is the number of
%   targets in a trial. The returned matrix stim is 
%   sorted on trial number. See Index for layout of 
%   stim matrix.
%
%   Jeroen Goossens

TrialNr = sort(TrialNr);
frewind(fpLOG);
i=1;
while ~feof(fpLOG), 
  LineNr = sprintf('#%d',TrialNr(i));
  ident=fscanf(fpLOG,'%s',1);
  if strcmp(ident,LineNr),
    % read line of targets 
    t1=[]; %#ok<NASGU>
    t2=[TrialNr(i),Ntarget];
    for j=1:Ntarget,
      t1=fscanf(fpLOG,'%f',7)';
      t2=[t2, t1(1:3), pa_rphi2azel([t1(2) t1(3)]), t1(4:7)]; %#ok<*AGROW> %PH PATCH
    end;
    stim(i,:)=t2;
    if i==length(TrialNr), 
      return; 
     end;
    i=i+1;
  end;
end;