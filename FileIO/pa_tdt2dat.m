function pa_tdt2dat(fnin,fnout)
%          PA_TDT2HUM(FNIN,FNOUT)
%                or
%          PA_TDT2HUM(FNIN)
%               or
%          PA_TDT2HUM
%
% Converts TDT-files (FNIN, optional) to DAT-files (FNOUT, optional). Data
% and file formats are used in the Biophysics department of the Donders
% Institute. This is only used for backwards-compatibility.
%
% See also PA_LOG2STIM, PA_READSTIM, PA_READLOG

% 2011 Marc van Wanrooij
% e-mail:marcvanwanrooij@neural-code.com

%% Initialization

if nargin==0
    fnin = pa_fcheckexist([],'*.tdt');
    fnout = pa_fcheckext(fnin,'dat');
    makehum(fnin,fnout);
elseif nargin==1
    fnin = pa_fcheckext(fnin,'tdt');
    fnout = pa_fcheckext(fnin,'dat');
    makehum(fnin,fnout);
elseif nargin == 2
    fnin = pa_fcheckext(fnin,'tdt');
    fnout = pa_fcheckext(fnout,'dat');
    makehum(fnin,fnout);
end


function makehum(fnin,fnout)
disp(['  Converting: ' fnin ' to ' fnout ]);

%% Opening TDT file
fid           = fopen(fnin,'r');
y             = fread(fid,'short');
fclose(fid);

%% Convert TDT data to DAT data
nsample = 1000;
X       = y(1:3:end);
ntrials = round(length(X)/nsample);
Y       = y(2:3:end);
Z       = y(3:3:end);
X       = reshape(X,nsample,ntrials);
Y       = reshape(Y,nsample,ntrials);
Z       = reshape(Z,nsample,ntrials);
DAT = NaN(size(y));
for ii = 1:ntrials
    indx = 1:nsample;
    indxD = (1:(3*nsample))+(ii-1)*3*nsample;
    DAT(indxD) = [X(indx,ii) Y(indx,ii) Z(indx,ii)];
end

%% Writing DAT and closing
fid1          = fopen(fnout,'wb');
fwrite(fid1,DAT,'float');
fclose(fid1);

