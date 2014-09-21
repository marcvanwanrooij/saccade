function pa_tdtcal2datcal(fnin,fnout)
%          PA_TDTCAL2DATCAL(FNIN,FNOUT)
%                or
%          PA_TDTCAL2DATCAL(FNIN)
%               or
%          PA_TDTCAL2DATCAL
%
% Converts TDT-files (FNIN, optional) to DAT-files (FNOUT, optional). Data
% and file formats are used in the Biophysics department of the Donders
% Institute. This is only used for backwards-compatibility.
%
% See also PA_TDT2DAT, PA_LOGCAL2CSVCAL

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
tonset =800;
toffset = tonset+99;

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
close all
subplot(131);
plot(X)
pa_verline([tonset toffset]);

subplot(132);
plot(Y)
pa_verline([tonset toffset]);

subplot(133);
plot(Z)
pa_verline([tonset toffset]);

DAT = NaN(size(y));
for ii = 1:ntrials
    indx = tonset:toffset;
    nindx = length(indx);
    indxD = (1:(3*nindx))+(ii-1)*3*nindx;
    DAT(indxD) = [X(indx,ii) Y(indx,ii) Z(indx,ii)];
    figure(2)
    subplot(221)
    plot(X(indx,ii),Y(indx,ii),'.');
    hold on
    
    subplot(222)
    plot(X(indx,ii),Z(indx,ii),'.');
    hold on
        subplot(223)
    plot(Z(indx,ii),Y(indx,ii),'.');
hold on
end
DAT = DAT(~isnan(DAT));

%% Writing DAT and closing
fid1          = fopen(fnout,'wb');
fwrite(fid1,DAT,'float');
fclose(fid1);

