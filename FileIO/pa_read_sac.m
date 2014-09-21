function [DAT]=pa_read_sac(datfile,nchan,nsample,ntrialsin)
% Load data from DAT-files
%
% function [DAT]=loaddat(datfile,nchan,nsample,<ntrials>)
%
%   Reads data from DATAFILE, based on nchan and nsample
%
%
%       datfile =   Datafile
%                       (e.g. XX0101<.dat>)
%       nchan =     Number of channels
%                       (e.g. HEAD x, HEAD y, EYE x & EYE y, gives
%                       nchan = 4;)
%       nsample =   Number of samples
%                       nsample = trail duration (sec) * samplerate
%                       (e.g. trail duration is 2 sec. and the samplerate
%                       is 1000Hz, than nsample = 2000;)
%       ntrials =   Number of trials. If given loaddat checks if number
%                   of bytes corresponds to this number of trials.
%
%   See also: READCSV
%
% tomg Oct 2006
% Modified by: MarcW 2007

%% Initialization
nin             = 4;
if nargin<nin
    ntrialsin   = NaN;
end
ext             = '.dat';
datfile         = fcheckext(datfile,ext);
datfile         = fcheckexist(datfile);
nbytesperfloat  = 4;




%% Check for number of trials
fid             = fopen(datfile,'r');
fseek(fid,0,'eof');
nbytes          = ftell(fid);
ntrials         = nbytes./(nchan*nsample*nbytesperfloat);

% error messages for number of trials
if mod(ntrials,floor(ntrials))
    disp(['     Number of trials: ' sprintf('%8.8f',ntrials)])
    warning('LOADDAT:Integer','Number of trials is not integer');
elseif ~isnan(ntrialsin)
    if ntrialsin ~= ntrials %#ok<BDSCI>
        disp(['     Number of trials: ' sprintf('%8.8f',ntrials) ' ~= ' sprintf('%8.8f',ntrialsin) ])
        error('Number of trials is not equal to input');
    end
end

%% Read Data
frewind(fid)
if ispc
    DAT         = fread(fid,inf,'float');
else
    DAT         = fread(fid,inf,'float','l');
end
DAT             = reshape(DAT,nsample,nchan,ntrials);
DAT             = permute(DAT,[1 3 2]);
fclose(fid);

