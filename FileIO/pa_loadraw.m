function varargout          = pa_loadraw(fname,nchan,nsample)
% Load data from raw data file (DAT, HV)
%
% VARARGOUT = LOADRAW(FNAME,NCHAN,NSAMPLE)
%   Load a raw data file (.hum, .aap, .hv, .hvn) in which
%   NCHAN channels with NSAMPLE datapoints are stored interleaved. 
%   Returns a matrix for each channel (VARARGOUT) in which each trial is stored row-wise.
%
% See also ULTRADET, GETVEL
%
% Adaptation by: MarcW


if nargin<1
    fname = [];
end

%% find file type
fname                                   = pa_fcheckexist(fname);
[pathstr,name,ext]                      = fileparts(fname); %#ok<*ASGLU>
ext                                     = upper(ext);
Nbyte                                   = 0;
if strcmp(ext,'.HUM'), Nbyte             = 2;  end;
if strcmp(ext,'.AAP'), Nbyte             = 2;  end;
if strcmp(ext,'.HV'),  Nbyte             = 4;  end;
if strcmp(ext,'.HVN'), Nbyte             = 4;  end;
if Nbyte==0,
    disp(['invalid fname extension ',ext]);
    return;
end;

%% read file
fid                                     = fopen(fname,'r','l');
if fid==-1,
    disp(['   Error opening ' fname ]);
    return;
end;
if Nbyte==2,
    mtx                                 = fread(fid,[nchan*nsample,inf],'ushort');
end;
if Nbyte==4,
    mtx                                 = fread(fid,[nchan*nsample,inf],'float');
end;
fclose(fid);


%% separate channels
nout                                    = nargout;
if nout~=nchan
    error('Number of outputs does not equal number of channels');
end
varargout                               = cell(nout,1);
for i                                   = 1:nchan,
    ch                                  = mtx(i:nchan:nchan*nsample,:);
    varargout(i)                        = {ch};
end;
