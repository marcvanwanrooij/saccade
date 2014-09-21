function pa_hvfilt(DATfiles,OUTfiles, Fcutoff)
% Digital low-pass filtering on HVFILES
%
% HVFILT(HVFILES,OUTFILES,FCUTOFF)
%
% Digital low-pass filtering on HVFILES (with extension HV).
% File names have to be entered as ['AAA';'BBB';'CCC'], with or without the
% '.hv'-extension. If HVFILES are not entered, HVFILT will by default
% filter all hvfiles in the current directory.
% Optional input includes:
% OUTFILES - by default the HVFILES are overwritten!!! If you do
% not want this to occur you have to supply the new name in OUTFILES.
% FCUTOFF - Cut-off frequency (Hz), by default this will amount to 80 Hz.
%
% See also ULTRADET, HVTOMAT, FILTFILT
%
% To do: logfile
%
% Author: Jeroen Goossens & Marc van Wanrooij
% Date: 05-03-2007
%

%% Initialization
if nargin<1
    d                   = dir('*.hv');
    DATfiles            = char(d.name);
end
if nargin <2,
    OUTfiles            = DATfiles;
end;
nchan                   = 2;        % Azimuth channel and elevation channel = 2 channels
h                       = 1;        % Channel number 1 = azimuth
v                       = 2;        % Channel number 2 = elevation

%% Design filter
if nargin<3
    Fcutoff             = 80;       % Cutoff Frequency  (Hz)
end
Order                   = 50;
Fsample                 = 1000;      % Sample Frequency
Wn                      = Fcutoff/(Fsample/2);
N                       = Order;
F                       = [0.0 Wn  Wn  1.0];
M                       = [1.0 1.0 0.0 0.0];
B                       = fir2(N,F,M);
ext                     = '.hv';

%% Looping all files
for i                   = 1 : size(DATfiles,1),
    fname               = DATfiles(i,:);
    fname               = pa_fcheckext(fname,ext);
    csvfile             = pa_fcheckext(fname,'csv');
    [expinfo,chaninfo]  = pa_readcsv(csvfile);
    nsample             = chaninfo(1,6);
    if nsample/3>N
        disp(['   Reading ' fname]);
        fname               = pa_fcheckexist(fname);
        fid                 = fopen(fname,'r','l');
        if fid==-1,
            disp(['   Error reading ' fname ]);
            return;
        end;
        [mtx,n]             = fread(fid,[nchan,inf],'float');
        fclose(fid);
        % filter the signals trial by trial
        disp(['   Low-Pass Fc = ' int2str(Fcutoff) ' Hz   N = ' int2str(N)] );
        nblock      = n/(nsample*nchan);
        for j       = 1:nblock,
            indx        = (j-1)*nsample+1:j*nsample;
            fs          = filtfilt(B,1,mtx(h,indx));
            mtx(h,indx) = fs;
            fs          = filtfilt(B,1,mtx(v,indx));
            mtx(v,indx) = fs;
        end;
        % save filtered data
        fname       = OUTfiles(i,:);
        fname       = pa_fcheckext(fname,ext);
        disp(['   Writing ' fname]);
        fid         = fopen(fname,'w','l');
        if fid==-1,
            disp(['   Error writing ' fname ]);
            return;
        end;
        MTX         = [mtx(h,:); mtx(v,:)];
        MTX         = MTX(:);
        fwrite(fid,MTX,'float');
        fclose(fid);
    elseif nsample/3<N
        disp(['skip ' fname])
        disp('Data is not filtered due to small number of samples');
%         warning('HVfilt:NotEnoughSamples','Data is not filtered due to small number of samples');
    end
end;
