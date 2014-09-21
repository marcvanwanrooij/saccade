function pa_dmicalfiles(DatFiles,CalFile)
% Calibrate all raw data contained in one or several files
%
% CALIBRATE(<DATFILES>,<CALFILE>)
%
%  Function to calibrate raw data into azimuth/elevation angles. These
%  angles are stored in hv-files.
%
%
%
%       DATFILES:       Raw data. 1 or more files in rows.
%                   eg. ['XX-XX-2000-01-0101<.dat>'; ...
%                        'XX-XX-2000-01-0102<.dat>'; ...
%                        'XX-XX-2000-01-0103<.dat>'];
%                       default: all dat-files in current directory
%
%       CALFILE:        Neural network file.
%                   eg. 'XX-XX-2000-01-01<.net>';
%                       default: user input
%
%
%       output:         hv-files named according to the DatFiles.
%                   eg. 'XX-XX-2000-01-0101.hv' ... etc.
%
%  See also ULTRADET, TRAINCAL
%
%  Author: Marcus
%  Date: 11-04-07


%% Initialization
if nargin<1
    d                        = dir('*.dat');
    DatFiles                = char(d.name);
end
if nargin<2
    d                        = dir('*.net');
    CalFile                = char(d.name);
    CalFile                 = pa_fcheckexist(CalFile,'*.net');
end

Hchan                       = 5;
Vchan                       = 6;
Fchan                       = 7;
Hheadchan                       = 1;
Vheadchan                       = 2;
Fheadchan                       = 3;
S                           = load(CalFile,'-mat');

%% Calibrate all DATfiles
for i                       = 1:size(DatFiles,1),
    % Loading file
    fname                   = pa_fcheckext(DatFiles(i,:),'.dat');
    fname                   = pa_fcheckexist(fname,'*.dat');
    csvname                 = pa_fcheckext(DatFiles(i,:),'.csv');
    [expinfo,chaninfo]      = pa_readcsv(csvname);
    nchan                   = expinfo(1,8);
    nsample                 = chaninfo(1,6);
    DAT                     = pa_loaddat(fname,nchan,nsample);
    H                       = squeeze(DAT(:,:,Hchan));
    H                       = H(:);
    V                       = squeeze(DAT(:,:,Vchan));
    V                       = V(:);
    F                       = squeeze(DAT(:,:,Fchan));
    F                       = F(:);
    Hhead                       = squeeze(DAT(:,:,Hheadchan));
    Hhead                       = Hhead(:);
    Vhead                       = squeeze(DAT(:,:,Vheadchan));
    Vhead                       = V(:);
    Fhead                       = squeeze(DAT(:,:,Fheadchan));
    Fhead                       = Fhead(:);
    DAT                     = [H V F Hhead Vhead Fhead]';
    [AZ,EL]                 = pa_calib(DAT,S);

    % Saving calibrated data
    fname                   = pa_fcheckext(DatFiles(i,:),'.hv');
    fid                     = fopen(fname,'w','l');
    AZEL                    = [AZ;EL];
    fwrite(fid,AZEL,'float');
    fclose(fid);
end


