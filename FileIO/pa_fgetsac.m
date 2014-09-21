function [ch1,ch2,ch3,ch4,ch5,...
    ch6,ch7,ch8,ch9,ch10]...
    =pa_fgetsac(DATfile,Nchan,ChNr,SACfile,Type,Pre,Post,Dur) %#ok<STOUT>

% Load saccade data from a data file
% [CH1, CH2, ..., CHN] = FGETSAC(DATFILE,NCHAN,CHNR,SACFILE,TYPE,PRE,POST,DUR)
%   load saccade data from a data file (DAT,AAP,HV,HVN) using a sac file.
%   Nchan is the number of channels in the file and ChNr provide the
%   indices to the desired channels.
%
%   Pre, Post, Dur and  Type are optional:
%     Type = Reading type ('Fix' or 'Sac') default 'Sac'
%     Pre  = Nsamples prior to saccade onset, default 0
%     Post = Nsamples after saccade offset, default 0
%     Dur  = Nsamples of saccade, default max(Off-On+1)
%
%   Samples that lie outside the (On-Pre)-(Off+Post) interval are not
%   read if Type='Sac' the remaining samples are set to NaN.
%
%   Original: Jeroen Goossens
%   Adaptation: MW


%% check parameters
if nargin<5, Type ='SAC';
else Type = upper(Type ); end
if nargin<6, Pre  =0; end
if nargin<7, Post =0; end
if nargin<8, Dur  =0; end
% nout    = nargout;

%% find file type
[pathname fname ext]    = fileparts(DATfile);

Nbyte                   = 0;
if strcmpi(ext,'.hum'), Nbyte=2; end
if strcmpi(ext,'.dat'), Nbyte=4; end
if strcmpi(ext,'.aap'), Nbyte=2; end
if strcmpi(ext,'.hv'),  Nbyte=4; end
if strcmpi(ext,'.hvn'), Nbyte=4; end
if Nbyte==0,
    disp(['   Error : Invalid DATfile extension ',ext]);
    return;
end

%% load saccade marks
sac         = pa_loadsac(SACfile);
onset       = sac(1,:)-Pre;
offset      = sac(2,:)+Post;
Nsac        = length(onset);
Npnt        = max(offset-onset)+1;

if Dur~=0,
    if Dur<Npnt,
        disp('   Warning : Duration not long enough to include all data');
    end;
    Npnt      = Dur;
end;

%% set number of samples to read for each saccade
if strcmpi(Type,'sac')
    Nsample   = offset-onset+1;
end
if strcmpi(Type,'fix')
    Nsample   = Npnt*ones(Nsac,1);
end

%% open dos data file
disp(['   Loading saccades from ' DATfile]);

fid         = fopen(DATfile,'r','l');
if fid==-1,
    disp(['   Error : Unable to open ' DATfile]);
    return;
end;

%% init ch matrices
for i       = 1:length(ChNr),
    cmd     = [ 'ch' num2str(i) ' = NaN*ones(Nsac,Npnt);' ];
    eval(cmd);
end;

%% read saccades
for i       = 1:Nsac,
    % set file pointer at beginning of saccade
    fpos      = Nbyte*Nchan*onset(i);
    status    = fseek(fid,fpos,'bof');
    if status ~= 0 ,
        disp(['   Error : Invalid file position ' DATfile ': Sac #' num2str(i) ]);
        return;
    end;

    % read saccade
    switch Nbyte
        case 2,
            [mtx,N] = fread(fid,[Nchan,Nsample(i)],'ushort');
        case 4,
            [mtx,N] = fread(fid,[Nchan,Nsample(i)],'float');
    end;
    if N~=Nchan*Nsample(i),
        disp(['   Warning : Saccade ' num2str(i) ' truncated ']);
    end;

    % extract channels
    ns        = size(mtx,2); %#ok<NASGU>
    for j     = 1:length(ChNr),
        nr      = ChNr(j);  %#ok<NASGU>
        cmd     = [ 'ch' num2str(j) '(i,1:ns) = mtx(nr,:); ' ];
        eval(cmd);
    end;

end;

%% END
fclose(fid);
