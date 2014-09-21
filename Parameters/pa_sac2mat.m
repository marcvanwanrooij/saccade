function pa_sac2mat(DatFile,LogFile,SacFile)
% PA_SAC2MAT(HVFILE,CSVFILE,SACFILE)
%
% Extract saccade parameters from saccades in data file HVFILE (obtained
% through PA_CALIBRATE) using saccade marks (onset/offset) stored in
% SACFILE (obtained through PA_SACDET) and trial information stored in
% CSVFILE.
%
% PA_SAC2MAT
% 
% Asks for HV-file, and by default sets CSV and SAC filenames to HV
% filename.
%
% See also PA_CALIBRATE, PA_SACDET

% 2013 Marc van Wanrooij
%e: marcvanwanrooij@neural-code.com

%% Initialization
% MinLat                  = -400; % Mininum Latency in ms
if nargin<1
    DatFile             = pa_fcheckexist([],'*.hv');
    LogFile             = pa_fcheckext(DatFile,'csv');
    SacFile             = pa_fcheckext(DatFile,'sac');
end

% [pathstr,name,ext]      = fileparts(DatFile);

%% load log data 
disp(['   Loading from ' LogFile ]);
[expinfo,chaninfo,clog] = pa_readcsv(LogFile);
Stim                    = pa_log2stim(clog);

% get header information
Nsample                 = chaninfo(1,6);
Fsample                 = chaninfo(1,5);
Tsample                 = 1000/Fsample;     % ms
Ntest                   = expinfo(1,4);     % Number of trials
Npresent                = expinfo(1,3);     % number of repeats
Ntrial                  = max(clog(:,1));
if Ntrial ~= Ntest*Npresent
    disp('Error: Number of trials does not correspond to number of tests times number of repeats')
end
Nstim                   = size(Stim,1);
if Nstim == 0,
  disp(['   No stimuli found in LogFile ' LogFile]);
  return;
end;
if Nstim<Ntrial,
  fprintf('   %d trials missing in log file ',Ntrial-Nstim);
end;

% Onset of first target 
% First get rid of Trig0 and Acq "stimuli"
% and look purely for LED, SKY and SND stimuli
seltype                 = ismember(clog(:,5),0:3);
% stimnr                  = unique(clog(seltype,2));
% % and also look for 2nd stimulus
% if size(stimnr,1) == 1
    TarOn                   = Stim(seltype,8); % in ms
% elseif size(stimnr,1)>1
%     seltype                 = seltype & clog(:,2) == stimnr(2);
%     TarOn                   = Stim(seltype,8); % in ms    
% end
%% load sac 

OnOff                       = pa_loadsac(SacFile);
if size(OnOff,2)==0,
  disp(['Error in SAC2MAT: No Saccades found in SacFile ' SacFile]);
  return;
end;
On                          = OnOff(1,:)';
Off                         = OnOff(2,:)';
Nresp                       = 1+fix(max(On)/Nsample);   % number of trials
Nsac                        = length(On);               % number of responses

%% load data
Nchan                   = 2;
ChNr                    = [1,2];
[H,V]                   = pa_fgetsac(DatFile,Nchan,ChNr,SacFile,'Sac');
Ntrace                  = size(H,1);
if Ntrace~=Nsac,
  disp(['   Number of traces does not match number of saccades ' DatFile ]);
  return;
end; 

%% saccade parameters
% truncate matrices if not all trials are included in log file
if Nstim<Nresp,
  S                     = (On<Nstim*Nsample);
  H                     = H(S,:);
  V                     = V(S,:);
  On                    = On(S);
  Off                   = Off(S);
  Ntrial                = Nstim;
  Nsac                  = length(On);
end;
% initialize sac matrix
Sac                     = zeros(Nsac,20);
Sac(:,1)                = 1+fix(On/Nsample);             % 1 trial number
Sac(:,3)                = rem(On,Nsample);               % 3 Onset Time in sample number
Sac(:,4)                = rem(Off,Nsample);              % 4 Offset Time in sample number
TarOn                   = TarOn(Sac(:,1));
Sac(:,5)                = Sac(:,3)*Tsample-TarOn;        % 5 Saccade latency
for i                   = 1:Ntrial,
%   ts                    = (Sac(:,1)==i) & (Sac(:,5)>MinLat);
  ts                    = (Sac(:,1)==i);
  if sum(ts)>0,
    Sac(ts,2)           = (1:sum(ts))';                 % 2 saccade number in trial
  end;
end;
Sac(:,6)                = H(:,1);                       % 6 Hor Onset Position      
Sac(:,7)                = V(:,1);                       % 7 Ver Onset Position      
for i                   = 1:Nsac,
  Sac(i,8)              = H(i,Off(i)-On(i)+1);                       % 8 Hor Offset Position    
  Sac(i,9)              = V(i,Off(i)-On(i)+1);                       % 9 Ver Offset Position
end;
Sac(:,10)               = Sac(:,8)-Sac(:,6);                          % 10 Horizontal Displacement
Sac(:,11)               = Sac(:,9)-Sac(:,7);                          % 11 Vertical Displacement
[R,Phi]                 = pa_azel2pol(Sac(:,10:11));
Sac(:,12)               = R;                     % 12 Saccade Amplitude
Sac(:,13)               = Phi;                     % 13 Saccade direction
Sac(:,14)               = (Off-On).*Tsample;             % 14 Saccade Duration
Sac(:,15)               = 1000*R./Sac(:,14);     % 15 Mean Saccade Velocity    
for i                   = 1:Nsac,
  Rx                    = H(i,isnan(H(i,:))==0)';
  Ry                    = V(i,isnan(V(i,:))==0)';
  Rx                    = Rx-Rx(1);
  Ry                    = Ry-Ry(1);
  [R,Phi]               = pa_azel2pol([Rx,Ry]);
  Vel                   = diff(R) .* Fsample;
  [Vmax,Tmax]           = max(Vel);
  Tmax                  = Tmax(1)*Tsample;
  Sac(i,16)             = Vmax;                                    % 16 Peak Sac Vel.    
  Sac(i,17)             = Tmax;                                    % 17 Time To Peak
  Sac(i,19)             = pa_instdir(R,Phi,Sac(i,12));  % 19 init sac direction
end;
Sac(:,18)               = Sac(:,17)./Sac(:,14);          % 18 Skewness
for i                   = 1:Nsac,
  Sac(i,20)             = pa_curv(H(i,:),V(i,:));         % 20 trajectory curvature 
end;


%% save output in matlab file
MatFile                 = pa_fcheckext(DatFile,'mat');
MatFile                 = MatFile(1:end-4);
disp(['   Writing ' MatFile ]);
save(MatFile,'Stim','Sac')
