function varargout = pa_sacdet(varargin)
% PA_SACDET(HVFILE,DET)
%
% Saccade detection interface for Matlab.
% HVFILE - file containing calibrated horizontal and vertical eye/head
% movement traces.
% DET - array containing detection parameters:
%           Structure       Eenheid     Default
%           det.velocityon    (deg/s)     100
%           det.amplitude   (deg)       2
%           det.duration    (ms)        10
%           det.start       (ms)        [] = 80 ms after first stimulus
%           det.end         (ms)        [] = 1000 ms after last stimulus
%           det.smooth      (s)          0.01
%
% See also HVFILT, LOADSAC, HVTOMAT
%
% Author: Marc van Wanrooij
% Date: 1 April 2007
% Version: Matlab R2006b
%
%   OnOff-MTX(5,N) consists of:
%   1 - Onset
%   2 - Offset
%   3 - Trialnr
%   4 - SacN
%   5 - Checked
%
% Depends on:
% See also:
% LOADSAC
% LOADRAW
% pa_readcsv
% GETVEL
% AZEL2POL, AZEL2CART, CART2POL, FART2AZEL
% pa_fcheckext
% pa_fcheckexist
% pa_log2stim
%
%
%

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to menu_help pa_sacdet

% Last Modified by GUIDE v2.5 04-Sep-2011 15:16:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State      = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @pa_sacdet_OpeningFcn, ...
    'gui_OutputFcn',  @pa_sacdet_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before pa_sacdet is made visible.
function pa_sacdet_OpeningFcn(hObject, eventdata, handles, varargin) 
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pa_sacdet (see VARARGIN)
home;
disp('>>PA_SACDET <<');

% Get Data
handles                                = MAIN_Check_And_Load(handles, varargin{:});
handles                                = MAIN_detsaccade(handles);
handles                                = MAIN_detsacfile(handles);
handles                                = MAIN_plottrial(handles);
% sacsave(handles);
% Choose default command line output for pa_sacdet
handles.output                              = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pa_sacdet wait for user response (see UIRESUME)
% uiwait(handles.pa_sacdet);


% --- Outputs from this function are returned to the command line.
function varargout = pa_sacdet_OutputFcn(hObject, eventdata, handles) %#ok<*INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1}                                = handles.output;

%% MAIN_Check_And_Load(HANDLES,VARARGIN)
function handles                          = MAIN_Check_And_Load(handles, varargin)

%% Initialization
lv                                          = length(varargin);
if lv<1
    handles.fname                               = [];
    handles.fname                               = pa_fcheckexist(handles.fname,'*.hv');
    det.velocityon                                = 10;       %deg/s
    det.velocityoff                             = 10;      %deg/s
    det.smooth                                  = 0.01;
    det.duration                                = 10;       % ms
    det.amplitude                               = 2;        % deg
    det.start                                   = 1500;       % ms, starting-time of detection
    det.end                                     = 2300;       % ms, end-time of detection
    det.acc                                     = 0;
end

if length(varargin)==1
    handles.fname                               = varargin{1};
    handles.fname                               = pa_fcheckext(handles.fname,'.hv');
    handles.fname                               = pa_fcheckexist(handles.fname);
    det.velocityon                                = 20;       % deg/s
    det.velocityoff                             = 10;      % deg/s
    det.smooth                                  = 0.01;
    det.duration                                = 20;       % ms
    det.amplitude                               = 4;        % deg
    det.start                                   = 80;       % ms, starting-time of detection
    det.end                                     = 800;       % ms, end-time of detection
    det.acc                                     = 0;
end
if length(varargin)==2
    handles.fname                               = varargin{1};
    handles.fname                               = pa_fcheckext(handles.fname,'.hv');
    handles.fname                               = pa_fcheckexist(handles.fname);
    det                                         = varargin{2};
end

%% History Log
str1                                            = get(handles.display,'String');
str2                                            = ['   Loading ' upper(handles.fname)];
str                                             = char({str1;str2});
set(handles.display,'String',str);
str1                                            = get(handles.display,'String');
str2                                            = ['   Loading ' upper([handles.fname(1:end-3) '.log'])];
str                                             = char({str1;str2});
set(handles.display,'String',str);

%% Trial, stimulus and channel information
[expinfo,chaninfo,mLog]                         = pa_readcsv([handles.fname(1:end-3) '.csv']); %#ok<*ASGLU>
handles.Nsamples                                = chaninfo(1,6);    % # of samples
handles.Fsample                                 = chaninfo(1,5);    % Hz (samples/s)
handles.Ntrial                                  = max(mLog(:,1));   % #
handles.StimType                                = mLog(:,5);        % stimulus type
handles.StimOnset                               = mLog(:,8);        % msec
handles.Stim                                    = pa_log2stim(mLog);   % Stim-matrix

% Adjusting detection time
% SRT>80 and <1000 ms, else listener predictive and or lazy
if length(varargin)<2
%     sel                                         = handles.StimType == 1 | handles.StimType == 0 | handles.StimType == 4;
    sel = ismember(handles.StimType,[0 1 4 5]);
    if floor(mod(length(sel(sel))./handles.Ntrial,1))
        handles.Onset                           = reshape(handles.StimOnset(sel),length(handles.StimOnset(sel))/handles.Ntrial,handles.Ntrial);

        if isempty(det.end)
            det.end                             = max(max(handles.Onset))+600;
        end
        if isempty(det.start)
            sel                                 = handles.Onset>0;
            det.start                           = min(min(handles.Onset(sel)))+80;
        end
        if det.end > 1000*handles.Nsamples/handles.Fsample;
            det.end                             = round(1000*handles.Nsamples/handles.Fsample);
        end
    else
        det.start                               = 100;
        if isempty(det.end)        
            det.end                             = 1800;
        end
    end
end

%% Eye/Head/Pointer Traces
[handles.htrace, handles.vtrace]                 = pa_loadraw(handles.fname, 2, handles.Nsamples);
[handles.veltrace, handles.smvtrace, handles.acctrace, handles.smatrace] = getvel(handles.htrace,handles.vtrace,handles.Fsample,det.smooth);


%% The handle structure
% all the information (except for On and Offset) you will ever need
handles.det                                     = det;
handles.ctrial                                  = 1;
handles.csaccade                                = 1;
handles.del                                     = [];
handles.axset                                   = 'a'; % Axis Setting: d - FixeD [-80 80], m - Manual, a - Automatic
set(handles.menu_axauto,'Checked','on');
set(handles.menu_axdef,'Checked','off');
set(handles.menu_axmanual,'Checked','off');
handles.extraaxes                               = 'm'; % Radial plot - r, Main Sequence - m
handles.mfigcheck                               = false;
handles.rfigcheck                               = true;
set(handles.menu_opt_radial_plot,'Checked','off');
set(handles.menu_opt_mainseq_axes,'Checked','on');
handles.stimulusshown                           = 'b'; % Blind detection - d, Visible Stimuli - v
set(handles.menu_blind,'Checked','on');
set(handles.menu_visible,'Checked','off');
if handles.det.acc==1
    set(handles.chck_acc,'Value',1);
elseif  handles.det.acc==0
    set(handles.chck_acc,'Value',0);
end
handles.redetect = 'sub'; % Redetect current and subsequent trials - sub; current trial - cur; all trials - all

function handles                                = MAIN_detsaccade(handles)
% Detection based on smoothed Velocity
% then checked for amplitude and duration of saccade
% and finally finetune detection by looking for 0 acceleration
% Initialization
SacOn                                           = NaN*ones(1,handles.Ntrial*3); % Saccade Onset
SacOff                                          = SacOn;    % Saccade Offset
Trial                                           = SacOn;     % Trial Number
SacN                                            = SacOn;    % Saccade Number in Trial
k                                               = 0;

% Find Velocity
for i                                           = 1:size(handles.smvtrace,2);
    V                                           = handles.smvtrace(:,i);
    HOR                                         = handles.htrace(:,i);
    VER                                         = handles.vtrace(:,i);
    %     % Detect Velocity
    %     sel                                         = (V>handles.det.velocityon);  % Trace Index in trial where saccade velocity exceeds minimum velocity
    %     onoff                                       = zeros(size(V));        % If saccade velocity does not exceed velocity -> value = 0
    %     onoff(sel)                                  = ones(size(sel(sel)));  % if it does exceed velocity -> value = 1
    %     onoff                                       = [diff(onoff);-1];      % Determine Difference between succeeding samples to obtain:
    %     on                                          = find(onoff>0);         % onsets (+1)
    %     off                                         = find(onoff<0);         %  and offsets (-1)

    sel1                                        = (V>handles.det.velocityon);  	% Trace Index in trial where saccade velocity exceeds minimum velocity
    sel2                                        = (V>handles.det.velocityoff);  % Trace Index in trial where saccade velocity exceeds minimum velocity
    tmpon                                       = find([0;diff(sel1)]==1); % onsets (+1)
    on                                          = tmpon(tmpon~=handles.Nsamples); % QUE?
    off                                         = NaN(size(on));
    tmpoff                                      = [find([0;diff(sel2)]==-1);numel(V)];                %  and offsets (-1)
    for I_on = 1:numel(on)
        off(I_on) = tmpoff(find(tmpoff-on(I_on)>0,1,'first'));
    end

    % Check if first saccade offset is not due to slowness of subject
    if ~isempty(on) && ~isempty(off)
        if off(1)<on(1);
            off                         = off(2:end);
        end
        % Check ending - End of trial is detected as offset of saccade
        if length(on)~=length(off);
            if round(off(end)) == handles.Nsamples;
                off                     = off(1:end-1);
            end
        end

        % Check start & ending of detection
        sel                            = (on>handles.det.start./1000.*handles.Fsample) & (off>handles.det.start./1000.*handles.Fsample);
        on                             = on(sel);
        off                            = off(sel);

        if ~isempty(on) && ~isempty(off)
            sel                        =(off>on(1));
            off                        = off(sel);
            sel                        = (on<handles.det.end/1000.*handles.Fsample);
            on                         = on(sel);
            off                        = off(sel);
            % Check Duration
            sel                         = ((off-on)./handles.Fsample>handles.det.duration/1000);
            on                          = on(sel);
            off                         = off(sel);
            % check Amplitude
            if ~isempty(on)
                ontmp                   = zeros(size(on));
                offtmp                  = zeros(size(off));
                for j                   = 1:length(on)
                    Azon                = HOR(on(j));
                    Elon                = VER(on(j));
                    Azoff               = HOR(off(j));
                    Eloff               = VER(off(j));
                    AZ                  = Azoff-Azon;
                    EL                  = Eloff-Elon;
                    R                   = pa_azel2pol(AZ,EL);
                    if R > handles.det.amplitude;
                        ontmp(j)        = on(j);
                        offtmp(j)       = off(j);
                    end
                end
                indx                   = find(ontmp);
                on                     = on(indx);
                off                    = off(indx);
            end
            % Combining Trials
            on                         = on';
            off                        = off';
            % Check End - Offset should not equal # samples. Why?
            sel			= off==handles.Nsamples;
            off(sel)	= handles.Nsamples-1;
            % Continue
            lOn                        = length(on);
            Trialtmp                   = ones(1,lOn).*i;
            NSac                       = (1:lOn);
            if ~isempty(NSac)
                SacOn(:,NSac+k)        = on;
                SacOff(:,NSac+k)       = off;
                Trial(:,NSac+k)        = Trialtmp;
                SacN(:,NSac+k)         = NSac;
                k                      = k+lOn;
            end
        end
    end
end
%% Remove NaN
sel                                     = ~isnan(SacOn);
SacOn                                   = SacOn(sel);
SacOff                                  = SacOff(sel);
Trial                                   = Trial(sel);
SacN                                    = SacN(sel);

%% Remove "double saccades"
% when onset of one saccade falls inbetween the on- and offset of another

uTrial = unique(Trial);
for i = 1:length(uTrial)
	indx = find(Trial==uTrial(i));
	if length(indx)>1
		A = [Trial(indx); SacOn(indx); SacOff(indx)];
		for j = 2:length(indx)
			if A(2,j)<A(3,j-1)
				Trial(indx(j)) = NaN;
				SacOn(indx(j)) = NaN;
				SacOff(indx(j)) = NaN;
			end
		end
	end
end
%% Remove NaN
sel                                     = ~isnan(SacOn);
SacOn                                   = SacOn(sel);
SacOff                                  = SacOff(sel);
Trial                                   = Trial(sel);
SacN                                    = SacN(sel);

% And now do some finetuning
% Acceleration can be negative and positive, and also 0 (radial velocity
% cannot). So search for 0 acceleration!

% And finally a check to see whether the saccades have been visually
% checked - 0 = no check, 1 = check
if ~isempty(SacOn)
    chck                                    = zeros(size(SacOn));
    chck(1)                                 = 1;
    OnOff                                   = [SacOn ; SacOff ;Trial ;SacN ;chck];
    if handles.det.acc == 1
        for i                                   = 1:size(OnOff,2)
            k                                   = OnOff(1,i);
            Acc                                 = handles.smatrace(k,OnOff(3,i));
            while Acc > 0
                k                               = k-1;
                if k>0
                    Acc                         = handles.smatrace(k,OnOff(3,i));
                else
                    Acc                         = 0;
                    k                           = k+1;
                end
            end
            OnOff(1,i)                          = k;
        end
        for i                                   = 1:size(OnOff,2)
            k                                   = OnOff(2,i);
            Acc                                 = handles.smatrace(k,OnOff(3,i));
            while Acc < 0
                k                               = k+1;
                if k<=handles.Nsamples
                    Acc                         = handles.smatrace(k,OnOff(3,i));
                else
                    Acc                         = 0;
                    k                           = k-1;
                end
            end
            OnOff(2,i)                          = k;
        end
    end
    sel = OnOff(2,:)==handles.Nsamples;
    OnOff(2,sel) = handles.Nsamples-1;
    handles.onoff = OnOff;
else
    handles.onoff = [];
end


% Check BugWare
if size(on)~=size(off)
    disp('Uh-0oh - number of saccade-starts does not match number of saccade-ends');
    return
end
handles.onoffdet                        = repmat([handles.det.velocityon handles.det.amplitude handles.det.duration handles.det.smooth handles.det.start handles.det.end],size(handles.onoff,2),1);
handles.onoffdet                        = handles.onoffdet';

function handles                        = MAIN_detsacfile(handles)
sacfile                                 = pa_fcheckext(handles.fname,'.sac');
if exist(sacfile,'file')
    OnOffExist                          = pa_loadsac(sacfile);
    str1                                = get(handles.display,'String');
    str2                                = ['Loading ' upper(sacfile)];
    str                                 = char({str1;str2});
    set(handles.display,'String',str);
    OnExist                             = OnOffExist(1,:);
    TrialExist                          = floor(OnExist./handles.Nsamples)+1;
    OnExist                             = OnExist - ((TrialExist-1).*handles.Nsamples);
    OffExist                            = OnOffExist(2,:);
    OffExist                            = OffExist - ((TrialExist-1).*handles.Nsamples);
    SacNrExist                          = my_unique(TrialExist)';
    CheckExist                          = ones(size(OnExist));
    OnOffExist                          = [OnExist; OffExist; TrialExist; SacNrExist; CheckExist];
    trialn                              = max(TrialExist);
    indx                                    = find(handles.onoff(3,:)==trialn); %#ok<MXFND>
    indx                                    = max(indx);
    handles.onoff                           = [OnOffExist handles.onoff(:,indx+1:end)];
    handles.ctrial                          = trialn;
    handles.csaccade                   = OnOffExist(4,end);
end

function handles                       = MAIN_plottrial(handles)
linkaxes([handles.axes_trace handles.axes_velocity],'x');
%% Get Variables from Handles Structure
ctrial                                 = handles.ctrial;
csaccade                               = handles.csaccade;
det                                    = handles.det;
trial                                  = handles.onoff(3,:);
on                                     = handles.onoff(1,:);
off                                    = handles.onoff(2,:);
sacnr                                  = handles.onoff(4,:);
% checksac                               = handles.onoff(5,:);
%% Get Current Trial Variables
htrace                                 = handles.htrace(:,handles.ctrial);
vtrace                                 = handles.vtrace(:,handles.ctrial);
veltrace                               = handles.veltrace(:,handles.ctrial);
smvtrace                               = handles.smvtrace(:,handles.ctrial);
smatrace                               = handles.smatrace(:,handles.ctrial);
TimeTrace                              = (1:length(htrace))./handles.Fsample;

%% Setting all Strings in GUI to corresponding values
[pstr,fname]                           = fileparts(handles.fname);
set(handles.fname_txt,'String',['File Name: ' upper(fname)]);
set(handles.trialnr_txt,'String',['Trial: ' num2str(handles.ctrial) ' of ' num2str(handles.Ntrial)]);
set(handles.txt_jumpsac,'String',num2str(handles.csaccade));
set(handles.txt_jumptrial,'String',num2str(handles.ctrial));
set(handles.sacnr_txt,'String',['Saccade # in Trial: ' num2str(handles.csaccade)]);
set(handles.totalsac_nr,'String',['Saccade # in Experiment: ' num2str(handles.csaccade)]);
set(handles.txt_vel,'String',num2str(handles.det.velocityon));
set(handles.txt_veloff,'String',num2str(handles.det.velocityoff));
set(handles.txt_amp,'String',num2str(handles.det.amplitude));
set(handles.txt_duration,'String',num2str(handles.det.duration));
set(handles.txt_smooth,'String',num2str(handles.det.smooth));
set(handles.txt_detstart,'String',num2str(handles.det.start));
set(handles.txt_detend,'String',num2str(handles.det.end));
set(handles.axman_txt,'Visible','off');
set(handles.maxax_txt,'Visible','off');
set(handles.minax_txt,'Visible','off');

%% And get Current Saccade Variables
indxsac                                  = find(trial == handles.ctrial & sacnr == handles.csaccade);
if ~isempty(indxsac)
    curOn                                = on(indxsac);             % index in sample numbers
    curOff                               = off(indxsac);            % index in sample numbers
    curOnTime                            = on(indxsac)/handles.Fsample;     % time in sec
    curOffTime                           = off(indxsac)/handles.Fsample;    % time in sec
    curHTrace                            = htrace(curOn:curOff);
    curVTrace                            = vtrace(curOn:curOff);
    curTTrace                            = (curOn:curOff)./handles.Fsample;
    cursmvtrace                          = smvtrace(curOn:curOff);
    cursmatrace                          = smatrace(curOn:curOff);
    curveltrace                          = veltrace(curOn:curOff);
    indx                                 = find(trial==ctrial);
    set(handles.sacnr_txt,'String',['Saccade # in Trial: ' num2str(csaccade) ' of ' num2str(sacnr(indx(end)))]);
    set(handles.totalsac_nr,'String',['Saccade # in Experiment: ' num2str(indxsac) ' of ' num2str(length(on))]);
elseif isempty(indxsac)
    set(handles.sacnr_txt,'String','Saccade # in Trial: 0 of 0');
    set(handles.totalsac_nr,'String',['Saccade # in Experiment: 0 of ' num2str(length(on))]);
end

%% Some Graphics Globals
indxonoff                                = find(trial==ctrial);

%% Graphics: Horizontal and Vertical Traces
axes(handles.axes_trace); %#ok<*NASGU>
cla;
% zoom on;
plot(TimeTrace,htrace,'b-');hold on;
plot(TimeTrace,vtrace,'r-');
% Axes
if strcmp(handles.axset,'d');
    ax                                   = axis;
    if ax(3)>=-80; ax(3)                 = -80; end;
    if ax(4)<=80; ax(4)                  = 80; end; axis(ax);
elseif strcmp(handles.axset,'a');
    axis auto;
    ax                                   = axis;
    %     sel                                  = TimeTrace>det.start/1000 & TimeTrace < det.end/1000;
    ax(4)                                = ceil(max([htrace;vtrace])/10)*10;
    ax(3)                                = floor(min([htrace;vtrace])/10)*10;
    axis(ax)
elseif strcmp(handles.axset,'m');
    set(handles.axman_txt,'Visible','on');
    set(handles.maxax_txt,'Visible','on');
    set(handles.minax_txt,'Visible','on');
    axis auto;
    ax                                   = axis;
    maxax                                = get(handles.maxax_txt,'String');
    ax(4)                                = str2double(maxax);
    minax                                = get(handles.minax_txt,'String');
    ax(3)                                = str2double(minax);
    axis(ax);
    ax                                   = axis;
end

%% Time periods in which no saccades are detected
xpatch                                   = [ax(1) ax(1) det.start/1000 det.start/1000];
ypatch                                   = [ax(3) ax(4) ax(4) ax(3)];
hp                                       = patch(xpatch,ypatch,'r'); set(hp,'FaceColor',[230 230 255]./255,'EdgeColor','b');
xpatch                                   = [det.end/1000 det.end/1000 ax(2) ax(2)];
ypatch                                   = [ax(3) ax(4) ax(4) ax(3)];
hp                                       = patch(xpatch,ypatch,'r'); set(hp,'FaceColor',[230 230 255]./255,'EdgeColor','b');
% 0-deg line
horline(0,'k-');
% Plotting traces over the previous patches
plot(TimeTrace,htrace,'b-');
plot(TimeTrace,vtrace,'r-');

% Plot current saccade if present
if ~isempty(indxsac)
    % including the time period of the current saccade
    xpatch                               = [curOnTime curOnTime curOffTime curOffTime];
    ypatch                               = [ax(3) ax(4) ax(4) ax(3)];
    hp                                   = patch(xpatch,ypatch,'g');
    set(hp,'FaceColor',[230 255 230]./255,'EdgeColor','g');
    hp                                   = plot(curTTrace,curHTrace,'b-'); set(hp,'LineWidth',3);
    hp                                   = plot(curTTrace,curVTrace,'r-'); set(hp,'LineWidth',3);
end;
ylabel('P');
% Plotting on- and offsets of all saccades in current trial
for i                                    = 1:length(indxonoff);
    pa_verline(on(indxonoff(i))./handles.Fsample,'g');
    pa_verline(off(indxonoff(i))./handles.Fsample,'g');
end
% 'Boldly' plotting on- and offset of current saccade if present
if ~isempty(indxsac)
    for i                                = 1:length(indxsac);
        h                                = pa_verline(on(indxsac(i))./handles.Fsample,'g');set(h,'LineWidth',2);
        h                                = pa_verline(off(indxsac(i))./handles.Fsample,'g');set(h,'LineWidth',2);
    end
end
% Stimulus Parameters
if strcmp(handles.stimulusshown,'v')
    sel                                 = handles.Stim(:,1) == ctrial;
    %     StimNumber      = handles.Stim(sel,2);
    StimOn                              = handles.Stim(sel,8);
    StimOn(StimOn<0)                    = 0;
    StimOff                             = handles.Stim(sel,9);
    StimOff(isnan(StimOff))             = StimOn(isnan(StimOff))+100;
    StimAz                              = handles.Stim(sel,4);
    StimEl                              = handles.Stim(sel,5);
    %     for i=1:length(StimNumber) % 1st stimulus = fixation spotlight
    for i=2 % 1st stimulus = fixation spotlight
        x1                              = StimOn(i);
        x2                              = StimOff(i);
        yaz                             = StimAz(i);
        yel                             = StimEl(i);
        horline(yaz,'b--');
        horline(yel,'r--');
        h                               = plot([x1 x2],[yaz yaz],'b-');set(h,'LineWidth',3);
        h                               = plot([x1 x2],[yel yel],'r-');set(h,'LineWidth',3);
    end
elseif  strcmp(handles.stimulusshown,'b')
    % do nothing - blind detection
end
horline(0,'k-');
% YTick                                   = [floor(ax(3)/10)*10 (ax(3)+ax(4))/2 ceil(ax(4)/10)*10];
% if ax(3) == 0
%     YTick(1)                            = 0.2*ax(4);
% end
% YTick                                   = round(YTick);
% set(gca,'XTickLabel',[],'YTick',YTick);
set(gca,'XTickLabel',[]);

%% Graphics: Velocity Traces
axes(handles.axes_velocity);
cla;
% zoom on;
plot(TimeTrace,veltrace,'g'); hold on;
plot(TimeTrace,smvtrace,'b');
ax                                       = axis;
sel                                      = TimeTrace>det.start/1000 & TimeTrace < det.end/1000;
ax(4)                                    = 1.1*max(veltrace(sel));
ax(3)                                    = 1.1*min(veltrace(sel));
axis(ax);
ylabel('V');
for i                                    = 1:length(indxonoff);
    pa_verline(on(indxonoff(i))./handles.Fsample,'g');
    pa_verline(off(indxonoff(i))./handles.Fsample,'g');
end
% Minimum Detection Time Patch - Blue
xpatch                                   = [ax(1) ax(1) det.start/1000 det.start/1000];
ypatch                                   = [ax(3) ax(4) ax(4) ax(3)];
hp                                       = patch(xpatch,ypatch,'r'); set(hp,'FaceColor',[230 230 255]./255,'EdgeColor','b');
% Maximum Detection Time Patch - Blue
xpatch                                   = [det.end/1000 det.end/1000 ax(2) ax(2)];
ypatch                                   = [ax(3) ax(4) ax(4) ax(3)];
hp                                       = patch(xpatch,ypatch,'r'); set(hp,'FaceColor',[230 230 255]./255,'EdgeColor','b');
% Minimum Velocity Patch - Red
xpatch                                   = [ax(1) ax(1) ax(2) ax(2)]; ypatch = [ax(3) det.velocityon det.velocityon ax(3)];
hp                                       = patch(xpatch,ypatch,'r'); set(hp,'FaceColor',[255 230 230]./255,'EdgeColor','r');

% Plotting Traces Again
plot(TimeTrace,veltrace,'g');plot(TimeTrace,smvtrace,'b');
if ~isempty(indxsac)
    % Current Saccade Patch - Green
    xpatch                               = [curOnTime curOnTime curOffTime curOffTime];
    ypatch                               = [ax(3) ax(4) ax(4) ax(3)];
    hp                                   = patch(xpatch,ypatch,'g'); set(hp,'FaceColor',[230 255 230]./255,'EdgeColor','g');
    % Minimum Velocity Patch - Red
    xpatch                               = [ax(1) ax(1) ax(2) ax(2)];
    ypatch                               = [ax(3) det.velocityon det.velocityon ax(3)];
    hp                                   = patch(xpatch,ypatch,'r'); set(hp,'FaceColor',[255 230 230]./255,'EdgeColor','r');
    % Current Saccade Velocity and Smoothed Velocity
    plot(curTTrace,curveltrace,'g-');
    hl                                   = plot(curTTrace,cursmvtrace,'b');set(hl,'LineWidth',3);
    plot(TimeTrace(1:curOn),veltrace(1:curOn),'g');
    plot(TimeTrace(curOff:end),veltrace(curOff:end),'g');
    plot(TimeTrace(1:curOn),smvtrace(1:curOn),'b');
    plot(TimeTrace(curOff:end),smvtrace(curOff:end),'b');
    % Current Saccade Markers
    hl                                      = pa_verline(curOn/handles.Fsample,'g-'); set(hl,'LineWidth',2);
    hl                                      = pa_verline(curOff/handles.Fsample,'g-'); set(hl,'LineWidth',2);
end
% YTick                                       = [0.8*ax(3) (ax(3)+ax(4))/2 0.8*ax(4)];
% if ax(3)==0
%     YTick(1)                                = 0.2*ax(4);
% end
% YTick                                       = round(YTick);
set(gca,'XTickLabel',[]);
horline(handles.det.velocityoff,'m--');

%% Graphics: Acceleration Traces
axes(handles.axes_acceleration);
cla;
plot(TimeTrace,smatrace,'b');
hold on;
ax                                          = axis;
sel                                         = TimeTrace>det.start/1000 & TimeTrace < det.end/1000;
ax(4)                                       = 1.1*max(smatrace(sel));
ax(3)                                       = 1.1*min(smatrace(sel));
axis(ax)
for i                                       = 1:length(indxonoff);
    pa_verline(on(indxonoff(i))./handles.Fsample,'g');
    pa_verline(off(indxonoff(i))./handles.Fsample,'g');
end
% Minimum Detection Time Patch - Blue
xpatch                                   = [ax(1) ax(1) det.start/1000 det.start/1000];
ypatch                                   = [ax(3) ax(4) ax(4) ax(3)];
hp                                       = patch(xpatch,ypatch,'r'); set(hp,'FaceColor',[230 230 255]./255,'EdgeColor','b');
% Maximum Detection Time Patch - Blue
xpatch                                   = [det.end/1000 det.end/1000 ax(2) ax(2)];
ypatch                                   = [ax(3) ax(4) ax(4) ax(3)];
hp                                       = patch(xpatch,ypatch,'r'); set(hp,'FaceColor',[230 230 255]./255,'EdgeColor','b');
% Plotting Traces Again
plot(TimeTrace,smatrace,'b');
if ~isempty(indxsac)
    % Current Saccade Patch - Green
    xpatch                               = [curOnTime curOnTime curOffTime curOffTime];
    ypatch                               = [ax(3) ax(4) ax(4) ax(3)];
    hp                                   = patch(xpatch,ypatch,'g'); set(hp,'FaceColor',[230 255 230]./255,'EdgeColor','g');
    % Current Saccade Velocity and Smoothed Velocity
    hl                                  = plot(curTTrace,cursmatrace,'b-'); set(hl,'LineWidth',3);
    % Current Saccade Markers
    hl                                  = pa_verline(curOnTime,'g-'); set(hl,'LineWidth',2);
    hl                                  = pa_verline(curOffTime,'g-'); set(hl,'LineWidth',2);
end
% Labeling
xlabel('Time (ms)');
ylabel('A (x10^4)');
horline(0,'r');
YTick                                       = [0.8*ax(3) (0.8*ax(3)+0.8*ax(4))/2 0.8*ax(4)];
if ax(3)==0
    YTick(1)                                = 0.2*ax(4);
end
YTick                                       = round(YTick);
set(gca,'YTick',YTick,'YTickLabel',round(YTick/1000));



function handles        = MAIN_redetect(handles)
switch handles.redetect
    case 'sub'
        OnOffExist              = handles.onoff;
        sel                     = OnOffExist(3,:)<handles.ctrial;
        OnOffExist              = OnOffExist(:,sel);
        handles                 = MAIN_detsaccade(handles);
        sel                     = handles.onoff(3,:)>=handles.ctrial;
        OnOff                   = handles.onoff(:,sel);
        OnOff                   = [OnOffExist OnOff];
        handles.onoff           = OnOff;
    case 'all'
        handles                 = MAIN_detsaccade(handles);
    case 'cur'
        OnOffExist              = handles.onoff;
        selpre                  = OnOffExist(3,:)<handles.ctrial;
        selpost                 = OnOffExist(3,:)>handles.ctrial;
        OnOffExistPre           = OnOffExist(:,selpre);
        OnOffExistPost          = OnOffExist(:,selpost);
        handles                 = MAIN_detsaccade(handles);
        sel                     = handles.onoff(3,:)==handles.ctrial;
        OnOff                   = handles.onoff(:,sel);
        if ~isempty(OnOff)
            OnOff               = [OnOffExistPre OnOff OnOffExistPost];
            handles.onoff       = OnOff;
        else
            OnOff               = [OnOffExistPre OnOffExistPost];
            handles.onoff       = OnOff;
        end
end
handles.mfigcheck           = false;
handles                     = MAIN_plottrial(handles);


% --- Executes on button press in btn_movestart.
function btn_movestart_Callback(hObject, eventdata, handles)
% hObject    handle to btn_movestart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

csaccade                                    = handles.csaccade;
ctrial                                      = handles.ctrial;
indx                                        = find(handles.onoff(3,:)== ctrial & handles.onoff(4,:) == csaccade);
axes(handles.axes_trace);
onnew                                       = ginput(1).*handles.Fsample;
onnew                                       = round(onnew(1));
if onnew<handles.onoff(2,indx)
    handles.onoff(1,indx)                   = onnew;
    handles                                 = MAIN_plottrial(handles);
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in btn_moveend.
function btn_moveend_Callback(hObject, eventdata, handles)
% hObject    handle to btn_moveend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

csaccade                                    = handles.csaccade;
ctrial                                      = handles.ctrial;
indx                                        = find(handles.onoff(3,:)==ctrial & handles.onoff(4,:)==csaccade);
axes(handles.axes_trace);
offnew                                      = ginput(1).*handles.Fsample;
offnew                                      = round(offnew(1));
if offnew>handles.onoff(1,indx)
    handles.onoff(2,indx)                   = offnew;
    handles                                 = MAIN_plottrial(handles);
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in btn_delsac.
function btn_delsac_Callback(hObject, eventdata, handles)
% hObject    handle to btn_delsac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Detect Current Saccade in Current Trial
csaccade                                    = handles.csaccade;
ctrial                                      = handles.ctrial;
trialnr                                     = handles.onoff(3,:);
sacnr                                       = handles.onoff(4,:);
indx                                        = find(trialnr==ctrial&sacnr==csaccade);
if isempty(indx)
    str1                                    = get(handles.display,'String');
    str2                                    = 'No saccade present that can be deleted!';
    str                                     = char({str1;str2});
    set(handles.display,'String',str);
    return
elseif length(indx)==1
    sel                                     = trialnr==ctrial;
    maxsacnr                                = max(handles.onoff(4,sel));
    % And Remove Current Saccade from OnOff-matrix
    handles.onoff(:,indx)                   = [];
    % Then Go To Next Saccade in Trial OR Remain in trial when no Saccade
    % present
    trialnr                                 = handles.onoff(3,:);
    %     sacnr                                   = handles.onoff(4,:);
    indx                                    = find(trialnr==ctrial);
    if ~isempty(indx)
        handles.onoff(4,indx)               = 1:length(indx);
        if csaccade == maxsacnr;
            csaccade                        = csaccade-1;
        elseif csaccade <  maxsacnr
            % do nothing
        end
        handles.csaccade                    = csaccade;
    end
    handles = MAIN_plottrial(handles);
elseif length(indx)>1
    str1                                    = get(handles.display,'String');
    str2                                    = 'Uh-Ooooh! BugWare!';
    str                                     = char({str1;str2});
    set(handles.display,'String',str);
    return
end
handles.mfigcheck                           = false;
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in btn_inssac.
function btn_inssac_Callback(hObject, eventdata, handles)
% hObject    handle to btn_inssac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

onoff                           = handles.onoff;
x                               = ginput(1);
onnew                           = x(1).*handles.Fsample;
onnew                           = round(onnew);
axes(handles.axes_trace);
pa_verline(onnew,'g');

x                               = ginput(1);
offnew                          = x(1).*handles.Fsample;
offnew                          = round(offnew);
ctrial                          = handles.ctrial;
indx                            = find(handles.onoff(3,:)==ctrial);
if isempty(indx)
    csaccade                    = 1;
    onoff_insert                = [onnew;offnew;ctrial;csaccade;1];
    indxprev                    = find(handles.onoff(3,:)<ctrial);
    if ctrial == 1
        onoff                   = [onoff_insert onoff];
	else
        onoff                   = [onoff(:,1:max(indxprev)) onoff_insert onoff(:,max(indxprev)+1:end)];
    end
elseif ~isempty(indx)
    csaccade                    = length(indx)+1;
    onoff_insert                = [onnew;offnew;ctrial;csaccade;1];
    onoff                       = [onoff(:,1:max(indx)) onoff_insert onoff(:,max(indx)+1:end)];
end
handles.onoff                   = onoff;
handles.csaccade                = csaccade;
handles                         = MAIN_plottrial(handles);

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_save_Callback(hObject, eventdata, handles) %#ok<INUSL>
% hObject    handle to menu_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the Tag of the menu selected
Tag         = get(hObject, 'Tag');
OnOff       = handles.onoff;
% Only menu_save detected trials
sel         = OnOff(5,:) == 1;
% ctrial      = handles.ctrial;
% sel         = OnOff(3,:)<=ctrial;
OnOff       = OnOff(:,sel);
% Produce OnOff-markers to insert in Sac-file
OnOff       = [OnOff(1,:)+(OnOff(3,:)-1)*handles.Nsamples OnOff(2,:)+(OnOff(3,:)-1)*handles.Nsamples];

% Based on the item selected, take the appropriate action
switch Tag
    case 'menu_save'
        File = handles.fname;
        File = pa_fcheckext(File,'.sac');
        str1 = get(handles.display,'String');
        str2 = ['Save ' upper(File)];
        str = char({str1;str2});
        set(handles.display,'String',str);
        fid = fopen(File,'w','l');
        fwrite(fid,OnOff,'ulong');
        fclose(fid);
    case 'menu_save_as'
        % Allow the user to select the file name to menu_save to
        [filename, pathname] = uiputfile( ...
            {'*.sac';'*.*'}, ...
            'Save as');
        % If 'Cancel' was selected then return
        if isequal([filename,pathname],[0,0])
            return
        else
            % Construct the full path and menu_save
            File    = filename;
            fid     = fopen(File,'w','l');
            fwrite(fid,OnOff,'ulong');
            fclose(fid);
        end
    case 'menu_exit'
        menu_exitsave(handles,OnOff);
    case 'ultradet'
        h           = get(handles.ultradet);
        CC          = h.CurrentCharacter;
        if strcmpi(CC,'s');
            File    = handles.fname;
            File    = [File(1:end-2) 'sac'];
            str1    = get(handles.display,'String');
            str2    = ['Save ' upper(File)];
            str     = char({str1;str2});
            set(handles.display,'String',str);
            fid     = fopen(File,'w','l');
            fwrite(fid,OnOff,'ulong');
            fclose(fid);
        elseif  strcmpi(CC,'x') || strcmpi(CC,'n') ;
            menu_exitsave(handles,OnOff);
        end
    case 'btn_nexttrial'
        menu_exitsave(handles,OnOff);
    case 'btn_nextsac'
        menu_exitsave(handles,OnOff);
    case 'load_new'
        Fname = handles.fname;
        Fname = [Fname(1:end-2) 'sac'];
        if exist(Fname,'file')
            txt = ['Saving SAC-file ' Fname ', overwrite existing one?'];
            name = 'Loading New HV-File, Saving Existing SAC-file';
            numl = 1;
            da = {'YES / no / new sacfile-name'};
            answer=inputdlg(txt,name,numl,da);
            if ~isempty(answer)
                if strcmpi(answer(1),'y') || strcmpi(answer(1:end),'yes');
                    File = handles.fname;
                    File = [File(1:end-2) 'sac'];
                    delete(File);
                    fid = fopen(File,'w','l');
                    fwrite(fid,OnOff,'ulong');
                    fclose(fid);
                    str1 = get(handles.display,'String');
                    str2 = ['   Save ' upper(File)];
                    str = char({str1;str2});
                    set(handles.display,'String',str);
                elseif strcmpi(answer,'n') || strcmpi(answer,'no')
                    return
                elseif strcmpi(answer(end-2:end),'sac')
                    File = answer;
                    warndlg(['Pressing OK will overwrite ' File '!']);
                    if exist(File,'file')
                        delete(File);
                    end
                    fid = fopen(File,'w','l');
                    fwrite(fid,OnOff,'ulong');
                    fclose(fid);
                end
            end
        else
            File = Fname;
            str1 = get(handles.display,'String');
            str2 = ['   Save ' upper(File)];
            str = char({str1;str2});
            set(handles.display,'String',str);
            fid = fopen(File,'w','l');
            fwrite(fid,OnOff,'ulong');
            fclose(fid);
        end
end

function menu_exitsave(handles,OnOff)
SacFile             = pa_fcheckext(handles.fname,'sac');
if exist(SacFile,'file')
    txt             = ['Overwrite existing sac-file, ' SacFile '?'];
    name            = 'Exit Ultradet';
    answer          = questdlg(txt,name,'Yes','Cancel','New Sacfile','Yes');
    if ~isempty(answer)
        if strcmpi(answer,'yes');
            File    = pa_fcheckext(handles.fname,'sac');
            delete(File);
            fid     = fopen(File,'w','l');
            fwrite(fid,OnOff,'ulong');
            fclose(fid);
            str1    = get(handles.display,'String');
            str2    = ['   Save ' upper(File)];
            str     = char({str1;str2});
            set(handles.display,'String',str);
        elseif strcmpi(answer,'cancel') ;
            return
        else
            txt             = 'Choose a new name';
            name            = 'Exit Ultradet';
            numl            = 1;
            da              = {SacFile};
            answer          = inputdlg(txt,name,numl,da);
            File            = answer{1};
            if exist(File,'file');
                delete(File);
            end
            fid     = fopen(File,'w','l');
            fwrite(fid,OnOff,'ulong');
            fclose(fid);
        end
        close(handles.ultradet);
    end
else
    File            = SacFile;
    str1            = get(handles.display,'String');
    str2            = ['   Save ' upper(File)];
    str             = char({str1;str2});
    set(handles.display,'String',str);
    fid             = fopen(File,'w','l');
    fwrite(fid,OnOff,'ulong');
    fclose(fid);
    close(handles.ultradet);
end


% --------------------------------------------------------------------
function menu_exit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

menu_save_Callback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function menu_load_new_Callback(hObject, eventdata, handles)
% hObject    handle to menu_load_new (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

menu_save_Callback(hObject, eventdata, handles);
[filename, pathname] = uigetfile( ...
    {'*.hv', 'All HV-Files (*.hv)'; ...
    '*.*','All Files (*.*)'}, ...
    'Select Trace File');
% If "Cancel" is selected then return
if isequal([filename,pathname],[0,0])
    return
    % Otherwise construct the fullfilename and Check and load the file
else
    str1 = get(handles.display,'String');
    str2 = ['   Exiting ' upper(handles.fname)];
    str = char({str1;str2});
    set(handles.display,'String',str);
    str1 = get(handles.display,'String');
    str2 = ['CD to ' pathname 'to load ' upper(filename)];
    str = char({str1;str2});
    set(handles.display,'String',str);
    cd(pathname);
    handles=MAIN_Check_And_Load(handles,filename);
    handles = MAIN_plottrial(handles);
end
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in btn_nexttrial.
function btn_nexttrial_Callback(hObject, eventdata, handles)
% hObject    handle to btn_nexttrial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


ctrial                          = handles.ctrial;
ntrial                          = handles.Ntrial;
if ctrial < ntrial;
    ctrial                      = ctrial+1; % Go to next trial
    handles.ctrial              = ctrial; % Make it current
    indx                        = find(handles.onoff(3,:)==ctrial); % Find Saccades
    if ~isempty(indx);
        handles.onoff(5,indx)   = ones(size(indx)); % Check Saccades as Detected
        csaccade                = handles.onoff(4,indx);
        handles.csaccade        = csaccade(1); % Make first saccade current
    end
    handles = MAIN_plottrial(handles); % Plot Again
    % Update handles structure
    guidata(hObject, handles);
elseif ctrial == ntrial;
    rndsong     = floor(2*rand(1));
    Fs          = 50000;
    if rndsong
        load handel;
    else
        fname   = 'VoodooChild.wav';
        y       = wavread(fname);
    end
    sound(y,Fs)
    menu_save_Callback(hObject, eventdata, handles);
    return
end


% --- Executes on button press in btn_prevtrial.
function btn_prevtrial_Callback(hObject, eventdata, handles)
% hObject    handle to btn_prevtrial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ctrial = handles.ctrial;
if ctrial > 1
    ctrial = ctrial-1;
    handles.ctrial = ctrial;
    indx = find(handles.onoff(3,:)==ctrial);
    if ~isempty(indx);
        csaccade = handles.onoff(4,indx);
        handles.csaccade = csaccade(1);
    end
    handles = MAIN_plottrial(handles);
end

% Update handles structure
guidata(hObject, handles);


%% --- Executes on button press in btn_nextsac.
function btn_nextsac_Callback(hObject, eventdata, handles)
% hObject    handle to btn_nextsac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ctrial                              = handles.ctrial;
csaccade                            = handles.csaccade;
trialnr                             = handles.onoff(3,:);
sacnr                               = handles.onoff(4,:);
indx                                = find(trialnr==ctrial);
nsaccade                            = max(sacnr(indx));

if ~(ctrial == max(handles.Ntrial) && csaccade == nsaccade)     % unless current saccade is last saccade of experimental run
    if csaccade < max(nsaccade);                                % if current saccade is not last saccade in trial
        csaccade                    = csaccade+1;               % go to next saccade
        handles.csaccade            = csaccade;                 % and make it current
        handles.onoff(5,indx)       = 1;
    elseif  isempty(nsaccade) ||csaccade == max(nsaccade);      % however, if current saccade is the last saccade of the trial
        handles.ctrial              = ctrial+1;                 % go to next trial and make it current
        handles.csaccade            = 1;                        % and make first saccade in that trial current
        ctrial                      = handles.ctrial;
        sel                        = trialnr==ctrial;    % Find which saccade-indixes belong to current trial
%         nsaccade                    = sacnr(indx);              % and determine the number of saccades in current trial
        handles.onoff(5,sel)       = 1;
%         if handles.ctrial==handles.Ntrial;
%             nsaccade                = 1;
%         end
%         while isempty(nsaccade);                                % and if there are no saccades in current trial
%             handles.ctrial          = ctrial+1;                 % go to next trial and make it current
%             ctrial                  = handles.ctrial;
%             handles.csaccade        = 1;                        % and make first saccade current
%             indx                    = find(trialnr==ctrial);    % Find Current Trial in Trials containing saccades
%             nsaccade                = sacnr(indx);              % Number of saccades in Current Trial
%             handles.onoff(5,indx)   = 1;
%             if handles.ctrial==handles.Ntrial;
%                 nsaccade            = 1;
%             end
%         end
    end
elseif (ctrial == max(handles.Ntrial) && csaccade == nsaccade)
    load handel;
    sound(y,Fs);
    menu_save_Callback(hObject, eventdata, handles);
    return
end
handles = MAIN_plottrial(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in btn_prevsac.
function btn_prevsac_Callback(hObject, eventdata, handles)
% hObject    handle to btn_prevsac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
csaccade = handles.csaccade;
ctrial = handles.ctrial;
trialnr = handles.onoff(3,:);
sacnr = handles.onoff(4,:);
% indx = find(trialnr==ctrial);
% nsaccade = max(sacnr(indx));
if ~(ctrial == 1 && csaccade == 1);
    if csaccade > 1;
        csaccade = csaccade-1;
        handles.csaccade = csaccade;
    elseif csaccade == 1;
        handles.ctrial = ctrial - 1;
        indx = find(trialnr==handles.ctrial);
        while isempty(indx);            % and if there are no saccades in current trial
            handles.ctrial = ctrial - 1;        % go to next trial and make it current
            ctrial = handles.ctrial;
            indx = find(trialnr==ctrial);
            if handles.ctrial==1;
                indx = 1;
            end
        end
        nsaccade = sacnr(indx);
        handles.csaccade = max(nsaccade);
    end
    handles = MAIN_plottrial(handles);
end

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function Tools_Callback(hObject, eventdata, handles)
% hObject    handle to menu_tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_about_Callback(hObject, eventdata, handles)
% hObject    handle to menu_about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
message = {'PA_SACDET';'version 1.0';'Matlab version: R2009B';'Date: September 2011';'Author: Marc van Wanrooij';...
    'Co-Authors: John van Opstal / Jeroen Goossens / Paul Hofman';...
    'Short Info: Saving occurs up until CURRENT trial, Re-detecting starts for current trial'};
title= 'ULTRADET';
msgbox(message,title)
% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_help_Callback(hObject, eventdata, handles)
% hObject    handle to menu_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Update handles structure
guidata(hObject, handles);


function ultradet_KeyPressFcn(hObject, eventdata, handles)
% --- Executes on key press over pa_sacdet with no controls selected.
% P     = Previous Trial
% N     = Next Trial
% Space = Next Trial
% S     = menu_save
h                                      = get(handles.ultradet);
CC                                     = h.CurrentCharacter;
if strcmpi(CC,'P')
    btn_prevtrial_Callback(hObject, eventdata, handles);
elseif strcmpi(CC,'N')
    btn_nextsac_Callback(hObject, eventdata, handles);
elseif strcmpi(CC,'D')
    btn_delsac_Callback(hObject, eventdata, handles);
elseif strcmpi(CC,'I')
    btn_inssac_Callback(hObject, eventdata, handles);
elseif strcmpi(CC,'X')
    menu_exit_Callback(hObject, eventdata, handles);
elseif strcmpi(CC,'L')
    menu_load_new_Callback(hObject, eventdata, handles);
elseif strcmpi(CC,'B')
    btn_movestart_Callback(hObject, eventdata, handles);
elseif strcmpi(CC,'E') || strcmpi(CC,'M')
    btn_moveend_Callback(hObject, eventdata, handles);
elseif strcmpi(CC,'R')
    handles = MAIN_redetect(handles);
    % Update handles structure
    guidata(hObject, handles);
elseif strcmpi(CC,' ')
    btn_nexttrial_Callback(hObject, eventdata, handles);
elseif strcmpi(CC,'S')
    menu_save_Callback(hObject, eventdata, handles);
end


function txt_vel_Callback(hObject, eventdata, handles)
% Online Velocity Parameter Adjustment
% Get User Input
user_entry                             = str2double(get(hObject,'string'));
handles.det.velocityon                   = user_entry;
% MAIN_redetect Saccades
handles                                = MAIN_redetect(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function txt_vel_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% Update handles structure
guidata(hObject, handles);

function txt_amp_Callback(hObject, eventdata, handles)
% Online Amplitude Parameter Adjustment
user_entry                             = str2double(get(hObject,'string'));
handles.det.amplitude                  = user_entry;
% MAIN_redetect Saccades
handles                                = MAIN_redetect(handles);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_amp_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% Update handles structure
guidata(hObject, handles);

function txt_duration_Callback(hObject, eventdata, handles)
% Online Duration Parameter Adjustment
user_entry                             = str2double(get(hObject,'string'));
handles.det.duration                   = user_entry;
% MAIN_redetect Saccades
handles                                = MAIN_redetect(handles);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_duration_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% Update handles structure
guidata(hObject, handles);


function txt_detstart_Callback(hObject, eventdata, handles)
% Online Start of Detection Parameter Adjustment
user_entry                              = str2double(get(hObject,'string'));
handles.det.start                       = user_entry;
% MAIN_redetect Saccades
handles                                 = MAIN_redetect(handles);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_detstart_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% Update handles structure
guidata(hObject, handles);


function txt_jumpsac_Callback(hObject, eventdata, handles)
% Online Saccade Number Adjustment
user_entry                              = str2double(get(hObject,'string'));
handles.csaccade                        = user_entry;
handles = MAIN_plottrial(handles);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_jumpsac_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% Update handles structure
guidata(hObject, handles);


function txt_jumptrial_Callback(hObject, eventdata, handles)
% Online Trial Number Adjustment
user_entry                              = str2double(get(hObject,'string'));
if user_entry <= handles.Ntrial
    handles.ctrial                      = user_entry;
    ctrialindx                          = find(handles.onoff(3,:)==handles.ctrial);
    if ~isempty(ctrialindx);
        ctrialindx                      = min(ctrialindx);
        handles.csaccade                = handles.onoff(4,ctrialindx);
        handles.onoff(5,ctrialindx)     = 1;
        set(handles.txt_jumpsac,'String',num2str(handles.csaccade));
        handles = MAIN_plottrial(handles);
    elseif isempty(ctrialindx);
        handles = MAIN_plottrial(handles);
    end
else
    str1                                = get(handles.display,'String');
    str2                                = '   Warning !';
    str                                 = char({str1;str2});
    set(handles.display,'String',str);
    if (user_entry-handles.Ntrial)>1
        str1                            = get(handles.display,'String');
        str2                            = ['    User Entry exceeds number of trials in experiment/hv-file by ' num2str(user_entry-handles.Ntrial) ' trials!'];
        str                             = char({str1;str2});
        set(handles.display,'String',str);
    elseif (user_entry-handles.Ntrial)==1
        str1                            = get(handles.display,'String');
        str2                            = ['    User Entry exceeds number of trials in experiment/hv-file by ' num2str(user_entry-handles.Ntrial) ' trial!'];
        str                             = char({str1;str2});
        set(handles.display,'String',str);
    end
    str1                                = get(handles.display,'String');
    str2                                = '    Please use correct trial number.';
    str                                 = char({str1;str2});
    set(handles.display,'String',str);
    set(handles.txt_jumptrial,'String',num2str(handles.ctrial));  % And set the JumpToSaccadeText congruent with this current saccade
end
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function txt_jumptrial_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% Update handles structure
guidata(hObject, handles);

function txt_detend_Callback(hObject, eventdata, handles)
% Online End of Detection Parameter Adjustment
user_entry                              = str2double(get(hObject,'string'));
handles.det.end                         = user_entry;
% MAIN_redetect Saccades
handles                                 = MAIN_redetect(handles);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_detend_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% Update handles structure
guidata(hObject, handles);


function txt_smooth_Callback(hObject, eventdata, handles)
% Online Smoothing Parameter Adjustment
user_entry                              = str2double(get(hObject,'string'));
handles.det.smooth                      = user_entry;
[veltrace smvtrace]                     = getvel(handles.htrace,handles.vtrace,handles.det.smooth,1./handles.Fsample);
handles.veltrace                        = veltrace;    % deg/s
handles.smvtrace                        = smvtrace;  % deg/s
% MAIN_redetect Saccades
handles                                 = MAIN_redetect(handles);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_smooth_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% Update handles structure
guidata(hObject, handles);


function l                              = logtype
%
% FUNCTION L = LOGTYPE
%
%   Structure holding log information
%
% .. Paul ..
%
l                                       = struct ('Name', '', ...
    'NTest',   [], 'NPresent', [], 'NTarget',  [], ...
    'NTrial',  [], ...
    'FSample', [], 'NSample',  [], 'NChan', [], ...
    'ChanStr', [], 'ChanX',  [], 'ChanY',   []);


function param                          = fgetpar(fp, par, format)
% FUNCTION param = fgetpar(fp, par, format)
%
%   read parameter(s) "par" form text file fp
%   the search procedure for parameter names
%   is case sensitive.
%
%   Jeroen Goossens
% search from beginning
frewind(fp);
par                                     = ['#' par];
while ~feof(fp)
    % search par identifier
    ident                               = fscanf(fp,'%s',1);
    if strcmp(ident,par)==1
        % find '=' sign
        symb                            = fscanf(fp,'%s',1);
        if strcmp(symb,'=')==1
            % read value(s)
            param                       = fscanf(fp,format,1);
            return
        end
        return
    end
end
if ~exist('param', 'var'), param        = []; end

function c = isin (a, b)
% FUNCTION C = ISIN(A,B)
%
% Are the elements of A to be found in B?
% Answer in C: 0 = no
%              1 = yes
%
% Marcus
c = zeros(size(a));
for i=1:length(b)
    c = c | (a == b(i));
end;

function hhor                             = horline(ypos,marker)
% Plot horizontal line
if nargin<2
    marker                                = 'k--';
    if nargin<1
        ypos                              = 0;
    end
end
ax                                        = axis;
hhor                                      = plot([ax(1) ax(2)],[ypos ypos],marker);


function [c,b]                            = my_unique(a)
% Get amount of occurrences (c) for each unique value (b) in matrix a.
a                                         = sort(a);
b                                           = unique(a);
c                                         = [];
for i                                     = 1:length(b)
    indx                                  = find(a==b(i));
    length(indx);
    tmp                                   = (1:length(indx))';
    c                                     = [c;tmp]; %#ok<AGROW>
end

% --------------------------------------------------------------------
function menu_q_key_menu_Callback(hObject, eventdata, handles)
% hObject    handle to get_start_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpstr={'N - Next Saccade/Trial';'P - Previous Saccade/Trial';...
    'D - Delete Saccade';'I - Insert Saccade';'R - Re-detect Current and Subsequent Trials';...
    'B - Move Saccade Beginning';'E - Move Saccade Ending';...
    'S - Save/Overwrite Sac-file';'L - Load New File';'X - Exit'};
dlgname='Shortcuts / Quick Keys';
helpdlg(helpstr,dlgname);

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_axdef_Callback(hObject, eventdata, handles)
% hObject    handle to menu_axdef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(gcbo, 'Checked'),'off')
    set(gcbo, 'Checked', 'on');
    set(handles.menu_axmanual,'Checked','off');
    set(handles.menu_axauto,'Checked','off');
    handles.axset = 'd';
    handles = MAIN_plottrial(handles);
end
% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_axauto_Callback(hObject, eventdata, handles)
% hObject    handle to menu_axauto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(gcbo, 'Checked'),'off')
    set(gcbo, 'Checked', 'on');
    set(handles.menu_axdef,'Checked','off')
    set(handles.menu_axmanual,'Checked','off')
    handles.axset = 'a';
    handles = MAIN_plottrial(handles);
end
% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_axmanual_Callback(hObject, eventdata, handles)
% hObject    handle to menu_axmanual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(gcbo, 'Checked'),'off')
    set(gcbo, 'Checked', 'on');
    set(handles.menu_axdef,'Checked','off')
    set(handles.menu_axauto,'Checked','off')
    handles.axset = 'm';
    handles = MAIN_plottrial(handles);
end
% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_axis_set_Callback(hObject, eventdata, handles)
% hObject    handle to menu_axis_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_tools_Callback(hObject, eventdata, handles)
% hObject    handle to menu_tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Update handles structure
guidata(hObject, handles);


function minax_txt_Callback(hObject, eventdata, handles)
% hObject    handle to minax_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minax_txt as text
%        str2double(get(hObject,'String')) returns contents of minax_txt as a double
handles = MAIN_plottrial(handles);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function minax_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minax_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% Update handles structure
guidata(hObject, handles);

function maxax_txt_Callback(hObject, eventdata, handles)
% hObject    handle to maxax_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxax_txt as text
%        str2double(get(hObject,'String')) returns contents of maxax_txt as a double
handles = MAIN_plottrial(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function maxax_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxax_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% Update handles structure
guidata(hObject, handles);



% --- Executes on selection change in display.
function display_Callback(hObject, eventdata, handles)
% hObject    handle to display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns display contents as cell array
%        contents{get(hObject,'Value')} returns selected item from display
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_delsacfile_Callback(hObject, eventdata, handles)
% hObject    handle to menu_delsacfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sacfile     = pa_fcheckext(handles.fname,'.sac');
if exist(sacfile,'file')
    delete(sacfile);
    str1        = get(handles.display,'String');
    str2        = ['   Deleting ' upper(sacfile)];
    str         = char({str1;str2});
    set(handles.display,'String',str);
else
    str1        = get(handles.display,'String');
    str2        = [upper(sacfile) ' does not exist, and cannot be deleted'];
    str         = char({str1;str2});
    set(handles.display,'String',str);
end

% --------------------------------------------------------------------
function menu_optional_axes_Callback(hObject, eventdata, handles)
% hObject    handle to menu_optional_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_opt_radial_plot_Callback(hObject, eventdata, handles)
% hObject    handle to opt_polar_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(gcbo, 'Checked'),'off')
    set(gcbo, 'Checked', 'on');
    set(handles.menu_opt_mainseq_axes,'Checked','off');
    handles.extraaxes = 'r';
    handles = MAIN_plottrial(handles);
end
% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_opt_mainseq_axes_Callback(hObject, eventdata, handles)
% hObject    handle to menu_opt_mainseq_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(gcbo, 'Checked'),'off')
    set(gcbo, 'Checked', 'on');
    set(handles.menu_opt_radial_plot,'Checked','off');
    handles.extraaxes = 'm';
    handles = MAIN_plottrial(handles);
end
% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_blind_Callback(hObject, eventdata, handles)
% hObject    handle to menu_blind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(gcbo, 'Checked'),'off')
    set(gcbo, 'Checked', 'on');
    set(handles.menu_visible,'Checked','off');
    handles.stimulusshown = 'b';
    handles = MAIN_plottrial(handles);
end
% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_visible_Callback(hObject, eventdata, handles)
% hObject    handle to menu_visible (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(gcbo, 'Checked'),'off')
    set(gcbo, 'Checked', 'on');
    set(handles.menu_blind,'Checked','off');
    handles.stimulusshown = 'v';
    handles = MAIN_plottrial(handles);
end
% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_stimulus_Callback(hObject, eventdata, handles)
% hObject    handle to menu_stimulus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Update handles structure
guidata(hObject, handles);



% --------------------------------------------------------------------
function menu_load_sac_Callback(hObject, eventdata, handles)
% hObject    handle to menu_load_sac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[sacfile, pathname] = uigetfile( ...
    {'*.sac', 'All Sac-Files (*.sac)'; ...
    '*.*','All Files (*.*)'}, ...
    'Select SAC File');
% If "Cancel" is selected then return
if isequal([sacfile,pathname],[0,0])
    return
    % Otherwise construct the fullfilename and Check and load the file
else
    str1 = get(handles.display,'String');
    str2 = ['CD to ' pathname 'to load ' upper(sacfile)];
    str = char({str1;str2});
    set(handles.display,'String',str);
    cd(pathname);
    if exist(sacfile,'file')
        OnOffExist                              = loadsac(sacfile);
        str1                                    = get(handles.display,'String');
        str2                                    = ['Loading ' upper(sacfile)];
        str                                     = char({str1;str2});
        set(handles.display,'String',str);
        OnExist                                 = OnOffExist(1,:);
        TrialExist                              = floor(OnExist./handles.Nsamples)+1;
        OnExist                                 = OnExist - ((TrialExist-1).*handles.Nsamples);
        OffExist                                = OnOffExist(2,:);
        OffExist                                = OffExist - ((TrialExist-1).*handles.Nsamples);
        SacNrExist                              = my_unique(TrialExist)';
        CheckExist                              = ones(size(OnExist));
        OnOffExist                              = [OnExist; OffExist; TrialExist; SacNrExist; CheckExist];
        trialn                                  = max(TrialExist);
        indx                                    = find(handles.onoff(3,:)==trialn); %#ok<MXFND>
        indx                                    = max(indx);
        handles.onoff                           = [OnOffExist handles.onoff(:,indx+1:end)];
        handles.ctrial                          = trialn;
        handles.csaccade                        = OnOffExist(4,end);
    end
end

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function panel_redetect_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to panel_redetect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch get(hObject,'Tag')   % Get Tag of selected object
    case 'radiobtn_currentandsub'
        handles.redetect = 'sub';
    case 'radiobtn_all'
        handles.redetect = 'all';
    case 'radiobtn_current'
        handles.redetect = 'cur';
end

% Update handles structure
guidata(hObject, handles);



function txt_veloff_Callback(hObject, eventdata, handles)
% hObject    handle to txt_veloff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_veloff as text
%        str2double(get(hObject,'String')) returns contents of txt_veloff as a double
% Online Velocity Parameter Adjustment
% Get User Input
user_entry                             = str2double(get(hObject,'string'));
handles.det.velocityoff                = user_entry;
% MAIN_redetect Saccades
handles                                = MAIN_redetect(handles);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_veloff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_veloff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in chck_acc.
function chck_acc_Callback(hObject, eventdata, handles)
% hObject    handle to chck_acc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chck_acc

user_entry                             = get(hObject,'Value');
handles.det.acc                        = user_entry;
% MAIN_redetect Saccades
handles                                = MAIN_redetect(handles);

% Update handles structure
guidata(hObject, handles);

function [veltrace,smvtrace,acctrace,smatrace]          = getvel(htrace,vtrace,Fsample,sd)
% Obtain radial velocity from horizontal HOR and vertical VER traces.
%
% [VEL, SVEL] = GETVEL(HOR,VER,SMOOTHFACTOR,FSAMPLE)
%
% Obtain radial velocity from horizontal HOR and vertical VER traces.
%
% See also GSMOOTH, PA_SACDET
%
% MarcW 2007

Rx                                      = htrace;
Ry                                      = vtrace;
R                                       = NaN*Rx;
veltrace                                = R;
acctrace                                = R;
smvtrace                                = R;
smatrace                                = R;
for i                                   = 1:size(htrace,2)
%     Rx(:,i)                             = Rx(:,i)-Rx(1,i);
%     Ry(:,i)                             = Ry(:,i)-Ry(1,i);
    Rx(:,i)                             = gradient(Rx(:,i),1);
    Ry(:,i)                             = gradient(Ry(:,i),1);
    R(:,i)                              = hypot(Rx(:,i),Ry(:,i));
    R(:,i)                              = cumsum(R(:,i));
    


    veltrace(:,i)                       = gradient(R(:,i),1./Fsample);
    smvtrace(:,i)                       = pa_gsmooth(veltrace(:,i),Fsample,sd);
    acctrace(:,i)                       = gradient(smvtrace(:,i),1./Fsample);
    smatrace(:,i)                       = pa_gsmooth(acctrace(:,i),Fsample,sd);
end

