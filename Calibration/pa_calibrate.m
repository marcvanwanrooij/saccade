function varargout = pa_calibrate(varargin)
% Check fixation and train calibration networks
%
% PA_CALIBRATE(FNAME)
%
%  PA_CALIBRATE
%
%   Train backpropagation networks with the calibration data in file FNAME.
%
%   NOTE: a good calibrated network will produce a coherent calibration rose,
%   and has little errors with mu~0 and std<2.0.
%
%
%  See also CALIBRATE, ULTRADET
%
%  Author: Marcus
%  Date: 06-04-07


% Edit the above text to modify the response to help fixcheck

% Last Modified by GUIDE v2.5 17-Jan-2012 08:05:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
	'gui_Singleton',  gui_Singleton, ...
	'gui_OpeningFcn', @pa_calibrate_OpeningFcn, ...
	'gui_OutputFcn',  @pa_calibrate_OutputFcn, ...
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


% --- Executes just before pa_calibrate is made visible.
function pa_calibrate_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pa_calibrate (see VARARGIN)

% Choose default command line output for pa_calibrate
handles.output                      = hObject;

handles                             = Check_And_Load(handles, varargin{:});
handles                             = trainnet(handles);
home;
disp('NEW CALL');

handles                             = plotXYZvsAZEL(handles);
handles                             = plottrain(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pa_calibrate wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%% CHECK and LOAD
function handles                    = Check_And_Load(handles,varargin)
% remove some warnings
warning('off','MATLAB:plot:IgnoreImaginaryXYPart');

lv                                  = length(varargin);
if lv<1
	handles.fname                   = [];
	handles.fname                   = pa_fcheckexist(handles.fname,'*.dat');
end
if length(varargin)==1
	handles.fname                   = varargin{1};
	handles.fname                   = pa_fcheckext(handles.fname,'.dat');
	handles.fname                   = pa_fcheckexist(handles.fname);
end

% Check for button-calibration
% Check for number of channels
%

% Trial, stimulus and channel information
handles.csvfile                     = pa_fcheckext(handles.fname,'csv');
[expinfo,chaninfo,mLog]             = pa_readcsv(handles.csvfile);
handles.Nsamples                    = chaninfo(1,6);
handles.Fsample                     = chaninfo(1,5);
handles.Ntrial                      = max(mLog(:,1));
handles.StimType                    = mLog(:,5);
handles.StimOnset                   = mLog(:,8);
handles.Stim                        = pa_log2stim(mLog);
handles.NChan                       = expinfo(1,8);
DAT                                 = pa_loaddat(handles.fname, handles.NChan, handles.Nsamples, handles.Ntrial);
% Traces
handles.hortrace                    = squeeze(DAT(:,:,1));
handles.vertrace                    = squeeze(DAT(:,:,2));
handles.fronttrace                  = squeeze(DAT(:,:,3));
% Average of traces
handles.H                           = mean(handles.hortrace);
handles.V                           = mean(handles.vertrace);
handles.F                           = mean(handles.fronttrace);


% Select Fixation LED (LED = type 0)
sel                                 = ismember(handles.Stim(:,3),0:1);
S                                   = handles.Stim(sel,:);
mn                                  = min(unique(S(:,2)));
sel                                 = ismember(S(:,2),mn);
handles.TarAz                       = S(sel,4);
handles.TarEl                       = S(sel,5);


% Network Properties
handles.NhiddenHor                  = 5;
set(handles.edit_nhidden_hor,'Value',handles.NhiddenHor);
set(handles.edit_nhidden_hor,'String',num2str(handles.NhiddenHor));
handles.NhiddenVer                  = 5;
set(handles.edit_nhidden_ver,'Value',handles.NhiddenVer);
set(handles.edit_nhidden_ver,'String',num2str(handles.NhiddenVer));
handles.netfname                       = pa_fcheckext(handles.fname,'net');
set(handles.edit_netfname,'String',handles.netfname);
set(handles.popupmenu_algorithm,'Value',2);
index_algorithm                      = get(handles.popupmenu_algorithm,'Value');
switch index_algorithm
	case 1
		handles.netalgorithm        = 'trainlm';
	case 2
		handles.netalgorithm        = 'trainbr';
end
handles.FixRem                      = [];
handles.Fix                         = 1:length(handles.H);
handles.zoom                        = false;
zoom off;
set(handles.check_front,'Value',1);

%%
% Obtain target locations and mean responses in vector style
C                                   = handles.Fix;
TarAz                               = handles.TarAz(C);
TarEl                               = handles.TarEl(C);
Az                                  = handles.H;
El                                  = handles.V;
Tar                                 = unique([round(TarAz) round(TarEl)],'rows');
Res                                 = ones(size(Tar));
for i                               = 1:length(Tar)
	sel                             = round(TarAz) == Tar(i,1) & round(TarEl) == Tar(i,2);
	Res(i,:)                        = [mean(Az(sel)) mean(El(sel))];
end
% create grids for each potential unique Target location
uEl                                 = unique(Tar(:,2));
nEl                                 = length(uEl);
lAz                                 = 1:nEl;
% obtain maximum number of azimuths for one elevation
for i                               = 1:nEl
	sel                             = Tar(:,2) == uEl(i);
	lAz(i)                          = length(Tar(sel,1));
end
% define NaN-grid based on maximum number of locations (Az*El)
mlAz                                = max(lAz);
AzGrid                              = NaN*ones(mlAz,nEl);
ElGrid                              = AzGrid;
% define response NaN-grid
AzRes                               = AzGrid;
ElRes                               = ElGrid;
% replace NaNs in Target Grids with actual locations
for i                               = 1:nEl
	sel                             = Tar(:,2) == uEl(i);
	Az                              = Tar(sel,1);
	lAz                             = length(Az);
	El                              = uEl(i)*ones(size(Az));
	indx                            = floor((mlAz+1)/2-(lAz-1)/2):((mlAz+1)/2+(lAz-1)/2);
	AzGrid(indx,i)                  = Az;
	ElGrid(indx,i)                  = El;
end
% replace NaNs in Response Grids with actual locations
for i                               = 1:size(AzGrid,1)
	for j                           = 1:size(AzGrid,2)
		if ~isnan(AzGrid(i,j))
			sel                         = Tar(:,1) == AzGrid(i,j) & Tar(:,2) == ElGrid(i,j) ;
			AzRes(i,j)                  = Res(sel,1);
			ElRes(i,j)                  = Res(sel,2);
		end
	end
end
handles.AzRes                       = AzRes;
handles.ElRes                       = ElRes;
handles.AzGrid                      = AzGrid;
handles.ElGrid                      = ElGrid;


%% Plot graphics
function handles                    = plotXYZvsAZEL(handles)
C                               = handles.Fix;
FixAz                           = handles.H(C);
FixEl                           = handles.V(C);
FixF							= handles.F(C);
TarAz                           = handles.TarAz(C);
TarEl                           = handles.TarEl(C);

%% XY AXES
axes(handles.axes_xy); %#ok<*MAXES>
cla; hold on;
Inipar(1)                              = -90;
Inipar(2)                              = 1/((max(FixEl)-min(FixEl))/2);
Inipar(3)                              = mean(FixEl);
Inipar(4)                              = 0;
Vpar                                   = pa_fitasin(FixEl,TarEl,Inipar);
CalV                                   = pa_asin(FixEl,Vpar);
CalV                                   = CalV(:);

Inipar(1)                               = -90;
Inipar(2)                               = 1/((max(FixAz)-min(FixAz))/2);
Inipar(3)                               = mean(FixAz);
Inipar(4)                               = 0;
Hpar                                    = pa_fitasin(FixAz,TarAz,Inipar);
CalH                                    = pa_asin(FixAz,Hpar);
CalH                                    = CalH(:);

% Plot AZEL boundaries
h                                         = plot([-90 0],[0 90],...
	[0 90],[90 0],...
	[90 0],[0 -90],...
	[0 -90],[-90 0]);set(h,'LineWidth',2,'Color','k')
hold on
h                                          = plot([0 0],[-90 90],'k');set(h,'LineWidth',2,'Color',[0.7 0.7 0.7])
h                                          = plot([-90 90],[0 0],'k');set(h,'LineWidth',2,'Color',[0.7 0.7 0.7])
% Plot Target
h=plot(handles.AzGrid,handles.ElGrid,'ko');set(h,'MarkerFaceColor','w','Color',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7]);
h=plot(handles.AzGrid',handles.ElGrid','ko');set(h,'MarkerFaceColor','w','Color',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7]);

h                                          = plot(TarAz,TarEl,'ko');set(h,'MarkerFaceColor','w');
plot(CalH,CalV,'k*');
box on; axis square;
xlabel('Azimuth (deg)');
ylabel('Elevation (deg)');
Az = [CalH TarAz]';
El = [CalV TarEl]';
plot(Az,El,'k-');
ax = [min(TarAz)-2 max(TarAz)+2 min(TarEl)-2 max(TarEl)+2];
axis(ax);

axes(handles.axes_theta);
cla; hold on;
plot(TarAz,FixAz,'k.');
hold on
plot(TarAz,FixF,'k.','Color',[.7 .7 .7]);
xlabel('\alpha_T (deg)');
ylabel('H and F Field (V)');
box on;
xlim([-90 90]);

axes(handles.axes_phi);
cla; hold on;
cla; hold on;
plot(TarEl,FixEl,'k.');
hold on
plot(TarEl,FixF,'k.','Color',[.7 .7 .7]);
xlabel('\epsilon_T (deg)');
ylabel('V and F Field (V)');
set(gca,'YAxisLocation','right');
box on;
xlim([-90 90]);

%% NETWORK GRAPHICS
function handles    = plottrain(handles)
% Initialization
C					= handles.Fix; % Selected targets
TarAz               = handles.TarAz(C);
TarEl               = handles.TarEl(C);
H					= handles.H(C);
V					= handles.V(C);
F					= handles.F(C);

fro                 = get(handles.check_front,'Value');
if fro
	AD              = [H;V;F];
else
	AD              = [H;V];
	Hfit	= linspace(min(H),max(H),20);
	ADfit	= [Hfit; 2*std(V)*randn(size(Hfit))+mean(V)];
	Afit    = sim(handles.hnet,ADfit)';
	
	Vfit	= linspace(min(V),max(V),20);
	ADfit	= [2*std(H)*randn(size(Vfit))+mean(H); Vfit];
	Efit    = sim(handles.vnet,ADfit)';
	axes(handles.axes_theta);
	plot(Afit,Hfit,'r-','LineWidth',2,'Color',[1 .5 .5]);
	
	axes(handles.axes_phi);
	plot(Efit,Vfit,'r-','LineWidth',2,'Color',[1 .5 .5]);
end
%% TEST
[A,E]	= calib(AD,handles);

axes(handles.axes_theta);
[Atmp,indx] = sort(A);
plot(A(indx),H(indx),'b-','LineWidth',2);

axes(handles.axes_phi);
[Etmp,indx] = sort(E);
plot(E(indx),V(indx),'r-','LineWidth',2);

mfixa           = mean(A-TarAz);
stfa            = std(A-TarAz);
mfixe           = mean(E-TarEl);
stfe            = std(E-TarEl);

axes(handles.axes_net);
cla;
plot(TarAz,A-TarAz,'b.');
hold on;
plot(TarEl,E-TarEl,'r.');
lsline; axis square; pa_horline(0,'k:');
xlabel('Input (deg)','FontSize',12); ylabel('Error (deg)','FontSize',12);
str = {sprintf('Mean \\alpha %0.2f +/- %0.2f',mfixa,stfa),sprintf('Mean \\epsilon %0.2f +/- %0.2f',mfixe,stfe)};
title(str);
axis([-90 90 min([A-TarAz;E-TarEl]) max([A-TarAz;E-TarEl])]);

%% NET2
fro                             = get(handles.check_front,'Value');
if fro
	
	h       = AD(1,:);
	v       = AD(2,:);
	f       = AD(3,:);
	X       = linspace(min(h),max(h),100);
	Y       = linspace(min(v),max(v),100);
	[XI,YI] = meshgrid(X,Y);
	ZI      = griddata(h,v,f,XI,YI,'cubic'); % Z is interpolated from H and V
	% Simulate network
	[m,n]     = size(XI);
	Hsim           = zeros(size(XI));
	Vsim           = zeros(size(YI));
	for i       = 1:m
		[Hsim(i,:),Vsim(i,:)]	= calib([XI(i,:)' YI(i,:)' ZI(i,:)']',handles);
% 
% 		Hsim(i,:)  = sim(handles.hnet,[XI(i,:)' YI(i,:)' ZI(i,:)']');
% 		Vsim(i,:)  = sim(handles.vnet,[XI(i,:)' YI(i,:)' ZI(i,:)']');
	end;
	% Conversion to AD numbers
	ADpost    = [XI(:)';YI(:)';ZI(:)'];
	adh         = ADpost(1,:);
	adh         = reshape(adh,m,n);
	adv         = ADpost(2,:);
	adv         = reshape(adv,m,n);
	% Next values should be changed according to your own set-up
	c           = -90:3:90;
	ctheta      = -90:10:90;
	cphi = -90:15:90;
	
	axes(handles.axes_net2);
	cla;
	
	% General iso-azimuth and iso-elevation contour-lines
	contour(adh,adv,Hsim,c);
	hold on;
	contour(adh,adv,Vsim,c);
	
	[c2,h2]     = contour (adh,adv,Vsim,cphi,'b-'); set(h2,'LineWidth',2); %#ok<*ASGLU>
	[c2,h2]     = contour (adh,adv,Hsim,ctheta,'b-'); set(h2,'LineWidth',2);
	[c2,h2]     = contour (adh,adv,Vsim,[-90 -90],'b-'); set(h2,'LineWidth',4);
	[c2,h2]     = contour (adh,adv,Vsim,[90 90],'b-'); set(h2,'LineWidth',4);
	[c2,h2]     = contour (adh,adv,Hsim,[-90 -90],'b-'); set(h2,'LineWidth',4);
	[c2,h2]     = contour (adh,adv,Hsim,[90 90],'b-'); set(h2,'LineWidth',4);
	[c3,h3]     = contour (adh,adv,Hsim,[0 0],'r-'); set(h3,'LineWidth',2);
	[c4,h4]     = contour (adh,adv,Vsim,[0 0],'r-'); set(h4,'LineWidth',2);
	plot(H,V,'k.','MarkerSize', 15);
	xlabel ('X');
	ylabel ('Y');
	axis square;
end
%% NET3

axes(handles.axes_net3);
cla;
h           = plot([-90 0],[0 90],...
	[0 90],[90 0],...
	[90 0],[0 -90],...
	[0 -90],[-90 0]);set(h,'LineWidth',2,'Color','k')
hold on
% Plot Target
h=plot(handles.AzGrid,handles.ElGrid,'ko');set(h,'MarkerFaceColor','w','Color',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7]);
h=plot(handles.AzGrid',handles.ElGrid','ko');set(h,'MarkerFaceColor','w','Color',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7]);
h           = plot([0 0],[-90 90],'k');set(h,'LineWidth',2)
h           = plot(TarAz,TarEl,'ko');set(h,'MarkerFaceColor','w');
plot(A,E,'k*');
Az          = [A TarAz]';
El          = [E TarEl]';
plot(Az,El,'k-');
box on; axis square;
ax = [min(TarAz)-5 max(TarAz)+5 min(TarEl)-5 max(TarEl)+5];
axis(ax);

xlabel('Azimuth (deg)');
ylabel('Elevation (deg)');
set(gca,'YaxisLocation','right');


%% Train Calibration Network
function handles                = trainnet(handles)
handles                         = fithor(handles);
handles                         = fitver(handles);

%% TRAIN horizontal NETWORK
function handles            = fithor(handles)
warning('off','MATLAB:singularMatrix');
% Input and target vectors
FixC                             = handles.Fix;
p(1,:)                           = handles.H(FixC);        % First row = horizontal response
p(2,:)                           = handles.V(FixC);        % Second row = vertical response
fro                              = get(handles.check_front,'Value');
if fro
	p(3,:)                       = handles.F(FixC);        % Second row = vertical response
end
t(1,:)                                = handles.TarAz(FixC)';    % Horizontal target position

[p,pmap]					= mapminmax(p);
[t,tmap]					= mapminmax(t);

% Initialize feedforward-network
% Default learning procedure for ff-networks is LM with BR
% this can be changed in handles.netalgorithm
net                         = newff(p,t,[handles.NhiddenHor],{'tansig','purelin'},handles.netalgorithm);
net.DivideFcn = '';

% Train network with (default) parameters
net.trainParam.epochs       = 500;
net.trainParam.goal         = 0.00001;
net.trainParam.min_grad     = 0.00001;
% net.trainParam.show         = 50;
net.trainParam.show = NaN;

net                         = train(net,p,t);

handles.hnet                = net;
handles.hpmap				= pmap;
handles.htmap				= tmap;

%% TRAIN vertical NETWORK
function handles            = fitver(handles)
warning('off','MATLAB:singularMatrix');

% Input and target vectors
FixC                        = handles.Fix;
p(1,:)                      = handles.H(FixC);        % First row = horizontal response
p(2,:)                      = handles.V(FixC);        % Second row = vertical response
fro = get(handles.check_front,'Value');
if fro
	p(3,:)                           = handles.F(FixC);        % Second row = vertical response
end
t                            = handles.TarEl(FixC)';    % Horizontal target position

[p,pmap]					= mapminmax(p);
[t,tmap]					= mapminmax(t);

% Initialize feedforward-network
net                         = newff(p,t,[handles.NhiddenVer],{'tansig','purelin'},handles.netalgorithm);

net.DivideFcn = '';
% Train network with (default) Levenberg-Maquard method
net.trainParam.epochs       = 500;
net.trainParam.goal         = 0.0001;
net.trainParam.min_grad     = 0.0001;
% net.trainParam.show         = 50;
net.trainParam.show = NaN;

net                         = train(net,p,t); % Default learning procedure for ff-networks is LM

handles.vnet                = net;
handles.vpmap				= pmap;
handles.vtmap				= tmap;


%% --- Outputs from this function are returned to the command line.
function varargout = pa_calibrate_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btn_train.
function btn_train_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to btn_train (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes_net);cla;
axes(handles.axes_net2);cla;
axes(handles.axes_net3);cla;
handles = trainnet(handles);
plottrain(handles);



function edit_nhidden_hor_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nhidden_hor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nhidden_hor as text
%        str2double(get(hObject,'String')) returns contents of edit_nhidden_hor as a double
handles.NhiddenHor = str2double(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function edit_nhidden_hor_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to edit_nhidden_hor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end



function edit_nhidden_ver_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nhidden_ver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nhidden_ver as text
%        str2double(get(hObject,'String')) returns contents of edit_nhidden_ver as a double
handles.NhiddenVer = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function edit_nhidden_ver_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nhidden_ver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu_algorithm.
function popupmenu_algorithm_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_algorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_algorithm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_algorithm
index_algorithm                      = get(handles.popupmenu_algorithm,'Value');
switch index_algorithm
	case 1
		handles.netalgorithm         = 'trainlm';
	case 2
		handles.netalgorithm         = 'trainbr';
end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_algorithm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_algorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in btn_savenet.
function btn_savenet_Callback(hObject, eventdata, handles)
% hObject    handle to btn_savenet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hnet			= handles.hnet; %#ok<NASGU>
vnet			= handles.vnet; %#ok<NASGU>

net.hor.net		=  handles.hnet;
net.hor.pmap	=  handles.hpmap;
net.hor.tmap	=  handles.htmap;
net.ver.net		=  handles.vnet;
net.ver.pmap	=  handles.vpmap;
net.ver.tmap	=  handles.vtmap; %#ok<STRNU>

netfname		= handles.netfname;
save(netfname,'net');
str1	= get(handles.display,'String');
str2	= ['Save ' upper(handles.netfname)];
str		= char({str1;str2});
set(handles.display,'String',str);
		
function edit_netfname_Callback(hObject, eventdata, handles)
% hObject    handle to edit_netfname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_netfname as text
%        str2double(get(hObject,'String')) returns contents of edit_netfname as a double

handles.netfname = get(hObject,'String');
handles.netfname = fcheckext(handles.netfname,'.net');
set(handles.edit_netfname,'String',handles.netfname);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_netfname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_netfname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_browse.
function btn_browse_Callback(hObject, eventdata, handles)
% hObject    handle to btn_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





function edit_nspeaker_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nspeaker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nspeaker as text
%        str2double(get(hObject,'String')) returns contents of edit_nspeaker as a double
handles.nspeaker = str2double(get(hObject,'String'));
handles.speakernr                  = repmat(1:handles.nspeaker:29,7,1);
for i = 1:7
	handles.speakernr(i,:)                  = handles.speakernr(i,:)+29*(i-1);
end
handles.speakernr = handles.speakernr(:);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_nspeaker_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nspeaker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in display.
function display_Callback(hObject, eventdata, handles)
% hObject    handle to display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns display contents as cell array
%        contents{get(hObject,'Value')} returns selected item from display


% --- Executes during object creation, after setting all properties.
function display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_calfiles.
function btn_calfiles_Callback(hObject, eventdata, handles)
% hObject    handle to btn_calfiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str1	= get(handles.display,'String');
str2	= 'Calibrate data-files';
str		= char({str1;str2});
set(handles.display,'String',str);
drawnow;
pa_calfiles;
str1	= get(handles.display,'String');
str2	= 'Calibrating finished';
str		= char({str1;str2});
set(handles.display,'String',str);
drawnow;

str1	= get(handles.display,'String');
str2	= 'Filtering calibrated traces';
str		= char({str1;str2});
set(handles.display,'String',str);
drawnow;
pa_hvfilt;
str1	= get(handles.display,'String');
str2	= 'Filtering finished';
str		= char({str1;str2});
set(handles.display,'String',str);
drawnow;

function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function traincal_KeyPressFcn(hObject, eventdata, handles)
% h       = get(handles.figure1);
% a       = h.CurrentCharacter;

function traincal_BtnDownFcn(hObject, eventdata, handles)
h       = get(handles.axes_xy);
ax      = axis;
pnt     = h.CurrentPoint;
pnt     = pnt(1,[1 2]);
if pnt(1)< ax(2) && pnt(1) > ax(1) && pnt(2)< ax(4) && pnt(2) > ax(3)
	axes(handles.axes_xy)
	plot(pnt(1),pnt(2),'ro');
	[mndist,indx]       = min(hypot(handles.TarAz-pnt(1),handles.TarEl-pnt(2)));
	
	handles.FixRem      = [handles.FixRem;indx];
	handles.FixRem      = unique(handles.FixRem);
	indx                = 1:length(handles.H);
	handles.Fix         = setxor(indx,handles.FixRem);
	plotXYZvsAZEL(handles);
end


% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in btn_zoom.
function btn_zoom_Callback(hObject, eventdata, handles)
% hObject    handle to btn_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes_xy);
if ~handles.zoom
	set(handles.btn_zoom,'String','Zoom off');
	zoom on;
	handles.zoom = true;
elseif handles.zoom
	set(handles.btn_zoom,'String','Zoom on');
	zoom off;
	handles.zoom = false;
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on selection change in popupmenu_azimuth.
function popupmenu_azimuth_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_azimuth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_azimuth contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_azimuth
val             = get(handles.popupmenu_azimuth,'Value');
string_list     = get(handles.popupmenu_azimuth,'String');
selected_string = string_list(val,:);
handles.azimuthspoke = str2double(selected_string);
handles = plotXYZvsAZEL(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_azimuth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_azimuth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_elevation.
function popupmenu_elevation_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_elevation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_elevation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_elevation

val             = get(handles.popupmenu_elevation,'Value');
string_list     = get(handles.popupmenu_elevation,'String');
selected_string = string_list(val,:);
handles.elevationspoke = str2double(selected_string);
handles = plotXYZvsAZEL(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_elevation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_elevation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in check_front.
function check_front_Callback(hObject, eventdata, handles)
% hObject    handle to check_front (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_front


function [A,E]	= calib(AD,handles)
ADh             = mapminmax('apply',AD,handles.hpmap);
A               = sim(handles.hnet,ADh)';
A               = mapminmax('reverse',A,handles.htmap);

ADv             = mapminmax('apply',AD,handles.vpmap);
E               = sim(handles.vnet,ADv)';
E               = mapminmax('reverse',E,handles.vtmap);
