function Log = pa_readlog (fname)
%
% LOG = READLOG(FNAME)
%
% Extract experimental paramters from a log file FNAME.log. File formats
% are used at the Biophysics department of the Donders Institute. This is
% used for backwards-compatibility.
%
% See also PA_TDT2DAT, PA_READSTIM

% 2011 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com


%% Initialization
if nargin<1
    fname = [];
end
LogFile    = pa_fcheckexist(fname);
LogFile    = pa_fcheckext(LogFile,'.log');

%% Load Data
FpLog       = fopen(LogFile,'rt','l');
if FpLog == -1,
  error(['Cannot open: ' LogFile]);
end;

Log = logtype;

% get header information
Log.Name     = fgetpar(FpLog,'Exp.ExpName', '%s');

Log.NPresent = fgetpar(FpLog,'Exp.Npresent','%d');
Log.NTest    = fgetpar(FpLog,'Exp.Ntest','%d');
Log.NTarget  = fgetpar(FpLog,'Exp.Ntarget','%d');
Log.NTrial   = Log.NTest*Log.NPresent;

Log.FSample  = fgetpar(FpLog,'Sig.Fsample','%d');
Log.NSample  = fgetpar(FpLog,'Sig.Nsample','%d');
Log.NChan    = fgetpar(FpLog,'Sig.Nchannel','%d');
Log.ChanStr  = fgetpar(FpLog,'Sig.SigType','%s');
Log.TSample  = 1000/Log.FSample;

Log.ChanEh = itempos (Log.ChanStr, 'Eh', '-');
Log.ChanEv = itempos (Log.ChanStr, 'Ev', '-');
Log.ChanHh = itempos (Log.ChanStr, 'Hh', '-');
Log.ChanHv = itempos (Log.ChanStr, 'Hv', '-');
Log.ChanBh = itempos (Log.ChanStr, 'Bh', '-');
Log.ChanBv = itempos (Log.ChanStr, 'Bv', '-');
Log.ChanS1 = itempos (Log.ChanStr, 'S1', '-');
Log.ChanS2 = itempos (Log.ChanStr, 'S2', '-');

fclose(FpLog);


function param = fgetpar(fp, par, format)

% FUNCTION param = fgetpar(fp, par, format)
%
%   read parameter(s) "par" form text file fp   
%   the search procedure for parameter names
%   is case sensitive.
%   
%   Jeroen Goossens

% search from beginning 
frewind(fp);
par = ['#' par];

while ~feof(fp)
  % search par identifier 
  ident=fscanf(fp,'%s',1);
  if strcmp(ident,par)==1
    % find '=' sign 
    symb=fscanf(fp,'%s',1);
    if strcmp(symb,'=')==1
      % read value(s)
      param = fscanf(fp,format,1);
      return
    end
    return
  end
end

if ~exist('param', 'var'), param=[]; end



function l = logtype
% FUNCTION L = LOGTYPE
%
%   Structure holding log information 
%   
% .. Paul ..
%

l = struct ('Name', '', ...
            'NTest',   [], 'NPresent', [], 'NTarget',  [], ...
            'NTrial',  [], ...
            'FSample', [], 'NSample',  [], 'NChan', [], ...
            'ChanStr', [], 'ChanX',  [], 'ChanY',   []);

  

function n = itempos (Str, ItemStr, SepChar)
%
% function n = itempos (str, itemstr <,sepchr>)
%
% Return the entry index of itemstr in the 
% string str in which items are separated by 
% separation characters (default '-')
%
% Ex. >>  n = itempos ('Eh-S1-Ev-S2', 'S1')
%     yields n = 2
%     
% .. Paul ..

if (nargin < 3)
  SepChar = '-';
end;

SepPos = find (Str == SepChar);
NItem = length(SepPos) + 1;

ItemStrPos = findstr (Str, ItemStr);

if isempty (ItemStrPos)
  n = 0;
else
  n = 1 + sum (SepPos < ItemStrPos);
end;


