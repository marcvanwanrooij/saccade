function [SAC,N] = pa_loadsac(fname)
% Load saccade markers
%
% [SAC,N]=LOADSAC(fname)
%   load saccade marks (on and offsets) from SAC files
%   SAC(1,:) = saccade onsets
%   SAC(2,:) = saccade offsets
%   N        = Number of saccades
%
%
% See also ULTRADET, HVTOMAT
%
% Jeroen Goossens
% Marc van Wanrooij 2007
if nargin<1
    fname = [];
end
%% Check
fname       = pa_fcheckexist(fname,'.sac');
fname       = pa_fcheckext(fname,'.sac');

%% Load
fid         = fopen(fname,'r','l');
[SAC,N]     = fread(fid,inf,'ulong');
SAC         = reshape(SAC,length(SAC)/2,2)';

% [SAC,N]     = fread(fid,[2,inf],'ulong');

fclose(fid);

%% 
if N==0,
   SAC      = [];
   return;
end;

%% Sort
if size(SAC,2),
  [dum,indx]= sort(SAC(1,:));
  SAC       = SAC(:,indx);
end;

%% Delete SACcades with size < 0
SAC         = SAC(:, SAC(1,:)<SAC(2,:) );
N           = size(SAC,2);

