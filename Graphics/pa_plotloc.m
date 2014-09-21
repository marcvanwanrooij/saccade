function pa_plotloc(input,varargin)

% PA_PLOTLOC(INPUT)
%   Plots the response vs stimulus for azimuth and elevation.
%   INPUT should be either a (split-off of a) SupSac-file or
%   Location Data ([Stim_az Resp_az Stim_el Resp_el]).
%
% Optional input arguments:
% PA_PLOTLOC(...,'WHICH',WHICH)
%   Plot either azimuth ('az') or elevation ('el') is plotted. Default is
%   both.
% PA_PLOTLOC(...,'RANGE',RANGE)
%   RANGE is used to indicate the range of the stimuli. 
%  Default value is [-60 60 -60 60].
%
% See also PA_LOC


% 2011  Marcus
% e-mail: marcvanwanrooij@neural-code.com

%% Optional arguments:
% simple regression (default) or a robust fitting procedure
which         = pa_keyval('which',varargin);
range         = pa_keyval('range',varargin);
if isempty(range)
	range = [-90 90 -90 90];
end

%%
if (size(input,2)==4)
    LocData           = input;
else
    LocData          = input(:,[23 8 24 9]);
end

%% Plot
if isempty(which)
    if numel(unique(LocData(:,1)))>2
        if numel(unique(LocData(:,3)))>2
            subplot(121);
        end
        pa_loc(LocData(:,1), LocData(:,2),'az',range);
        xlabel('Stimulus azimuth (deg)');
        ylabel('Response azimuth (deg)');
        axis square
    end
    if numel(unique(LocData(:,3)))>2
        if numel(unique(LocData(:,1)))>2
            
            subplot(122);
        end
        pa_loc(LocData(:,3), LocData(:,4),'el',range);
        if numel(unique(LocData(:,1)))>2
            set(gca,'YAxisLocation','right');
        end
        xlabel('Stimulus elevation (deg)');
        ylabel('Response elevation (deg)');
        axis square
    end
elseif strcmpi(which,'az')
    pa_loc(LocData(:,1), LocData(:,2),which,range);
elseif strcmpi(which,'el')
    pa_loc(LocData(:,3), LocData(:,4),which,range);
end


