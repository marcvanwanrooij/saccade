function Sac = lastsac(Sac)
% SEL = lastsac(Sac)
%   Select vector for last saccades
%
%   Marc van Wanrooij

%% Initialization
Tnr 	= Sac(:,1);
sel		= zeros(size(Sac,1),1);

% find index for saccades in certain trial
for i=Tnr',

	indx	= find( Sac(:,1)==i,1,'last');
	sel(indx)	= 1;
end;

%% make logical
sel = logical(sel);
Sac = Sac(sel,:);

