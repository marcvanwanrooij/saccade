function [AZ,EL] = pa_calib(AD1,AD2,AD3,S)
% Calibrate raw HVF-data traces contained in one matrix
%
% [AZ EL] = CALIB(AD,S)
%
% Calibrate the 3 raw data traces in AD (horizontal H, vertical V, frontal F). 
% S is either the NET file containing the neural network structure or the
% structure itself.
%
% Optionally, you can supply the various traces separately:
% [AZ EL] = CALIB(H,V,F,S)
%
%   See also PA_CALIBRATE
%
% Author: MW

%% Initialization
if nargin<3
    S           = AD2;
    AD          = AD1;
    [M,N]       = size(AD); %#ok<ASGLU>
    M           = 1;
else
    [M,N]       = size(AD1);
    AD          = [AD1(:) AD2(:) AD3(:)]';
end
if ischar(S)
    S           = fcheckexist(S,'.net');
    S           = load(S,'-mat');
end


%% NETWORK
ADh             = mapminmax('apply',AD,S.net.hor.pmap);
AZ              = sim(S.net.hor.net,ADh)';
AZ              = mapminmax('reverse',AZ,S.net.hor.tmap);

ADv             = mapminmax('apply',AD,S.net.ver.pmap);
EL              = sim(S.net.ver.net,ADv)';
EL              = mapminmax('reverse',EL,S.net.ver.tmap);

%% Reshape back
AZ              = reshape(AZ,M,N);
EL              = reshape(EL,M,N);


