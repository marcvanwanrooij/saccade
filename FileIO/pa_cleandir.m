% Clean data-directory
%
% PA_CLEANDIR
%
% Zips files in current data-directory XX-SS-YYYY-MM-DD to:
% XX-SS-YYYY-MM-DDdat.zip
% XX-SS-YYYY-MM-DDsac.zip
% XX-SS-YYYY-MM-DDmat.zip
% And deletes all non-mat and non-zip files.
% Essentially 'cleans' the data-directory of unnecessary files


% 2012 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

%% Directory
wd          = cd;
str         = findstr(wd,'\');
str         = max(str);
session     = wd(str+1:length(wd));
disp(['   Current directory is: ' upper(session)])
d = dir;

%% Zip interesting stuff
eval(['!zip ' session 'dat *.dat *.csv *.log *.cfg *.exp'])
eval(['!zip ' session 'mat *.mat'])
eval(['!zip ' session 'sac *.sac *.fix *.net'])
if isunix
    eval(['!unzip -t ' session])
    eval(['!unzip -t ' session 'mat'])
    eval(['!unzip -t ' session 'sac'])
end

%% Clean directory
decision        = input('Continue with removing the non-essential files? (y or n): ','s');
if strcmp(decision,'n')
    disp('Too bad')
    return
elseif strcmp(decision,'y')
    disp('   Removing sac, net, fix, dat, csv, hv, and log-files...')
    disp('   Keeping zip, txt and mat-files')
    if isunix
        !rm *.dat *.sac *.fix *.log *.hv *.net
    elseif ispc
        !del *.dat *.sac *.fix *.log *.hv *.net *.cfg *.exp
    end
end

clear wd decision
