% Example of how to run the OPS
% Load data(X and y).
load('test.mat')

% Generating options
options.window = 10; % number of initial variables in OPS (window)
options.increment = 5; % increment of variable in OPS
options.percentage = 100; % percent of variables investigated in OPS
options.crossvalidation = 'vblinds'; % cross-validation (CV) option
options.split = 10; % splits of CV
options.preprocessing = {'mean' 'mean'}; % X and y pre-processing
options.criteria = [0.02 10]; % stop criteria 
options.var = 50; % interval size
options.iOPS_type = 'autoOPS'; % run autoOPS in iOPS
options.calhops = 'no'; % hops will not be calculated in intervals

% Choose vector option
vector = 'main';

% Calcuates hmod and hops
% The user plays a major role in the choice of hmod
[hmod,hops] = nlvOPS(X,y,[],vector,options);
% If hmod is already knew:
% hmod=6;
% [hmod,hops] = nlvOPS(X,y,hmod,vector,options);

% Choose OPS approaches
type = {'autoOPS' 'feedOPS' 'iOPS'};

% Run OPS
OPS = newOPS(X,y,hmod,hops,type,vector,options);
