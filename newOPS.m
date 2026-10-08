function OPSresult = newOPS(X,y,hmod,hops,type,vector,options)
%---------------------------------------------------------------------------
%  The NEWOPS algorithm performs the variable selection using the new three
%  approaches of OPS: autoOPS, feedOPS and iOPS.
%
%  autoOPS is an automatic method to perform the variable selection using
%  several informative vectors and their combinations, providing the best
%  result.
%
%  feedOPS is a method to apply the variable selection using a feedback OPS,
%  wherein pre-selected variables can return to a new selection until
%  achieving specific criteria.
%
%  iOPS is a method to perform the variable selection using autoOPS or
%  feedOPS in each subdivision of the full array.
%---------------------------------------------------------------------------
% Input:
%
%    X (m,n): Independent variable.
%
%    y (m,1): Dependent variable.
%       hmod: (OPTIONAL) Number of latent variables of model.
%       hops: (OPTIONAL) Number of latent variables used by OPS algorithm.
%              Vector with 1 to 5 elements. [reg sqr vip nas wght] are the
%              vector options that require hops. If one of them will be
%              used, hops will be one element. If more than one vector will
%              be used, please, enter with five hops.
%
%       type: (OPTIONAL) String or cell with the OPS option. Default option run
%       only autoOPS algorithm.
%          Default: 'autoOPS'
%          OBS: Please, use a cell to enter with more than one option.
%          Type options:
%           'autoOPS': run autoOPS;
%           'feedOPS': run feedOPS;
%           'iOPS': run feedOPS;
%           {'autoOPS','feedOPS'}: run autoOPS and feedOPS;
%           {'autoOPS','iOPS'}: run autoOPS and iOPS;
%           {'feedOPS','iOPS'}: run feedOPS and iOPS;
%           {'autoOPS','feedOPS','iOPS'}: run autoOPS, feedOPS and iOPS;
%
%     vector: (OPTIONAL) String with the vector option. Default option run
%     all vectors individually: 'main' option.
%      Individual vectors:
%         'reg': regression.
%         'cor': correlation.
%         'sqr': residual.
%         'vip': variable influence on projection.
%         'nas': net analyte signal.
%        'urxy': univariate regression xy.
%        'wght': weights.
%         'cov': covariance procedures.
%        'main': all vectors individually.
%       'inter': bynary combinations of vectors and product of all vectors.
%         'all': 'main' option and 'inter' option.
%    OBS -> Any binary combination of vectors is possible to be performed
%    individually. Please, attention to pairwise these vectors. The
%    combination follows a logic: from the vectors list above, 'reg' and
%    'cor' must be combined as 'regcor' and not as 'correg'. This example
%    shows that the logic follows the list order, with the first element of
%    the combination always been the first to appear in the list.
%
%    options: (OPTIONAL)
%      Window: {10}  window of initial variables.
%      Increment: {5} increment of variables to be added over the window.
%      Percentage: {100} percentage of the variables to be considered in
%                  the selection.
%      cv: [ 'full' | 'cblocks' | {'vblinds'} | 'random' ]
%          A cell containing standard cross-validation (CV).
%          cross-validation options:
%          'full': leave-one-out;
%          'cblocks': contiguous block;
%          'vblinds': venetian blinds;
%          'random': random subset).
%      Split: {10} number of subsets to divide data into for CV;
%      Preprocessing: A cell containing standard preprocessing structures
%      for the X- and Y- blocks.
%          { 'X preprocessing' 'y preprocessing' }
%          Default: {'mean' 'mean'}
%          Preprocessing options:
%           'mean': mean center;
%           'auto': autoscale;
%           'none': no preprocessing.
%      Criteria: Vector with 2 inputs [p nloop] used by feedOPS.
%             p: percentage difference between two consecutive RMSECVS.
%             nloop: number of loops.
%      Var: Number of variables included in each interval in iOPS.
%              Intervals less than 50 variables are not allowed.
%              Default option: 50 variables.
%      iOPS_type: String or cell with the OPS option. Default option run
%       only autoOPS algorithm.
%          Default: 'autoOPS'
%          OBS: Please, use a cell to enter with more than one option.
%          Type options:
%           {'autoOPS','feedOPS'}: run autoOPS and feedOPS;
%           'autoOPS': run autoOPS;
%           'feedOPS': run feedOPS;
%      Calhops: [ {'no'} | 'yes' ] The option to perform (or not) the
%            calculation of hOPS for each interval. If 'yes', the running
%            time is higher.
%
%      Options e.g.:
%        options.window=10;
%        options.increment=5;
%        options.percentage=100;
%        options.crossvalidation='random';
%        options.split=10;
%        options.preprocessing={'mean','mean'};
%        options.criteria = [0.02 10];
%        options.var = 50;
%        options.iOPS_type = 'autoOPS';
%        options.calhops = 'no';
%
%   Output:
%
%      Auto : Struct with the results: best informative vector, selected
%      variables,quality parameters (RMSECV and Rcv), and predict values of
%      Ycv.
%
%   I/O:  OPS = newOPS(X,y);
%
%         OPS = newOPS(X,y,hmod);
%
%         OPS = newOPS(X,y,hmod,hops);
%
%         OPS = newOPS(X,y,hmod,hops,type);
%
%         OPS = newOPS(X,y,hmod,hops,type,vector);
%
%         OPS = newOPS(X,y,hmod,hops,type,vector,options);
%
% See also: NLVOPS, SELOPS, CREATEOPTIONS
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%% OPTIONS

validation

if nargin < 7
    options.window = 10;
    options.increment = 5;
    options.percentage = 100;
    options.crossvalidation = 'vblinds';
    options.split = 10;
    options.preprocessing = {'mean' 'mean'};
    options.criteria = [0.02 10];
    options.var = 50;
    options.iOPS_type = 'autoOPS';
    options.calhops = 'no';
else
    if ~isfield(options,'window')
        options.window = 10;
    end
    if ~isfield(options,'increment')
        options.increment = 5;
    end
    if ~isfield(options,'percentage')
        options.percentage = 100;
    end
    if ~isfield(options,'crossvalidation')
        options.crossvalidation = 'vblinds';
    end
    if ~isfield(options,'split')
        options.split = 10;
    end
    if ~isfield(options,'preprocessing')
        options.preprocessing = {'mean' 'mean'};
    end
end

%%  VECTOR

if nargin < 6 || isempty(vector)
    display('Default VECTOR option is "main".');
    vector = 'main';
end

%% OPS APPROACHES

if nargin < 5 || isempty(type)
    display('AutoOPS will be performed!');
    type = {'auto'};
end

%% HMOD and HOPS

if nargin < 3
    hmod = calchmod(X,y,options,1);
    NhOPS = calchops(X,y,hmod,vector,options);
    hops = NhOPS.hops;
else
    if isempty(hmod)
        hmod = calchmod(X,y,options,1);
        if nargin < 4
            NhOPS = calchops(X,y,hmod,vector,options);
            hops = NhOPS.hops;
        else
            if isempty(hops)
                switch vector
                    case {'reg','sqr','vip','nas','wght','REG','SQR','VIP','NAS','WGHT'}
                        NhOPS = calchops(X,y,hmod,vector,options);
                        hops = NhOPS;
                    case {'main','MAIN','inter','INTER','all','ALL',}
                        NhOPS = calchops(X,y,hmod,vector,options);
                        hops = NhOPS.hops;
                    case {'regcor','regsqr','regvip','regnas','regurxy','regwght','regcov','corsqr','corvip','cornas','corurxy','corwght','corcov','sqrvip','sqrnas','sqrurxy','sqrwght','sqrcov','vipnas','vipurxy','vipwght','vipcov','nasurxy','naswght','nascov','urxywght','urxycov','wghtcov','prodall'}
                        NhOPS = calchopscomb(X,y,hmod,vector,options);
                        hops = NhOPS.hops;
                    case {'cor','urxy','cov','COR','URXY','COV','custom','Custom','CUSTOM'}
                        hops = hmod;
                    otherwise
                        display('Invalid option: "main" option is used.');
                        vector = 'main';
                        NhOPS = calchops(X,y,hmod,vector,options);
                        hops = NhOPS.hops;
                end
            else
                switch vector
                    case {'reg','sqr','vip','nas','wght','REG','SQR','VIP','NAS','WGHT'}
                        l = length(hops);
                        if l > 1 && l < 5
                            warning('Wrong length of hOPS, it will be recalculated.')
                            NhOPS = calchops(X,y,hmod,vector,options);
                            hops = NhOPS;
                        end
                    case {'main','MAIN','inter','INTER','all','ALL'}
                        l = length(hops);
                        if l ~= 5
                            warning('Some values are missing in the hOPS, it will be recalculated.')
                            NhOPS = calchops(X,y,hmod,vector,options);
                            hops = NhOPS.hops;
                        end
                    case {'regcor','regsqr','regvip','regnas','regurxy','regwght','regcov','corsqr','corvip','cornas','corurxy','corwght','corcov','sqrvip','sqrnas','sqrurxy','sqrwght','sqrcov','vipnas','vipurxy','vipwght','vipcov','nasurxy','naswght','nascov','urxywght','urxycov','wghtcov','prodall'}
                        l = length(hops);
                        if l ~= 5
                            warning('Some values are missing in the hOPS, it will be recalculated.')
                            NhOPS = calchopscomb(X,y,hmod,vector,options);
                            hops = NhOPS.hops;
                        end
                    case {'cor','urxy','cov','COR','URXY','COV','custom','Custom','CUSTOM'}
                        hops = hmod;
                    otherwise
                        display('Invalid option: "main" option is used.');
                        vector = 'main';
                        l = length(hops);
                        if l ~= 5
                            warning('Some values are missing in the hOPS, it will be recalculated.')
                            NhOPS = calchops(X,y,hmod,vector,options);
                            hops = NhOPS.hops;
                        end
                end
            end
        end
    else
        if nargin < 4
            NhOPS = calchops(X,y,hmod,vector,options);
            hops = NhOPS.hops;
        else
            if isempty(hops)
                switch vector
                    case {'reg','sqr','vip','nas','wght','REG','SQR','VIP','NAS','WGHT'}
                        NhOPS = calchops(X,y,hmod,vector,options);
                        hops = NhOPS;
                    case {'main','MAIN','inter','INTER','all','ALL',}
                        NhOPS = calchops(X,y,hmod,vector,options);
                        hops = NhOPS.hops;
                    case {'regcor','regsqr','regvip','regnas','regurxy','regwght','regcov','corsqr','corvip','cornas','corurxy','corwght','corcov','sqrvip','sqrnas','sqrurxy','sqrwght','sqrcov','vipnas','vipurxy','vipwght','vipcov','nasurxy','naswght','nascov','urxywght','urxycov','wghtcov','prodall'}
                        NhOPS = calchopscomb(X,y,hmod,vector,options);
                        hops = NhOPS.hops;
                    case {'cor','urxy','cov','COR','URXY','COV','custom','Custom','CUSTOM'}
                        hops = hmod;
                    otherwise
                        display('Invalid option: "main" option is used.');
                        vector = 'main';
                        NhOPS = calchops(X,y,hmod,vector,options);
                        hops = NhOPS.hops;
                end
            else
                switch vector
                    case {'reg','sqr','vip','nas','wght','REG','SQR','VIP','NAS','WGHT'}
                        l = length(hops);
                        if l > 1 && l < 5
                            warning('Wrong length of hOPS, it will be recalculated.')
                            NhOPS = calchops(X,y,hmod,vector,options);
                            hops = NhOPS;
                        end
                    case {'main','MAIN','inter','INTER','all','ALL'}
                        l = length(hops);
                        if l ~= 5
                            warning('Some values are missing in the hOPS, it will be recalculated.')
                            NhOPS = calchops(X,y,hmod,vector,options);
                            hops = NhOPS.hops;
                        end
                    case {'regcor','regsqr','regvip','regnas','regurxy','regwght','regcov','corsqr','corvip','cornas','corurxy','corwght','corcov','sqrvip','sqrnas','sqrurxy','sqrwght','sqrcov','vipnas','vipurxy','vipwght','vipcov','nasurxy','naswght','nascov','urxywght','urxycov','wghtcov','prodall'}
                        l = length(hops);
                        if l ~= 5
                            warning('Some values are missing in the hOPS, it will be recalculated.')
                            NhOPS = calchopscomb(X,y,hmod,vector,options);
                            hops = NhOPS.hops;
                        end
                    case {'cor','urxy','cov','COR','URXY','COV','custom','Custom','CUSTOM'}
                        hops = hmod;
                    otherwise
                        display('Invalid option: "main" option is used.');
                        vector = 'main';
                        l = length(hops);
                        if l ~= 5
                            warning('Some values are missing in the hOPS, it will be recalculated.')
                            NhOPS = calchops(X,y,hmod,vector,options);
                            hops = NhOPS.hops;
                        end
                end
            end
        end
    end
end

%% VARIABLE SELECTION

switch vector
    case {'custom','Custom','CUSTOM'}
        infvec = input('Insert the informative vector. Type the name here -> ','s');
        custom = evalin('base',char(infvec));
        options.custom = custom;
end
if ischar(type)
    type = cellstr(type);
end
s = size(type,2);

for j=1:s
    retype = type(j);
    retype = char(retype);
    switch retype
        case {'auto','AUTO','Auto','autoops','autoOPS','AUTOOPS'}
            Auto = autoOPS(X,y,hmod,hops,vector,options);
            OPSresult.autoOPS = Auto;
        case {'feed','FEED','Feed','feedOPS','FEEDOPS','feedops'}
            if ~isfield(options,'criteria')
                criteria = [0.02 10];
                options.criteria = criteria;
            end
            Feed = feedOPS(X,y,hmod,hops,vector,options);
            OPSresult.feedOPS = Feed;
        case {'iops','IOPS','iOPS'}
            if ~isfield(options,'criteria')
                criteria = [0.02 10];
                options.criteria = criteria;
            end
            if ~isfield(options,'var')
                options.var = 50;
            end
            if ~isfield(options,'iOPS_type')
                options.iOPS_type = 'auto';
            end
            if ~isfield(options,'calhops')
                options.calhops = 'no';
            end
            Inter = iOPS(X,y,hmod,hops,vector,options);
            OPSresult.iOPS = Inter;
        otherwise
            error('Invalid "type" option.');
    end
end
arrangefig

end

function [Xauto,meanX,stdX] = autoescale(X)

% CrossValidation Function OPS
%
% I/O: [Xauto,meanX,stdX] = autoescale(X);
%
% Copyright: Jussara V. Roque, 2018.

[m,~] = size(X);
meanX = mean(X);
stdX  = std(X);
a = stdX == 0;
Xauto = (X-meanX(ones(m,1),:))./stdX(ones(m,1),:);
Xauto(:,a)=0;

end

function [Xmean,meanX] = meancenter(X)

% CrossValidation Function OPS
%
% I/O: [Xmean,meanX] = meancenter(X);
%
% Copyright: Jussara V. Roque, 2017.

[m,~] = size(X);
meanX = mean(X);
Xmean = (X-meanX(ones(m,1),:));

end

function [RMSECV,RCV,Bias]=par(Yref,Ypred)

% CrossValidation Function OPS
%
% I/O: [[RMSECV,RCV,Bias]=par(Yref,Ypred);
%
% Copyright: Jussara V. Roque, 2017.

[n,m]=size(Yref);
RMSECV = sqrt( sum(sum((Ypred-Yref).^2))/(n*m) );
RCV = corrcoef(Yref,Ypred);
RCV = RCV(2,1);
Bias = sum(sum(Ypred-Yref))/(n*m);

end

function Xscaleback = scaleback(X,meanX,stdX)

% CrossValidation Function OPS
%
% I/O: Xscaleback = scaleback(X,meanX,stdX);
%
% Copyright: Jussara V. Roque, 2018.

[m,~] = size(X);
if nargin == 2
    Xscaleback = X + meanX(ones(m,1),:);
elseif nargin == 3
    Xscaleback = (X.*stdX(ones(m,1),:)) + meanX(ones(m,1),:);
end

end

function Xscalenew = scalenew(Xnew,meanXold,stdXold)

% CrossValidation Function OPS
%
% I/O: Xscalenew = scalenew(Xnew,meanXold,stdXold);
%
% Copyright: Jussara V. Roque, 2018.

[m,~] = size(Xnew);
if nargin == 2
    Xscalenew = (Xnew-meanXold(ones(m,1),:));
elseif nargin == 3
    Xscalenew = (Xnew-meanXold(ones(m,1),:))./stdXold(ones(m,1),:);
    a = stdXold == 0;
    Xscalenew(:,a)=0;
end

end

function [B,C,P,T,U,R,R2X,R2Y] = simplsorig(X,Y,A,S,XtX)

% SIMPLSORIG Full implementation of SIMPLS approach to PLS regression for (multivariate) Y.
%
%     Input:
%       X (m,n): Independent variable.
%       Y (m,p): Dependent variable.
%       A (1,1): Number of latent variables of model (hmod).
%
%    Output:
%       B (n,p): Regression coefficients.
%	    C (p,A): Y loadings
%	    P (n,A): X loadings
%	    T (m,A): X scores <standardized>
%	    U (m,A): Y scores
%	    R (n,A): X weights
%	    R2X (1,A): X-variance accounted for
%	    R2Y (1,A): Y-variance accounted for
%
% Copyright: Sijmen de Jong, 15/6/1997
% 	      Unilever Research Laboratorium, Vlaardingen, The Netherlands
%	      Copyright (c) 1997 for ChemoAC
%	      Dienst FABI, Vrije Universiteit Brussel
%	      Laarbeeklaan 103, 1090-Brussel Jette, Belgium
%
% Version: 1.1 (28/02/1998)
% Reference: S.de Jong, Chemom.Intell.Lab.Syst.,18 (1993) 251-263.

[~,px] = size(X);
[n,m] = size(Y); % size of the input data matrices

if nargin < 5
    S = [];
end

if isempty(S)  % if XtX not inputted, S=[]; always when S=[] then S=(Y'*X)'
    S = (Y'*X)'; 
end                
if nargin < 4 % if S is not inputted, XtX=[];
    XtX = [];
end				

if isempty(XtX) && n>3*px % when XtX=[] and X is very "tall", the booster XtX is calculated
    XtX = X'*X; 
end			

if nargin<3
    A=10;
end

A = min([A px n-1]); % if A is not inputted, then the defaul A is min[10 px n-1]
T = zeros(n ,A);
U = T; % initialization of variables
R = zeros(px,A); 
P = R; V = R;
C = zeros(m ,A);
R2Y = zeros(1,A);
z = zeros(m,1); 

StS = S'*S;	% SIMPLS algorithm
nm1 = n-1;
tol = 0;
for a = 1:A
    StS = StS-z*z';
    [Q,LAMBDA] = eig(StS);
    [lambda,j] = max(diag(LAMBDA));
    q = Q(:,j(1));
    r = S*q;
    t = X*r;
    if isempty(XtX)
        p = (t'*X)';
    else
        p = XtX*r;
    end
    if n>px
        d = sqrt(r'*p/nm1);
    else
        d = sqrt(t'*t/nm1);
    end
    tol=max(tol,d/1e5);
    v = p-V(:,1:max(1,a-1))*(p'*V(:,1:max(1,a-1)))';
    v = v/sqrt(v'*v);
    z = (v'*S)';
    S = S-v*z';
    % save results
    V(:,a) = v;
    R(:,a) = r/d; % X weights
    P(:,a) = p/(d*nm1); % X loadings
    T(:,a) = t/d; % X scores
    U(:,a) = Y*q; % Y scores
    C(:,a) = q*(lambda(1)/(nm1*d)); % Y loadings
    R2Y(1,a) =  lambda(1)/d; % Y-variance accounted for
end
clear StS V LAMBDA Q p q r t v z;
if d<tol,
    A=a-1;
    a=A;
    T=T(:,1:A);
    U=U(:,1:A);
    R=R(:,1:A);
    P=P(:,1:A);
    C=C(:,1:A);
end
while a>1
    U(:,a) = U(:,a)-T(:,1:a-1)*(U(:,a)'*T(:,1:a-1)/nm1)';
    a=a-1;
end
B = R*C'; % B-coefficients of the regression Y on X
if isempty(XtX)
    sumX2=sum(X.^2);
else
    sumX2 = sum(diag(XtX));
end
R2X = 100*nm1/sum(sumX2)*cumsum(sum(P.^2));
R2Y = 100/nm1/sum(sum(Y.^2))*cumsum(R2Y(1:A).^2);

end

function t=convtime(s)

T1=s/86400;
T2=floor(T1); %Days
T3=(T1-T2)*24;
T4=floor(T3); %Hours
T5=(T3-T4)*60;
T6=floor(T5); %Minutes
T7=(T5-T6)*60;
T8=floor(T7); %Seconds

t=[num2str(T2) ' days ' num2str(T4) ' hrs ' num2str(T6) ' min ' num2str(T8) ' s ']; 

end

function h = waitbarn(X,varargin)

%  WAITBAR a modified version of MATLAB's waitbar function.
%__________________________________________________________________________
%     H = WAITBAR(X,'message') creates and displays a waitbar of fractional
%           length X.  The handle to the waitbar figure is returned in H.
%           X should be between 0 and 1.
%
%     WAITBAR(X) will set the length of the bar in the most recently
%           created waitbar window to the fractional length X.
%
%     WAITBAR(X,H) will set the length of the bar in waitbar H
%           to the fractional length X.
%
%     WAITBAR(X,H,'message') will update the message text in
%           the waitbar figure, in addition to setting the fractional
%           length to X.
%
%     WAITBAR is typically used inside a FOR loop that performs a
%           lengthy computation.
%
%     Example:
%         h = waitbar(0,'Please wait...');
%         for i=1:1000,
%             % computation here %
%             waitbar(i/1000,h)
%         end
%
% NOTES:
% - This progarm produced with heavy modification of Chad English's timebar
% function.  The update was designed to recieve input identically to
% MATLAB's waitbar function to allow for interchangability.
%
% - This program does not apply the property values that the traditional
% waitbar allows.
%
%__________________________________________________________________________

% 1 - GATHER THE INPUT
if nargin == 1;
    h = findobj(allchild(0),'flat','Tag','waitbar');
    message = '';
elseif isnumeric(X) & ishandle(varargin{1}) & nargin == 2;
    h = varargin{1}; message = '';
elseif isnumeric(X) & ischar(varargin{1}) & nargin == 2;
    h = []; message = varargin{1};
elseif isnumeric(X) & ishandle(varargin{1}) & nargin == 3;
    h = varargin{1}; message = varargin{2};
else
    disp('Error defnining waitbar'); return;
end

% 2 - BUILD/UPDATE THE MESSAGE BAR
if  isempty(h) || ~ishandle(h(1)); h = buildwaitbar(X,message);
else updatewaitbar(h,X,message); end

end
%--------------------------------------------------------------------------
% SUBFUNCTION: buildwaitbar
function h = buildwaitbar(X,message)
% BUILDWAITBAR constructs the figure containing the waitbar

% 1 - SET WINDOW SIZE AND POSITION
% 1.1 - Gather screen information
screensize = get(0,'screensize');  % User's screen size
screenwidth = screensize(3);       % User's screen width
screenheight = screensize(4);      % User's screen height

% 1.2 - Define the waitbar position
winwidth = 350;           % Width of timebar window
winheight = 85;           % Height of timebar window
winpos = [0.5*(screenwidth-winwidth), ...
    0.5*(screenheight-winheight), winwidth, winheight];  % Position

% 2 - OPEN FIGURE AND SET PROPERTIES
wincolor = [1 1 1]; % Define window color

% 2.1 - Define the main waitbar figure
h = figure('menubar','none','numbertitle','off',...
    'name','0% Complete','position',winpos,'color',wincolor,...
    'tag','waitbar','IntegerHandle','off');

% 2.2 - Define the message textbox
userdata.text(1) = uicontrol(h,'style','text','hor','left',...
    'pos',[10 winheight-30 winwidth-20 20], 'string',message,...
    'backgroundcolor',wincolor,'tag','message');

% 2.3 - Build estimated remaining static text textbox
est_text = 'Estimated time remaining: ';
userdata.text(2) = uicontrol(h,'style','text','string',est_text,...
    'pos',[10 15 winwidth/2 20],'FontSize',7,...
    'backgroundcolor',wincolor,'HorizontalAlignment','right');

% 2.4 - Build estimated time textbox
userdata.remain = uicontrol(h,'style','text','string','',...
    'FontSize',7,'HorizontalAlignment','left',...
    'pos',[winwidth/2+10 14.5 winwidth-25 20], ...
    'backgroundcolor',wincolor);

% 2.5 - Build elapsed static text textbox
est_text = 'Total elapsed time: ';
userdata.text(3) = uicontrol(h,'style','text','string',est_text,...
    'pos',[10 3 winwidth/2 20],'FontSize',7,...
    'backgroundcolor',wincolor,'HorizontalAlignment','right');

% 2.6 - Build elapsed time textbox
userdata.elapse = uicontrol(h,'style','text','string','',...
    'pos',[winwidth/2+10 3.5 winwidth-25 20],'FontSize',7, ...
    'backgroundcolor',wincolor,'HorizontalAlignment','left');

% 2.7 - Build percent progress textbox
userdata.percent = uicontrol(h,'style','text','hor','right',...
    'pos',[winwidth-45 winheight-52 30 20],'string','',...
    'backgroundcolor',wincolor);

% 2.8 - Build progress bar axis
userdata.axes = axes('parent',h,'units','pixels','xlim',[0 1],...
    'pos',[10 winheight-45 winwidth-60 15],'box','on',...
    'color',[1 1 1],'xtick',[],'ytick',[]);

% 3 - INITILIZE THE PROGESS BAR
userdata.bar = ...
    patch([0 0 0 0 0],[0 1 1 0 0],[1 0.4 0.6]);  % Initialize  bar to zero area

tic

set(h,'userdata',userdata)               % Store data in thefigure
updatewaitbar(h,X,message);              % Updates waitbar if X~=0

set(userdata.remain,'string',convtime(0));

end
%--------------------------------------------------------------------------
% SUBFUNCTION: updatewaitbar
function updatewaitbar(h,progress,message)
% UPDATEWAITBAR changes the status of the waitbar progress

% 1 - GATHER WAITBAR INFORMATION
drawnow;                        % Needed for window to appear
h = h(1);                       % Only allow newest waitbar to update
userdata = get(h,'userdata');   % Get userdata from waitbar figure

% Check object tag to see if it is a timebar
if ~strcmp(get(h,'tag'), 'waitbar')
    error('Handle is not for a waitbar window')
end

% Update the message
if ~isempty(message);
    hh = guihandles(h);
    set(hh.message,'String',message);
end

% 2.2 - Calculate the estimated time remaining

sec_remain = convtime(toc*(1/progress-1));
e_mes = convtime(toc);
r_mes = sec_remain;

% 2.3 - Produce error if progress is > 1
if progress > 1; r_mes = 'Error, progress > 1'; end

% 2.4 - Update information
set(userdata.bar,'xdata',[0 0 progress progress 0]) % Update bar
set(userdata.remain,'string',r_mes); % Update remaining time string
set(userdata.elapse,'string',e_mes); % Update elapsed time string
set(userdata.percent,'string',...
    strcat(num2str(floor(100*progress)),'%')); % Update progress %
set(h,'Name',[num2str(floor(100*progress)),...
    '% Complete']); % Update figure name


set(gcf,'Resize','off')

end

function arrangefig

%% How many figures
FigHandle = findobj('Type','figure');
nFig = size(FigHandle,1);

% [x y w h]
P=get(gcf,'OuterPosition');

%% Monitor Size
scr_sz = get(0,'MonitorPositions');
scr_sz = scr_sz(1, :);

scn_h = (scr_sz(4) - P(4))/2;
scn_w = 50;


%% Arrange
for nCol=1:1:nFig
    
    fig_pos = [scn_w scn_h P(3) P(4)];
    
    set(FigHandle(nCol),'OuterPosition',fig_pos);
  
    scn_w = scn_w + (scr_sz(3)-P(3))/(nFig);
    
end

end

function autofig(X,resultauto)

% Graphic plots of AutoOPS!
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%
[~,n] = size(X);
var = 1:n;
figure
subplot(2,3,1) % Full DATA and selected variables
selOPS(X,resultauto.vs)
axis tight
xlabel ('Variables','FontSize',12)
ylabel ('Intensity','FontSize',12)
title('autoOPS','FontSize',12)

subplot(2,3,2) % Informative Vector
plot(var,resultauto.result.vect,'k')
axis tight
xlabel ('Variables','FontSize',12)
ylabel ('Intensity','FontSize',12)
title(char(resultauto.best_vec),'FontSize',12)

nsetb = resultauto.result.graph.nsetb;
rmsecv_b = resultauto.result.graph.rmsecv_b;
nlvs = resultauto.result.graph.nlvs;
nlv_b = resultauto.result.graph.nlv_b;
[~,h] = size(resultauto.result.graph.rmsecv_b);
s = 1:h;
subplot(2,3,3) % Variables subset
plot(nsetb,rmsecv_b,'-ko')
axis tight
text(nsetb+0.15,rmsecv_b,nlvs)
grid on
ylabel ('RMSECV','FontSize',12)
xlabel ('Variables subset','FontSize',12)
title('RMSECV values','FontSize',12)
hline(rmsecv_b(h),'-r');
text(min(s)+1,rmsecv_b(h),['RMSECV Full = ',num2str(rmsecv_b(h)),...
    ' nLV = ',num2str(nlv_b(h))],'verticalalignment','bottom','FontSize',10)

tabMP = table2array(resultauto.MP_results);
y = tabMP(:,2);
ycv = tabMP(:,3);
err = tabMP(:,5);
k = polyfit(y,ycv,1);
mi = min(min([y,ycv]));
ma = max(max([y,ycv]));

subplot(2,3,4.4) % Graphic y_measured x y_predicted in cross-validation
plot(y,ycv,'ro','MarkerFaceColor',[1 0 0])
hold on
plot(linspace(mi,ma,length(y)),linspace(mi,ma,length(y)),'--k')
axis tight
refplot(k);
xlabel ('Y Measured','FontSize',12)
ylabel ('Ycv Predicted','FontSize',12)
title('Correlation Y x Ycv','FontSize',12)
subplot(2,3,5.6) % Graphic of relative error (%)
bar(err,'k')
axis tight
hline([5,-5],'-r');
xlabel ('Samples','FontSize',12)
ylabel ('Relative Error (%)','FontSize',12)
title('Relative Error CV','FontSize',12)

end

function feedfig(X,resultfeed)

% Graphic plots of FeedOPS!
%
% Copyright: Jussara V. Roque, 2019.
% Check by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo. 
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019). 
% https://doi.org/10.1016/j.aca.2019.05.039

%%

figure
subplot(2,2,1) % Full DATA and selected variables
selOPS(X,resultfeed.vs)
axis tight
xlabel ('Variables','FontSize',12)
ylabel ('Intensity','FontSize',12)
title('feedOPS','FontSize',12)

%%
subplot(2,2,2) % Graphic FEED loops
result = table2array(resultfeed.Quality_parameters);
bar(result(2:end,4),'k')
axis tight
hline(result(1,4),'-r');
hold on
s = resultfeed.Loops.bestloop;
bar(s,result(s+1,4),'r')
z = length (resultfeed.Loops.vectors);
for i=1:z
    mvar = max(result(i+1,4));
    t = text(i-0.2,(mvar+(mvar*0.025)),char(resultfeed.Loops.vectors{i}));
    set(t,'rotation',45)
end
xlabel('Loops','FontSize',12)
ylabel('RMSECV','FontSize',12)
title('Feed Loops','FontSize',12)

%%
tabMP = table2array(resultfeed.MP_results);
y = tabMP(:,2);
ycv = tabMP(:,3);
err = tabMP(:,5);
k = polyfit(y,ycv,1);
mi = min(min([y,ycv]));
ma = max(max([y,ycv]));

subplot(2,2,3) % Graphic y_measured x y_predicted in cross-validation
plot(y,ycv,'ro','MarkerFaceColor',[1 0 0])
hold on
plot(linspace(mi,ma,length(y)),linspace(mi,ma,length(y)),'--k')
axis tight
refplot(k);
xlabel ('Y Measured','FontSize',12)
ylabel ('Ycv Predicted','FontSize',12)
title('Correlation Y x Ycv','FontSize',12)

%%
subplot(2,2,4) % Graphic of relative error (%)
bar(err,'k')
axis tight
hline([5,-5],'-r');
xlabel ('Samples','FontSize',12)
ylabel ('Relative Error (%)','FontSize',12)
title('Relative Error CV','FontSize',12)

end

function interfig(X,resultinter)

% Graphic plots of InterOPS!
%
% Copyright: Jussara V. Roque, 2019.
% Check by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%
if isfield(resultinter,'autoOPS')
    figure
    subplot(2,3,1) % Full DATA and selected variables of Ordered Interval Model
    selOPS(X,resultinter.autoOPS.vs)
    axis tight
    xlabel ('Variables','FontSize',12)
    ylabel ('Intensity','FontSize',12)
    title('Auto iOPS','FontSize',12)
    
    subplot(2,3,2) % RMSECV in intervals
    bar(resultinter.autoOPS.rmsecv_intervals,'k')
    axis tight
    hline(table2array(resultinter.autoOPS.Quality_parameters(1,4)),'-r');
    xlabel ('Intervals','FontSize',12)
    ylabel ('RMSECV','FontSize',12)
    title ('Intervals','FontSize',12)
    grid on
    
    nsetb = resultinter.autoOPS.auto.result.graph.nsetb;
    rmsecv_b = resultinter.autoOPS.auto.result.graph.rmsecv_b;
    nlvs = resultinter.autoOPS.auto.result.graph.nlvs;
    nlv_b = resultinter.autoOPS.auto.result.graph.nlv_b;
    [~,h] = size(resultinter.autoOPS.auto.result.graph.rmsecv_b);
    s = 1:h;
    
    subplot(2,3,3) % Variables subset
    rmsecvf = table2array(resultinter.autoOPS.Quality_parameters(1,4));
    plot(nsetb,rmsecv_b,'-ko')
    axis tight
    text(nsetb+0.15,rmsecv_b,nlvs)
    grid on
    ylabel ('RMSECV','FontSize',12)
    xlabel ('Variables subset','FontSize',12)
    title('RMSECV values','FontSize',12)
    hline(rmsecvf,'-r');
    text(min(s)+1,rmsecvf,['RMSECV Full = ',num2str(rmsecvf),...
        ' nLV = ',num2str(nlv_b(h))],'verticalalignment','bottom','FontSize',10)
    
    yOI = table2array(resultinter.autoOPS.MP_results(:,2));
    ycvOI = table2array(resultinter.autoOPS.MP_results(:,3));
    errOI = table2array(resultinter.autoOPS.MP_results(:,5));
    kOI = polyfit(yOI,ycvOI,1);
    mi = min(min([yOI,ycvOI]));
    ma = max(max([yOI,ycvOI]));
    
    % Graphic y_measured x y_predicted in cross-validation
    subplot(2,3,4.4)
    plot(yOI,ycvOI,'ro','MarkerFaceColor',[1 0 0])
    hold on
    plot(linspace(mi,ma,length(yOI)),linspace(mi,ma,length(yOI)),'--k')
    axis tight
    refplot(kOI);
    xlabel ('Y Measured','FontSize',12)
    ylabel ('Ycv Predicted','FontSize',12)
    title('Correlation Y x Ycv','FontSize',12)
    
    % Graphic of relative error (%)
    subplot(2,3,5.6)
    bar(errOI,'k')
    axis tight
    hline([5,-5],'-r');
    xlabel ('Samples','FontSize',12)
    ylabel ('Relative Error (%)','FontSize',12)
    title('Relative Error CV','FontSize',12)
else
    figure
    subplot(2,3,1) % Full DATA and selected variables of Ordered Interval Model
    selOPS(X,resultinter.feedOPS.vs)
    axis tight
    xlabel ('Variables','FontSize',12)
    ylabel ('Intensity','FontSize',12)
    title('Feed iOPS','FontSize',12)
    
    subplot(2,3,2) % RMSECV in intervals
    bar(resultinter.feedOPS.rmsecv_intervals,'k')
    axis tight
    hline(table2array(resultinter.feedOPS.Quality_parameters(1,4)),'-r');
    xlabel ('Intervals','FontSize',12)
    ylabel ('RMSECV','FontSize',12)
    title ('Intervals','FontSize',12)
    grid on

    subplot(2,3,3) % Final selection
    result = resultinter.feedOPS.feed.resultSEL;
    bar(result(:,4),'k')
    axis tight
    hline(table2array(resultinter.feedOPS.Quality_parameters(1,4)),'-r');
    hold on
    s = resultinter.feedOPS.feed.bestloop;
    bar(s,result(s,4),'r')
    z = length (resultinter.feedOPS.feed.vectors);
    for i=1:z
        mvar = max(result(i,4));
        t = text(i-0.2,(mvar+(mvar*0.025)),char(resultinter.feedOPS.feed.vectors{i}));
        set(t,'rotation',45)
    end
    xlabel('Loops','FontSize',12)
    ylabel('RMSECV','FontSize',12)
    title('Feed Loops','FontSize',12)
    
    yOI = table2array(resultinter.feedOPS.MP_results(:,2));
    ycvOI = table2array(resultinter.feedOPS.MP_results(:,3));
    errOI = table2array(resultinter.feedOPS.MP_results(:,5));
    kOI = polyfit(yOI,ycvOI,1);
    mi = min(min([yOI,ycvOI]));
    ma = max(max([yOI,ycvOI]));
    
    % Graphic y_measured x y_predicted in cross-validation
    subplot(2,3,4.4)
    plot(yOI,ycvOI,'ro','MarkerFaceColor',[1 0 0])
    hold on
    plot(linspace(mi,ma,length(yOI)),linspace(mi,ma,length(yOI)),'--k')
    axis tight
    refplot(kOI);
    xlabel ('Y Measured','FontSize',12)
    ylabel ('Ycv Predicted','FontSize',12)
    title('Correlation Y x Ycv','FontSize',12)
    
    % Graphic of relative error (%)
    subplot(2,3,5.6)
    bar(errOI,'k')
    axis tight
    hline([5,-5],'-r');
    xlabel ('Samples','FontSize',12)
    ylabel ('Relative Error (%)','FontSize',12)
    title('Relative Error CV','FontSize',12)
end

end

function PLSmodel = cval_pls(X,y,nlv,prep,cval,splits)

% CVAL_PLS .....
%
% Input:
%    X (m,n): Independent variable ordered according to the vector(output of ordx).
%    y (m,p): Dependent variable.
%    nlv: Maximum number of Latent Variables for the PLS model
%    prep: A cell containing standard preprocessing structures for the
%          X- and Y- blocks. { 'X preprocessing' 'y preprocessing' }
%          Default: {'mean' 'mean'}
%          Preprocessing options:
%           'mean': mean center;
%           'auto': autoscale;
%           'none': no preprocessing.
%    cval: [ 'full' | 'cblocks' | {'vblinds'} | 'random' ]
%          A cell containing standard cross-validation (CV).
%          cross-validation options:
%          'full': leave-one-out;
%          'cblocks': contiguous block;
%          'vblinds': venetian blinds;
%          'rnd': random subset).
%      splits: {10} number of subsets to divide data into for CV;
%
%  Output:
%     PLSmodel: structured array containing model and cvalidation information.
%
%  I/O:   PLSmodel = cval_pls(X,y,nlv,prep,cval,splits);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%

[n,~] = size(X);
[~,p] = size(y);

% if strcmpi(cval,'full')
%     if nargin==6 && ~isempty(splits)
%         disp('Splits are not necessary when applied full cross validation')
%     end
% end

%%
if nargin==6 && max(size(splits))==1
    no_sampl=fix(n/splits);
    left_over_samples=mod(n,splits);
end

%%
Ypred = zeros(n,p,nlv);
count=1;

if strcmpi(cval,'full')
    cval='vblinds';
    splits=n;
end

if strcmpi(cval,'random')
    ix=rand(n,1);
    [~,ix]=sort(ix);
end

switch cval
    case {'cblocks','vblinds','random'}
        %%
        for i=1:splits
            if strcmpi(cval,'cblocks')
                if left_over_samples==0
                    p_cvs=((i-1)*no_sampl+1+(count-1):i*no_sampl+(count-1))';
                else
                    p_cvs=((i-1)*no_sampl+1+(count-1):i*no_sampl+count)';
                    count=count+1;
                    left_over_samples=left_over_samples-1;
                end
            elseif strcmpi(cval,'vblinds')
                p_cvs=(i:splits:n)';
            elseif strcmpi(cval,'random')
                p_cvs=(i:splits:n)';
                p_cvs=ix(p_cvs)';
            end
            
            tot=(1:n)';
            tot(p_cvs)=[];
            m_cvs = tot;
            PLSmodel.cv{i}=p_cvs;
            Xseg=X(m_cvs,:);
            Yseg=y(m_cvs,:);
            Xpseg=X(p_cvs,:);
            
            if strcmpi(prep{1},'mean')
                [Xseg,mx]=meancenter(Xseg);
                Xpseg=scalenew(Xpseg,mx);
                Xtr=meancenter(X);
            elseif strcmpi(prep{1},'auto')
                [Xseg,mx,stdx]=autoescale(Xseg);
                Xpseg=scalenew(Xpseg,mx,stdx);
                Xtr=autoescale(X);
            elseif strcmpi(prep{1},'none')
                Xtr = X;
            end
            
            if strcmpi(prep{2},'mean')
                [Yseg,my]=meancenter(Yseg);
                Ytr=meancenter(y);
            elseif strcmpi(prep{2},'auto')
                [Yseg,my,stdy]=autoescale(Yseg);
                Ytr=autoescale(y);
            elseif strcmpi(prep{2},'none')
                Ytr = y;
            end
            
            if p == 1
                B = plsbdg(Xseg,Yseg,nlv);
                [~,~,~,~,VarX] = plsbdg(Xtr,Ytr,nlv);
                for k = 1:nlv
                    Ypred(p_cvs,1,k) = Xpseg*B(:,k);
                end
            else
                [~,~,~,~,~,~,VarX,~] = simplsorig(Xtr,Ytr,nlv);
                for k = 1:nlv
                    B(:,:,k) = simplsorig(Xseg,Yseg,k);
                    Ypred(p_cvs,:,k) = Xpseg*B(:,:,k);
                end
            end
            
            
            if strcmpi(prep{2},'mean')
                for j=1:nlv
                    Ypred(p_cvs,:,j)=scaleback(Ypred(p_cvs,:,j),my); % Subfunction
                end
            elseif strcmpi(prep{2},'auto')
                for j=1:nlv
                    Ypred(p_cvs,:,j)=scaleback(Ypred(p_cvs,:,j),my,stdy); % Subfunction
                end
                
            end
        end
        %%
        RMSECV = zeros(1,nlv);
        RCV = zeros(1,nlv);
        Bias = zeros(1,nlv);
        for i=1:nlv
            [RMSECV(i),RCV(i),Bias(i)]=par(y,Ypred(:,:,i));
        end
end

PLSmodel.B = B;
PLSmodel.Ypred = Ypred;
PLSmodel.RMSECV = RMSECV;
PLSmodel.RCV = RCV;
PLSmodel.Bias =  Bias;
PLSmodel.VarX = VarX;

end

function NhOPS = calchops(X,y,nlv,vector,options)

% Kernel Function OPS
%
% I/O: NhOPS = calchops(X,y,nlv,vector,options);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%

fn = 20;
in=nlv;
[~,n] = size(X);
if fn > n
    fn = n;
end
lvs = in:fn;
steps = length(lvs);

switch vector
    case {'reg','sqr','vip','nas','wght','REG','SQR','VIP','NAS','WGHT'}
        h = waitbarn(0,'Calculating hOPS...');
        Rp = zeros(steps,1);
        for i = 1:steps
            waitbarn(i/steps,h);
            resultvec = vecops(X,y,nlv,lvs(i),vector,options);
            res = resultvec.sel(4);
            Rp(i,:) = res;
        end
        close(h)
        nin = in:fn;
        [~,b] = min(Rp);
        NhOPS = nin(b);
        
    case{'main','MAIN','inter','INTER','all','ALL'}
        h = waitbarn(0,'Calculating hOPS...');
        Rpreg = zeros(steps,1);
        Rpsqr = zeros(steps,1);
        Rpvip = zeros(steps,1);
        Rpnas = zeros(steps,1);
        Rpwght = zeros(steps,1);
        for i = 1:steps
            waitbarn(i/steps,h);
            resultreg = vecops(X,y,nlv,lvs(i),'reg',options);
            resreg = resultreg.sel(4);
            Rpreg(i,:) = resreg;
            resultsqr = vecops(X,y,nlv,lvs(i),'sqr',options);
            ressqr = resultsqr.sel(4);
            Rpsqr(i,:) = ressqr;
            resultvip = vecops(X,y,nlv,lvs(i),'vip',options);
            resvip = resultvip.sel(4);
            Rpvip(i,:) = resvip;
            resultnas = vecops(X,y,nlv,lvs(i),'nas',options);
            resnas = resultnas.sel(4);
            Rpnas(i,:) = resnas;
            resultwght = vecops(X,y,nlv,lvs(i),'wght',options);
            reswght = resultwght.sel(4);
            Rpwght(i,:) = reswght;
        end
        close all force
        [~,b1] = min(Rpreg);
        [~,b2] = min(Rpsqr);
        [~,b3] = min(Rpvip);
        [~,b4] = min(Rpnas);
        [~,b5] = min(Rpwght);
        
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [nin(b1) nin(b2) nin(b3) nin(b4) nin(b5)];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
end

end

function NhOPS = calchopscomb(X,y,nlv,vector,options)

% Kernel Function OPS
%
% I/O: NhOPS = calchops(X,y,nlv,vector,options);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%
fn = 20;
in=nlv;
[~,n] = size(X);
if fn > n
    fn = n;
end
lvs = in:fn;
steps = length(lvs);

switch vector
    case {'regcor','regurxy','regcov'}
        h = waitbarn(0,'Calculating hOPS...');
        Rp = zeros(steps,1);
        for i = 1:steps
            waitbarn(i/steps,h);
            resultvec = vecops(X,y,nlv,lvs(i),'reg',options);
            res = resultvec.sel(4);
            Rp(i,:) = res;
        end
        close(h)
        nin = in:fn;
        [~,b] = min(Rp);
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        NhOPS.vec = vec;
        NhOPS.hops = [nin(b) 0 0 0 0];
    case {'corsqr','sqrurxy','sqrcov'}
        h = waitbarn(0,'Calculating hOPS...');
        Rp = zeros(steps,1);
        for i = 1:steps
            waitbarn(i/steps,h);
            resultvec = vecops(X,y,nlv,lvs(i),'sqr',options);
            res = resultvec.sel(4);
            Rp(i,:) = res;
        end
        close(h)
        nin = in:fn;
        [~,b] = min(Rp);
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        NhOPS.vec = vec;
        NhOPS.hops = [0 nin(b) 0 0 0];
    case {'corvip','vipurxy','vipcov'}
        h = waitbarn(0,'Calculating hOPS...');
        Rp = zeros(steps,1);
        for i = 1:steps
            waitbarn(i/steps,h);
            resultvec = vecops(X,y,nlv,lvs(i),'vip',options);
            res = resultvec.sel(4);
            Rp(i,:) = res;
        end
        close(h)
        nin = in:fn;
        [~,b] = min(Rp);
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        NhOPS.vec = vec;
        NhOPS.hops = [0 0 nin(b) 0 0];
    case {'cornas','nasurxy','nascov'}
        h = waitbarn(0,'Calculating hOPS...');
        Rp = zeros(steps,1);
        for i = 1:steps
            waitbarn(i/steps,h);
            resultvec = vecops(X,y,nlv,lvs(i),'nas',options);
            res = resultvec.sel(4);
            Rp(i,:) = res;
        end
        close(h)
        nin = in:fn;
        [~,b] = min(Rp);
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        NhOPS.vec = vec;
        NhOPS.hops = [0 0 0 nin(b) 0];
    case {'corwght','urxywght','wghtcov'}
        h = waitbarn(0,'Calculating hOPS...');
        Rp = zeros(steps,1);
        for i = 1:steps
            waitbarn(i/steps,h);
            resultvec = vecops(X,y,nlv,lvs(i),'wght',options);
            res = resultvec.sel(4);
            Rp(i,:) = res;
        end
        close(h)
        nin = in:fn;
        [~,b] = min(Rp);
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        NhOPS.vec = vec;
        NhOPS.hops = [0 0 0 0 nin(b)];        
    case{'prodall'}
        h = waitbarn(0,'Calculating hOPS...');
        Rpreg = zeros(steps,1);
        Rpsqr = zeros(steps,1);
        Rpvip = zeros(steps,1);
        Rpnas = zeros(steps,1);
        Rpwght = zeros(steps,1);
        for i = 1:steps
            waitbarn(i/steps,h);
            resultreg = vecops(X,y,nlv,lvs(i),'reg',options);
            resreg = resultreg.sel(4);
            Rpreg(i,:) = resreg;
            resultsqr = vecops(X,y,nlv,lvs(i),'sqr',options);
            ressqr = resultsqr.sel(4);
            Rpsqr(i,:) = ressqr;
            resultvip = vecops(X,y,nlv,lvs(i),'vip',options);
            resvip = resultvip.sel(4);
            Rpvip(i,:) = resvip;
            resultnas = vecops(X,y,nlv,lvs(i),'nas',options);
            resnas = resultnas.sel(4);
            Rpnas(i,:) = resnas;
            resultwght = vecops(X,y,nlv,lvs(i),'wght',options);
            reswght = resultwght.sel(4);
            Rpwght(i,:) = reswght;
        end
        close all force
        [~,b1] = min(Rpreg);
        [~,b2] = min(Rpsqr);
        [~,b3] = min(Rpvip);
        [~,b4] = min(Rpnas);
        [~,b5] = min(Rpwght);
        
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [nin(b1) nin(b2) nin(b3) nin(b4) nin(b5)];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
    case {'corurxy','corcov','urxycov'}
        NhOPS.hops = [0 0 0 0 0];
    case {'regsqr'}
        h = waitbarn(0,'Calculating hOPS...');
        Rpreg = zeros(steps,1);
        Rpsqr = zeros(steps,1);
        for i = 1:steps
            waitbarn(i/steps,h);
            resultreg = vecops(X,y,nlv,lvs(i),'reg',options);
            resreg = resultreg.sel(4);
            Rpreg(i,:) = resreg;
            resultsqr = vecops(X,y,nlv,lvs(i),'sqr',options);
            ressqr = resultsqr.sel(4);
            Rpsqr(i,:) = ressqr;
        end
        close all force
        [~,b1] = min(Rpreg);
        [~,b2] = min(Rpsqr);
        
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [nin(b1) nin(b2) 0 0 0];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
    case {'regvip'}
        h = waitbarn(0,'Calculating hOPS...');
        Rpreg = zeros(steps,1);
        Rpvip = zeros(steps,1);
        for i = 1:steps
            waitbarn(i/steps,h);
            resultreg = vecops(X,y,nlv,lvs(i),'reg',options);
            resreg = resultreg.sel(4);
            Rpreg(i,:) = resreg;
            resultvip = vecops(X,y,nlv,lvs(i),'vip',options);
            resvip = resultvip.sel(4);
            Rpvip(i,:) = resvip;
            
        end
        close all force
        [~,b1] = min(Rpreg);
        [~,b3] = min(Rpvip);
        
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [nin(b1) 0 nin(b3) 0 0];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
    case {'regnas'}
        h = waitbarn(0,'Calculating hOPS...');
        Rpreg = zeros(steps,1);
        Rpnas = zeros(steps,1);
        for i = 1:steps
            waitbarn(i/steps,h);
            resultreg = vecops(X,y,nlv,lvs(i),'reg',options);
            resreg = resultreg.sel(4);
            Rpreg(i,:) = resreg;
            resultnas = vecops(X,y,nlv,lvs(i),'nas',options);
            resnas = resultnas.sel(4);
            Rpnas(i,:) = resnas;
        end
        close all force
        [~,b1] = min(Rpreg);
        [~,b4] = min(Rpnas);
              
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [nin(b1) 0 0 nin(b4) 0];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
    case {'regwght'}
        h = waitbarn(0,'Calculating hOPS...');
        Rpreg = zeros(steps,1);
        Rpwght = zeros(steps,1);
        for i = 1:steps
            waitbarn(i/steps,h);
            resultreg = vecops(X,y,nlv,lvs(i),'reg',options);
            resreg = resultreg.sel(4);
            Rpreg(i,:) = resreg;
            resultwght = vecops(X,y,nlv,lvs(i),'wght',options);
            reswght = resultwght.sel(4);
            Rpwght(i,:) = reswght;
        end
        close all force
        [~,b1] = min(Rpreg);
        [~,b5] = min(Rpwght);
        
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [nin(b1) 0 0 0 nin(b5)];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
    case {'sqrvip'}
        h = waitbarn(0,'Calculating hOPS...');
        Rpsqr = zeros(steps,1);
        Rpvip = zeros(steps,1);
        for i = 1:steps
            waitbarn(i/steps,h);
            resultsqr = vecops(X,y,nlv,lvs(i),'sqr',options);
            ressqr = resultsqr.sel(4);
            Rpsqr(i,:) = ressqr;
            resultvip = vecops(X,y,nlv,lvs(i),'vip',options);
            resvip = resultvip.sel(4);
            Rpvip(i,:) = resvip;
        end
        close all force
        [~,b2] = min(Rpsqr);
        [~,b3] = min(Rpvip);
              
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [0 nin(b2) nin(b3) 0 0];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
    case{'sqrnas'}
        h = waitbarn(0,'Calculating hOPS...');
        Rpsqr = zeros(steps,1);
        Rpnas = zeros(steps,1);
        for i = 1:steps
            waitbarn(i/steps,h);
            resultsqr = vecops(X,y,nlv,lvs(i),'sqr',options);
            ressqr = resultsqr.sel(4);
            Rpsqr(i,:) = ressqr;
            resultnas = vecops(X,y,nlv,lvs(i),'nas',options);
            resnas = resultnas.sel(4);
            Rpnas(i,:) = resnas;
        end
        close all force
        [~,b2] = min(Rpsqr);
        [~,b4] = min(Rpnas);
        
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [0 nin(b2) 0 nin(b4) 0];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
    case {'sqrwght'}
        h = waitbarn(0,'Calculating hOPS...');
        Rpsqr = zeros(steps,1);
        Rpwght = zeros(steps,1);
        for i = 1:steps
            waitbarn(i/steps,h);
            resultsqr = vecops(X,y,nlv,lvs(i),'sqr',options);
            ressqr = resultsqr.sel(4);
            Rpsqr(i,:) = ressqr;
            resultwght = vecops(X,y,nlv,lvs(i),'wght',options);
            reswght = resultwght.sel(4);
            Rpwght(i,:) = reswght;
        end
        close all force
        [~,b2] = min(Rpsqr);
        [~,b5] = min(Rpwght);
        
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [0 nin(b2) 0 0 nin(b5)];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
    case {'vipnas'}
        h = waitbarn(0,'Calculating hOPS...');
        Rpvip = zeros(steps,1);
        Rpnas = zeros(steps,1);
        for i = 1:steps
            waitbarn(i/steps,h);
            resultvip = vecops(X,y,nlv,lvs(i),'vip',options);
            resvip = resultvip.sel(4);
            Rpvip(i,:) = resvip;
            resultnas = vecops(X,y,nlv,lvs(i),'nas',options);
            resnas = resultnas.sel(4);
            Rpnas(i,:) = resnas;
        end
        close all force
        [~,b3] = min(Rpvip);
        [~,b4] = min(Rpnas);
        
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [0 0  nin(b3) nin(b4) 0];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
    case {'vipwght'}
        h = waitbarn(0,'Calculating hOPS...');
        Rpvip = zeros(steps,1);
        Rpwght = zeros(steps,1);
        for i = 1:steps
            waitbarn(i/steps,h);
            resultvip = vecops(X,y,nlv,lvs(i),'vip',options);
            resvip = resultvip.sel(4);
            Rpvip(i,:) = resvip;
            resultwght = vecops(X,y,nlv,lvs(i),'wght',options);
            reswght = resultwght.sel(4);
            Rpwght(i,:) = reswght;
        end
        close all force
        [~,b3] = min(Rpvip);
        [~,b5] = min(Rpwght);
        
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [0 0 nin(b3) 0 nin(b5)];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
    case {'naswght'}
        h = waitbarn(0,'Calculating hOPS...');
        Rpnas = zeros(steps,1);
        Rpwght = zeros(steps,1);
        for i = 1:steps
            waitbarn(i/steps,h);
            resultnas = vecops(X,y,nlv,lvs(i),'nas',options);
            resnas = resultnas.sel(4);
            Rpnas(i,:) = resnas;
            resultwght = vecops(X,y,nlv,lvs(i),'wght',options);
            reswght = resultwght.sel(4);
            Rpwght(i,:) = reswght;
        end
        close all force
        [~,b4] = min(Rpnas);
        [~,b5] = min(Rpwght);
        
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [0 0 0 nin(b4) nin(b5)];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
end

end

function [T,W,b,Q,P,b0] = nipops(X,y,nlv)

% Kernel Function OPS: NIPALS algorithm!
%
%  Input:
%    X (m,n): Independent variable.
%    y (m,p): Dependent variable.
%        nlv: Number of latent variables of model.
%
%   Output:
%     T (m,nlv): Matrix of X scores.
%     P (n,nlv): Matrix of X loadings.
%     Q (p,nlv): Matrix of y scores.
%     W (n,nlv): Matrix of X weights.
%     B (n,nlv): Matrix with regression vector of each latent variable.
%
%   I/O: [T,W,b,Q,P,b0] = nipops(X,y,nlv);
%

%%
E = X;
F = y;
u = y;
xm = mean(X);
ym = mean(y);

for i = 1:nlv
    w = u'*E/(u'*u);
    w = w./norm(w);
    w = w';
    t = (E*w)./(w'*w);
    q = (t'*F)./(t'*t);
    q = q';
    p = (t'*E)./(t'*t);
    p = p';
    
    T(:,i) = t;
    Q(:,i) = q;
    W(:,i) = w;
    P(:,i) = p;
    
    E = E - t*p';
    F = F - t*q';
    
    b(:,i) = W*((P'*W)\Q');
    b0(i) = ym - xm*b(:,i);
end

end

function [Xor,ind] = ordx(X,vec)

% Kernel Function OPS
%
% I/O: [Xor,ind] = ordx(X,vec);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%

vec = reshape(vec,length(vec),1);
[~,i] = sort(vec,1,'descend');
ind = i';
Xor = X(:,ind);

end

function [sel,resultops] = perfops(Xor,y,ind,options,hmod)

% Kernel Function OPS
%
% I/O: [sel,resultops] = perfops(Xor,y,ind,options,nlv);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%
win = options.window;
incr = options.increment;
p = options.percentage;

[~,n] = size(Xor);
if n <= 2
    
    prep = options.preprocessing;
    cval = options.crossvalidation;
    splits = options.split;

    if n < hmod
        hmod = n;
    end
    
    PLSmodel = cval_pls(Xor,y,hmod,prep,cval,splits);
    rmsecvout = PLSmodel.RMSECV;
    rcvout = PLSmodel.RCV;
    err = rmsecvout;
    errn = err/norm(err,inf);
    q = rcvout;
    param = q./errn;
    [~,g] = sort(param);
    lv_opt1 = g(length(g));
    lv = lv_opt1;
    rmsecv_b = err(lv);
    rcv_b = q(lv);
    
    sel = [1 n lv rmsecv_b rcv_b];
    var_sel = {1:n};
    nlvs = num2str(lv);
    nlv_b = lv;
    nseta = 1;
    nsetb= 1;
else
    
    if win >= n
        win = n-1;
        incr = 1;
    end
    
    if p < 10
        nv = round((10*n)/100);
    elseif p > 100
        nv = n;
    else
        nv = round((p*n)/100);
    end
    
    prep = options.preprocessing;
    cval = options.crossvalidation;
    splits = options.split;
    
    X1 = Xor(:,1:win);
    [~,b] = size(X1);
    
    if b < hmod
        win = hmod;
        X1 = Xor(:,1:win);
    end
    
    PLSmodel = cval_pls(X1,y,hmod,prep,cval,splits);
    rmsecvout(1) = {PLSmodel.RMSECV};
    rcvout(1) = {PLSmodel.RCV};
    
    var_sel(1)={ind(1:win)};
    
    i = 2;
    ex = win + incr;
    
    while ex <= nv
        X1 = Xor(:,1:ex);

        PLSmodel = cval_pls(X1,y,hmod,prep,cval,splits);
        rmsecvout(i) = {PLSmodel.RMSECV};
        rcvout(i) = {PLSmodel.RCV};
        
        var_sel(i)={ind(1:ex)};
        ex = ex + incr;
        
        i = i+1;
    end
    
    if (ex-incr) < nv
        
        PLSmodel = cval_pls(X1,y,hmod,prep,cval,splits);
        rmsecvout(i) = {PLSmodel.RMSECV};
        rcvout(i) = {PLSmodel.RCV};
        
        var_sel(i)={ind(1:nv)};
        i = i+1;
    end
    
    % Full set
    if p < 100
        
        PLSmodel = cval_pls(X1,y,hmod,prep,cval,splits);
        rmsecvout(i) = {PLSmodel.RMSECV};
        rcvout(i) = {PLSmodel.RCV};
        
        var_sel(i)={ind(1:n)};
    end
    
    [~,h] = size(rmsecvout);
    s = 1:h;
    rmsecv_b = zeros(h,1);
    rcv_b = zeros(h,1);
    nlv_b = zeros(h,1);
    nvars = zeros(h,1);
    
    for j = 1:h
        err = rmsecvout{j};
        errn = err/norm(err,inf);
        q = rcvout{j};
        vars = size(var_sel{j},2);
        param = q./errn;
        [~,g] = sort(param);
        lv_opt1 = g(length(g));
        lv = lv_opt1;
        rmsecv_b(j) = err(lv);
        rcv_b(j) = q(lv);
        nlv_b(j) = lv;
        nvars(j) = vars;
    end
    
    resultab = [s' nvars nlv_b rmsecv_b rcv_b];
    T = resultab;
%     assignin('base','T',T)
    de1 = isnan(T);
    del1 = find(de1(:,4)==1);
    if ~isempty(del1)
        T(del1,:) = [];
        de2 = isnan(T);
        del2 = find(de2(:,4)==1);
        if ~isempty(del2)
            T(del2,:) = [];
        end
    end
    pn = ((T(:,5)*norm(T(:,4),inf))./T(:,4));
    pm = max(pn);
    L = find(pn == pm);
    sel = [T(L,1) T(L,2) T(L,3) T(L,4) T(L,5)];
    
    if size(sel,1) > 1
        r = find(min(sel(:,3)));
        sel = sel(r(1),:);
    end
    nlvs = num2str(nlv_b);
    nseta = length(nvars);
    nsetb = (1:nseta);
end


%%
% Saving results
resultops.setvar =  var_sel;
graph.nsetb = nsetb;
graph.rmsecv_b = rmsecv_b;
graph.rcv_b = rcv_b;
graph.nlvs = nlvs;
graph.nlv_b = nlv_b;
resultops.graph = graph;

end

function [B,t,s,p,varX] = plsbdg(x,y,n)

% Kernel Function OPS: PLSBDG Partial Least Squares using bidiagonal algorithm for OPS. 
%
%  Input:
%    X (m,n): Independent variable.
%    y (m,p): Dependent variable.
%        nlv: Number of latent variables of model.
%
%   Output:
%      B(n,nlv): Matrix with regression vector of each latent variable.
%	   t(m,nlv): Matrix of X scores.
%	 s(nlv,nlv): Diagonal matrix of singular values.
%      p(n,nlv): Matrix of X loadings.
%   varX(1,nlv): Explained variance of X.
%
%   I/O: [B,t,s,p,varX] = plsbdg(X,y,nlv);
%
%	Copyright © OPS_Toolbox 2.0.

[~,nc] = size(x);
if n == 1
    n = 2;
end

if n > nc
    n = nc;
end    

[t,w,rv1,p] = bdiag(x,y,n);
s = diag(w) + diag(rv1,1);
B = zeros(nc,n);
for i = 1:n
    U = t(:,1:i);
    S = s(1:i,1:i);
    V = p(:,1:i); % A matriz de loadings não está transposta
    c = U'*y;
    q = S\c;
    B(:,i) = V*q; % Matriz com vetores de regressão
end

T = t*s;
varX = cumsum(sum(T.^2)*100)/sum(sum(T.^2));
end

% Extra functions
function [u,w,rv1,v] = bdiag(a,y,n)

% BIDIAGONALIZATION of matrix using Lanczos algorithm
% [u,w,rv1,v] = bdiag(a,n,y);
%
% Inputs:
%		a:	Rectangular matrix
% 		n:	Numbers of latent variables
%  Optional Input
%		<y>:	Initial estimate of y (column space), OR 
%				y vector for PLS regression,  OR
%				if y omitted,  y = [1 1 .. 1]'
% Outputs:
%		u:	row eigenvectors
% 		w:	diagonal elements of the singular value matrix s = t' * a * p
%		rv1:	off-diagonal elements
%		v:	spectral eigenvectors

[nrow,~] = size(a);			%dimensions of matrix
if (nargin <= 2), 
	y = ones(nrow,1); 
end					%if y not supplied, create it

v = []; %initialization 

for i=1:n
	vi = normalize(a'*y);		%initial estimate of v vector
   	ui = a * vi;			%initial estimate of u vector
 
	if i==1,
		[ui,nor] = normalize(ui);	%from normalization of u...
        w(i) = nor;			% save factor as the weights
	else
		temp = u' * ui;		%temp off-diagonal, store as...
		rv1(i-1) = temp(i-1);		% the off-diagonal elements
		[ui,nor] = normalize(ui - u * temp);	%from normalization of u...
        w(i) = nor;			% save factor as the weights
	end

	v = [v vi];			%store the computed v vectors
	u(:,i) = ui;	
	if i<n,
		% and update y, subtracting influence of u vector
		y = y - ui * (ui' * y);		
	end	%if
end		%for i
end

function [vn,fn] = normalize(v)
[~,n] = size(v);

for i = 1:n
    fn(i) = norm(v(:,i));
    vn(:,i) = v(:,i)/fn;
end

end

function resultvec = vecops(X,y,hmod,hops,vector,options)

% Kernel Function OPS
%
% I/O: resultvec = vecops(X,y,hmod,hops,vector,options);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%
switch vector
    case {'reg','REG'}
        % Regression (REG) vector
        if length(hops)>1
            reg = vecreg(X,y,hops(1));
        else
            reg = vecreg(X,y,hops);
        end
        [Xor,ind] = ordx(X,reg);
        [sel,resultreg] = perfops(Xor,y,ind,options,hmod);
        resultvec.vect = reg;
        resultvec.graph = resultreg.graph;
        resultvec.vec = char('reg');
        resultvec.varsel = resultreg.setvar;
        resultvec.sel = sel;
        resultvec.namesel =  {'Subset' 'nVars' 'nlv' 'RMSECV' 'Rcv'};
    case {'cor','COR'}
        % Correlation (COR) vector
        cor = veccorr(X,y);
        [Xor,ind] = ordx(X,cor);
        [sel,resultcor] = perfops(Xor,y,ind,options,hmod);
        resultvec.vect = cor;
        resultvec.graph = resultcor.graph;
        resultvec.vec = char('cor');
        resultvec.varsel = resultcor.setvar;
        resultvec.sel = sel;
        resultvec.namesel =  {'Subset' 'nVars' 'nlv' 'RMSECV' 'Rcv'};
    case {'sqr','SQR'}
        % Residual (SQR) vector
        if length(hops)>1
            sqr = vecsqr(X,y,hops(2));
        else
            sqr = vecsqr(X,y,hops);
        end
        [Xor,ind] = ordx(X,sqr);
        [sel,resultsqr] = perfops(Xor,y,ind,options,hmod);
        resultvec.vect = sqr;
        resultvec.graph = resultsqr.graph;
        resultvec.vec = char('sqr');
        resultvec.varsel = resultsqr.setvar;
        resultvec.sel = sel;
        resultvec.namesel =  {'Subset' 'nVars' 'nlv' 'RMSECV' 'Rcv'};
    case {'vip','VIP'}
        % Variable influence on projection (VIP) vector
        if length(hops)>1
            vip = vecvip(X,y,hops(3));
        else
            vip = vecvip(X,y,hops);
        end
        [Xor,ind] = ordx(X,vip);
        [sel,resultvip] = perfops(Xor,y,ind,options,hmod);
        resultvec.vect = vip;
        resultvec.graph = resultvip.graph;
        resultvec.vec = char('vip');
        resultvec.varsel = resultvip.setvar;
        resultvec.sel = sel;
        resultvec.namesel =  {'Subset' 'nVars' 'nlv' 'RMSECV' 'Rcv'};
    case {'nas','NAS'}
        % Net analyte signal (NAS) vector
        if length(hops)>1
            nas = vecnas(X,y,hops(4));
        else
            nas = vecnas(X,y,hops);
        end
        [Xor,ind] = ordx(X,nas);
        [sel,resultnas] = perfops(Xor,y,ind,options,hmod);
        resultvec.vect = nas;
        resultvec.graph = resultnas.graph;
        resultvec.vec = char('nas');
        resultvec.varsel = resultnas.setvar;
        resultvec.sel = sel;
        resultvec.namesel =  {'Subset' 'nVars' 'nlv' 'RMSECV' 'Rcv'};
    case {'urxy','URXY'}
        % Error in y (urxy) vector
        urxy = vecerry (X,y);
        [Xor,ind] = ordx(X,urxy);
        [sel,resulturxy] = perfops(Xor,y,ind,options,hmod);
        resultvec.vect = urxy;
        resultvec.graph = resulturxy.graph;
        resultvec.vec = char('urxy');
        resultvec.varsel = resulturxy.setvar;
        resultvec.sel = sel;
        resultvec.namesel =  {'Subset' 'nVars' 'nlv' 'RMSECV' 'Rcv'};
    case {'wght','WGHT'}
        % Weight (WGHT) vector
        if length(hops)>1
            wght = vecweight (X,y,hops(5));
        else
            wght = vecweight (X,y,hops);
        end
        [Xor,ind] = ordx(X,wght);
        [sel,resultwght] = perfops(Xor,y,ind,options,hmod);
        resultvec.vect = wght;
        resultvec.graph = resultwght.graph;
        resultvec.vec = char('wght');
        resultvec.varsel = resultwght.setvar;
        resultvec.sel = sel;
        resultvec.namesel =  {'Subset' 'nVars' 'nlv' 'RMSECV' 'Rcv'};
    case {'cov','COV'}
        % Covariance procedures (COV) vector
        cov = veccov (X,y);
        [Xor,ind] = ordx(X,cov);
        [sel,resultcov] = perfops(Xor,y,ind,options,hmod);
        resultvec.vect = cov;
        resultvec.graph = resultcov.graph;
        resultvec.vec = char('cov');
        resultvec.varsel = resultcov.setvar;
        resultvec.sel = sel;
        resultvec.namesel =  {'Subset' 'nVars' 'nlv' 'RMSECV' 'Rcv'};
   case {'regcor','regsqr','regvip','regnas','regurxy','regwght','regcov','corsqr','corvip','cornas','corurxy','corwght','corcov','sqrvip','sqrnas','sqrurxy','sqrwght','sqrcov','vipnas','vipurxy','vipwght','vipcov','nasurxy','naswght','nascov','urxywght','urxycov','wghtcov','prodall'}
       vecopt = {'regcor';'regsqr';'regvip';'regnas';'regurxy';'regwght';'regcov';'corsqr';'corvip';'cornas';'corurxy';'corwght';'corcov';'sqrvip';'sqrnas';'sqrurxy';'sqrwght';'sqrcov';'vipnas';'vipurxy';'vipwght';'vipcov';'nasurxy';'naswght';'nascov';'urxywght';'urxycov';'wghtcov';'prodall'};
       vec = find(strcmp(vecopt,char(vector)));
       veccomb = combvec(X,y,hmod,hops,options,vec);
       resultvec = veccomb;
    case {'main','MAIN'}
        l = length(hops);
        if l < 5
            error('Incorrect length of hOPS. Please provide a vector with 5 elements.')
        end
        vecmain = mainvec(X,y,hmod,hops,options);
        resultvec = vecmain;
    case {'inter','INTER'}
        l = length(hops);
        if l < 5
            error('Incorrect length of hOPS. Please provide a vector with 5 elements.')
        end
        vecinter = intervec(X,y,hmod,hops,options);
        resultvec = vecinter;
    case {'all','ALL'}
        l = length(hops);
        if l < 5
            error('Incorrect length of hOPS. Please provide a vector with 5 elements.')
        end
        vecall = allvec(X,y,hmod,hops,options);
        resultvec = vecall;
    case {'custom','Custom','CUSTOM'}
        % Custom vector
        custom = options.custom;
        [Xor,ind] = ordx(X,custom);
        [sel,resultcustom] = perfops(Xor,y,ind,options,hmod);
        resultvec.vect = custom;
        resultvec.graph = resultcustom.graph;
        resultvec.vec = char('custom');
        resultvec.varsel = resultcustom.setvar;
        resultvec.sel = sel;
        resultvec.namesel =  {'Subset' 'nVars' 'nlv' 'RMSECV' 'Rcv'};
    otherwise
        display('Invalid option: "main" option will be used.');
        l = length(hops);
        if l < 5
            error('Incorrect length of hOPS. Please provide a vector with 5 elements.')
        end
        vecmain = mainvec(X,y,hmod,hops,options);
        resultvec = vecmain;
end

end

function Auto = autoOPS(X,y,hmod,hops,vector,options)

% AUTOOPS is an automatic method to perform the variable selection using
% several informative vectors and their combinations, providing the best
% result.
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%% VARIABLE SELECTION

Auto.name = 'autoOPS';
Auto.options = options;

switch vector
    case {'reg','cor','sqr','vip','nas','urxy','wght','cov','REG','COR','SQR','VIP','NAS','URXY','WGHT','COV'}
        resultvec = vecops(X,y,hmod,hops,vector,options);
        Auto.result.graph = resultvec.graph;
        Auto.result.vect = resultvec.vect';
        vect = resultvec.vect;
        set = resultvec.sel(1);
        varsels = resultvec.varsel{set};
        sel = resultvec.sel;
        b = {vector};
    case {'regcor','regsqr','regvip','regnas','regurxy','regwght','regcov','corsqr','corvip','cornas','corurxy','corwght','corcov','sqrvip','sqrnas','sqrurxy','sqrwght','sqrcov','vipnas','vipurxy','vipwght','vipcov','nasurxy','naswght','nascov','urxywght','urxycov','wghtcov','prodall'}
        resultvec = vecops(X,y,hmod,hops,vector,options);
        Auto.result.graph = resultvec.(char(vector)).graph;
        Auto.result.vect = resultvec.vect;
        vect = resultvec.vect;
        sel = resultvec.sel;
        varsels = resultvec.(char(vector)).setvar{1,sel(1)};
        b = {vector};
    case {'main','MAIN'}
        resultvec = vecops(X,y,hmod,hops,vector,options);
        rnames_sel = resultvec.vectors;
        colnames_sel = {'Subset' 'nVars' 'hMod' 'RMSECV' 'Rcv'};
        sel = array2table(resultvec.sel,'RowNames',rnames_sel,'VariableNames',colnames_sel);
        Auto.result.sel = sel;
        prod = (resultvec.sel(:,4))./resultvec.sel(:,5);
        [~,i] = sort(prod,1);
        vec = i(1);
        vect = resultvec.Xvet(i(1),:);
        set = resultvec.sel(vec);
        b = resultvec.vectors(vec);
        [~,varsels,sel] = setmain(vec,set,resultvec);
        Auto.result.graph = resultvec.(char(b)).graph;
        Auto.result.vect = vect;
    case {'inter','INTER'}
        resultvec = vecops(X,y,hmod,hops,vector,options);
        rnames_sel = resultvec.vectors;
        colnames_sel = {'Subset' 'nVars' 'hMod' 'RMSECV' 'Rcv'};
        sel = array2table(resultvec.sel,'RowNames',rnames_sel,'VariableNames',colnames_sel);
        Auto.result.sel = sel;
        prod = (resultvec.sel(:,4))./resultvec.sel(:,5);
        [~,i] = sort(prod,1);
        vec = i(1);
        vect = resultvec.Xvet(i(1),:);
        set = resultvec.sel(vec);
        b = resultvec.vectors(vec);
        [~,varsels,sel] = setinter(vec,set,resultvec);
        Auto.result.graph = resultvec.(char(b)).graph;
        Auto.result.vect = vect;
    case {'all','ALL'}
        resultvec = vecops(X,y,hmod,hops,vector,options);
        rnames_sel = resultvec.vectors;
        colnames_sel = {'Subset' 'nVars' 'hMod' 'RMSECV' 'Rcv'};
        sel = array2table(resultvec.sel,'RowNames',rnames_sel,'VariableNames',colnames_sel);
        Auto.result.sel = sel;
        prod = (resultvec.sel(:,4))./resultvec.sel(:,5);
        [~,i] = sort(prod,1);
        vec = i(1);
        vect = resultvec.Xvet(i(1),:);
        set = resultvec.sel(vec);
        b = resultvec.vectors(vec);
        [~,varsels,sel] = setall(vec,set,resultvec);
        Auto.result.graph = resultvec.(char(b)).graph;
        Auto.result.vect = vect;
    case {'custom','Custom','CUSTOM'}
        resultvec = vecops(X,y,hmod,hops,vector,options);
        Auto.result.graph = resultvec.graph;
        Auto.result.vect = resultvec.vect;
        vect = resultvec.vect';
        sel = resultvec.sel;
        varsels = resultvec.varsel{1,sel(1)};
        b = {vector};
end

%% ORGANIZING RESULTS...

Auto.vec_option = char(vector);

cval = options.crossvalidation;
splits = options.split;
prep = options.preprocessing;

hmod = sel(3);
PLSmodelF = cval_pls(X,y,hmod,prep,cval,splits);
Xs = X(:,varsels);
PLSmodelS = cval_pls(Xs,y,hmod,prep,cval,splits);
ycv = PLSmodelS.Ypred(:,:,hmod);
err = (((y-ycv)./y)*100);
res = y-ycv;
nu = (1:length(y))';
T = [nu y ycv res err];

full = [0 size(X,2) hmod   PLSmodelF.RMSECV(hmod) PLSmodelF.RCV(hmod)];
resulttab = [full;sel];
rnames_result = {'Full';'Sel'};
colnames_result = {'Subset' 'nVars' 'hMod' 'RMSECV' 'Rcv'};
fprintf(' - - - - - - - - - - - - autoOPS - - - - - - - - - - - -');
Quality_parameters = array2table(resulttab,'RowNames',rnames_result,'VariableNames',colnames_result)
colnames_MP = {'Sample' 'y' 'ycv' 'residual' 'Rel_error'};
MP_results = array2table(T,'VariableNames',colnames_MP)

%% SAVING RESULTS...

switch vector
    case {'custom','Custom','CUSTOM'}
        Auto.hmod = hmod;
        Auto.best_vec = 'custom';
    case {'main','MAIN','inter','INTER','all','ALL'}
        Auto.best_vec = char(b);
        Auto.hmod = hmod;
        Auto.hops = hops;
    otherwise
        Auto.best_vec = char(vector);
        Auto.hmod = hmod;
        Auto.hops = hops;
end
Auto.vs = varsels;
Auto.informative_vector = vect(varsels);
Auto.Quality_parameters = Quality_parameters;
Auto.MP_results = MP_results;

%% PLOT OF RESULTS
autofig(X,Auto)

end

function [hmod,rmsecv] = calchmod(X,y,options,~)
%---------------------------------------------------------------------------
%  CALCHMOD algorithm calculates the hmod, i.e., the number of latent
%  variables to build PLS models.
%---------------------------------------------------------------------------
% Input:
%
%    X (m,n): Independent variable ordered according to the vector(output of ordx).
%
%    y (m,1): Dependent variable.
%
%    options: (OPTIONAL)
%      Window: {10}  window of initial variables.
%      Increment: {5} increment of variables to be added over the window.
%      Percentage: {100} percentage of the variables to be considered in
%                  the selection.
%      cv: [ 'full' | 'cblocks' | {'vblinds'} | 'random' ]
%          A cell containing standard cross-validation (CV).
%          cross-validation options:
%          'full': leave-one-out;
%          'cblocks': contiguous block;
%          'vblinds': venetian blinds;
%          'random': random subset).
%      Split: {10} number of subsets to divide data into for CV;
%      Preprocessing: A cell containing standard preprocessing structures
%      for the X- and Y- blocks.
%          { 'X preprocessing' 'y preprocessing' }
%          Default: {'mean' 'mean'}
%          Preprocessing options:
%           'mean': mean center;
%           'auto': autoscale;
%           'none': no preprocessing.
%
%      Options e.g.:
%        options.window=10;
%        options.increment=5;
%        options.percentage=100;
%        options.crossvalidation='random';
%        options.split=10;
%        options.preprocessing={'mean','mean'};
%        options.criteria = [0.02 10];
%        options.var = 50;
%        options.iOPS_type = 'autoOPS';
%        options.calhops = 'no';
%
%   Output:
%
%       hmod: Number of latent variables of model.
%             
%   I/O:  hmod = calchmod(X,y);
%
%         hmod = calchmod(X,y,options);
%
% See also: NLVOPS, NEWOPS, SELOPS, CREATEOPTIONS
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%
if nargin < 3
    options.name = 'options';
    options.crossvalidation = 'vblinds';
    options.split = 10;
    options.preprocessing = {'mean' 'mean'};
    options.fn = 20;
end

prep = options.preprocessing;
cval = options.crossvalidation;
splits = options.split;
fn = 20;

l = size(X,2);
if fn > l
    fn = l;
end

PLSmodel = cval_pls(X,y,fn,prep,cval,splits);
[~,c] = sort(PLSmodel.RMSECV,2,'ascend');
hmod1 = c(1);
if hmod1 == 1
    hmod1 = c(2);
end
d = diff(PLSmodel.VarX);
r = find(d<1);
if isempty(r)
    r = find(d==min(d));
end
hmod2 = r(1);


hmod = min([hmod1 hmod2]);
rmsecv = PLSmodel.RMSECV;

if nargin < 4
figure
subplot(1,2,1)
plot(1:fn,PLSmodel.RMSECV,'-ko');
hold on
plot(hmod1,PLSmodel.RMSECV(hmod1),'rpentagram','MarkerSize',15)
xlabel ('Latent Variables','FontSize',12)
ylabel ('RMSECV','FontSize',12)
title ('hmod','FontSize',12)
grid on
subplot(1,2,2)
plot(1:fn,PLSmodel.VarX,'-ko');
hold on
plot(hmod2,PLSmodel.VarX(hmod2),'rpentagram','MarkerSize',15)
xlabel ('Latent Variables','FontSize',12)
ylabel ('% Explained Variance','FontSize',12)
title ('hmod','FontSize',12)
grid on
end

end

function Feed = feedOPS(X,y,hmod,hops,vector,options)

% FEEDOPS is a method to apply the variable selection using a feedback OPS, 
% wherein pre-selected variables can return to a new selection until
% achieving specific criteria.
%
% Copyright: Jussara V. Roque, 2019.
% Check by JVR: 20/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo. 
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019). 
% https://doi.org/10.1016/j.aca.2019.05.039 

%% FULL MODEL

cval = options.crossvalidation;
splits = options.split;
prep = options.preprocessing;

a = 1;
PLSmodel = cval_pls(X,y,hmod,prep,cval,splits);
e(a) = PLSmodel.RMSECV(hmod);
rcvf = PLSmodel.RCV(hmod);

%% VARIABLE SELECTION

Feed.name = 'feedOPS';
Feed.options = options;

criteria = options.criteria;
p = criteria(1);
nloop = criteria(2);

Feed.vec_option = char(vector);

compare = strcmp(char(vector),{'custom','Custom','CUSTOM'});
comp = find(compare == 1,1);
if ~isempty(comp)
    custom = options.custom;
    resultf = customvec(X,y,hmod,'auto',[p nloop],options,custom);
    result = {};
    result(1,1) = {resultf};
    Xa = resultf.Xs;
    a = 2;
    e(a) = resultf.resultSEL(4);
    select =  resultf.resultSEL;
    SELresult = [];
    SELresult (1,:) = select;
    vs = resultf.vs;
    VS = {};
    VS(1,1) = {vs};
    Vectors = {};
    Vectors(1,1) = {resultf.vec};
    InfVec = {};
    InfVec(1,1) = {resultf.informative_vector};
    fprintf('Number of loops: #%d\n', 1);
    per = 1;
    counter = 2;
    custom = custom(resultf.vs);
    while per > p
        resultauto = customvec(Xa,y,hmod,'auto',[p nloop],options,custom);
        Xa = resultauto.Xs;
        a = a+1;
        result(1,counter) = {resultauto};
        rmsecv = resultauto.resultSEL(4);
        e(a) = rmsecv;
        select =  resultauto.resultSEL;
        SELresult (counter,:) = select;
        vsn = resultauto.vs;
        vso = VS{1,(counter-1)};
        vs = vso(vsn);
        VS(1,counter) = {vs};
        Vectors(1,counter) = {resultauto.vec};
        InfVec(1,counter) = {resultauto.informative_vector};
        ne = length(e);
        de = e(ne) - e(ne-1);
        per = abs(de)/e(ne);
        fprintf('Number of loops: #%d\n', counter);
        custom = custom(vsn);
        if (de >= 0) || ((counter -1) == nloop)
            break
        end
        counter = counter + 1;
    end
else
    resultf = autofeed(X,y,hmod,hops,vector,options);
    result = {};
    result(1,1) = {resultf};
    Xa = resultf.Xs;
    a = 2;
    e(a) = resultf.resultSEL(4);
    select =  resultf.resultSEL;
    SELresult = [];
    SELresult (1,:) = select;
    vs = resultf.varsels;
    VS = {};
    VS(1,1) = {vs};
    Vectors = {};
    Vectors(1,1) = {resultf.vec};
    InfVec = {};
    InfVec(1,1) = {resultf.informative_vector};
    fprintf('Number of loops: #%d\n', 1);
    per = 1;
    counter = 2;
    while per > p
        [~,n] = size(Xa);
        f = hops > n;
        hops(f) = n;
        resultauto = autofeed(Xa,y,hmod,hops,vector,options);
        Xa = resultauto.Xs;
        a = a+1;
        result(1,counter) = {resultauto};
        rmsecv = resultauto.resultSEL(4);
        e(a) = rmsecv;
        select =  resultauto.resultSEL;
        SELresult (counter,:) = select;
        vsn = resultauto.varsels;
        vso = VS{1,(counter-1)};
        vs = vso(vsn);
        VS(1,counter) = {vs};
        Vectors(1,counter) = {resultauto.vec};
        InfVec(1,counter) = {resultauto.informative_vector};
        ne = length(e);
        de = e(ne) - e(ne-1);
        per = abs(de)/e(ne);
        fprintf('Number of loops: #%d\n', counter);
        counter = counter + 1;
        if (de >= 0) || ((counter -1) == nloop)
            break
        end
    end
end

%% SAVING RESULTS ...

[r,~] = size(SELresult);
if (r == nloop) || (r == 1)
    Feed.info = strcat(['The best variables were selected in loop ',num2str(r),'.']);
    vsopt = VS{1,r};
    infvecopt = InfVec{1,r};
    resultopt = result{1,r};
    T = resultopt.tabMP;
    hmod = SELresult(r,3);
    
    full = [0 size(X,2) hmod   e(1) rcvf];
    resulttab = [full;SELresult];
    rnames_result = {'Full'};
    for i = 1:r
        rnames_result{i+1,1} = strcat('Loop',num2str(i));
    end
    colnames_result = {'Subset' 'nVars' 'nlv' 'RMSECV' 'Rcv'};
    fprintf(' - - - - - - - - - - - - feedOPS - - - - - - - - - - - - -');
    Quality_parameters = array2table(resulttab,'RowNames',rnames_result,'VariableNames',colnames_result)
    colnames_MP = {'Sample' 'y' 'ycv' 'residual' 'Rel_error'};
    MP_results = array2table(T,'VariableNames',colnames_MP)
    
    Feed.nlvs = [hmod hops];
    if ~isempty(comp)
        Feed.nlvs = hmod;
    end
    Feed.vs = vsopt;
    switch vector
        case {'main','MAIN','inter','INTER','all','ALL'}
            Feed.best_vec = char(Vectors{1,r-1});
    end
    Feed.informative_vector = infvecopt;
    Feed.Loops.bestloop = r;
    Feed.Loops.informative_vector = InfVec;
    Feed.Loops.varsels = VS;
    Feed.Loops.vectors = Vectors;
    Feed.Quality_parameters = Quality_parameters;
    Feed.MP_results = MP_results;
else
    Feed.info = strcat(['The best variables were selected in loop ',num2str(r-1),'.']);
    vsopt = VS{1,r-1};
    infvecopt = InfVec{1,r-1};
    resultopt = result{1,r-1};
    T = resultopt.tabMP;
    hmod = SELresult(r-1,3);
    
    full = [0 size(X,2) hmod   e(1) rcvf];
    resulttab = [full;SELresult];
    rnames_result = {'Full'};
    for i=1:r
        rnames_result{i+1,1} = strcat('Loop',num2str(i));
    end
    colnames_result = {'Subset' 'nVars' 'nlv' 'RMSECV' 'Rcv'};
    fprintf(' - - - - - - - - - - - - feedOPS - - - - - - - - - - - -');
    Quality_parameters = array2table(resulttab,'RowNames',rnames_result,'VariableNames',colnames_result)
    colnames_MP = {'Sample' 'y' 'ycv' 'Residual' 'Rel_error'};
    MP_results = array2table(T,'VariableNames',colnames_MP)
    
    Feed.hmod = hmod;
    Feed.hops = hops;
    if ~isempty(comp)
        Feed.hmod = hmod;
    end
    Feed.vs = vsopt;
    switch vector
        case {'main','MAIN','inter','INTER','all','ALL'}
            Feed.best_vec = char(Vectors{1,r-1});
    end
    Feed.informative_vector = infvecopt;
    Feed.Loops.bestloop = r-1;
    Feed.Loops.informative_vector = InfVec;
    Feed.Loops.varsels = VS;
    Feed.Loops.vectors = Vectors;
    Feed.Quality_parameters = Quality_parameters;
    Feed.MP_results = MP_results;
end

%% PLOT OF RESULTS
feedfig(X,Feed)

end

function Inter = iOPS(X,y,hmod,hops,vector,options)

% IOPS is a method to perform the variable selection using autoOPS or
% feedOPS in each subdivision of the full array.
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 20/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%% FULL MODEL

cval = options.crossvalidation;
splits = options.split;
prep = options.preprocessing;

PLSmodel = cval_pls(X,y,hmod,prep,cval,splits);
rmsecvF = PLSmodel.RMSECV(hmod);
rcvf = PLSmodel.RCV(hmod);

[~,n] = size(X);
full = [0 n hmod   rmsecvF rcvf];

%% VARIABLE SELECTION

Inter.name = 'iOPS';
Inter.options = options;

itype = options.iOPS_type;
calhops = options.calhops;

criteria = options.criteria;
pe = criteria(1);
nloop = criteria(2);
var = options.var;

if n < 50
    error('Impossible to perform iOPS')
end
if var < 50
    disp('Intervals less than 50 variables are not allowed!')
    var = 50;
end

int = n/var;
ri = floor(int);
nu = ri*var;
dif = n - nu;

switch vector
    case {'custom','Custom','CUSTOM'}
        calhops = 'no';
        custom = options.custom;
        custom = reshape(custom,1,length(custom));
        interv = {};
        interv(1,1) = {custom(:,1:var)};
        for j = 2:ri
            interv(j,1) = {custom(:,((j-1)*var)+1:(j)*var)};
        end
end
Inter.var_interval = var;
Inter.vec_option = char(vector);

if ischar(itype)
    itype = cellstr(itype);
end
s = size(itype,2);

for j=1:s
    irtype = itype(j);
    irtype = char(irtype);
    switch irtype
        case {'feed','Feed','FEED','feedops','feedOPS','FEEDOPS'}
            display('FeedOPS will be apllied in each interval.');
            Inter.feedOPS.criteria = criteria;
            Inter.feedOPS.type = 'feed';

            % Intervals models
            Xint = {};
            Xint(1,1) = {X(:,1:var)};
            for i = 2:ri
                Xint(i,1) = {X(:,((i-1)*var)+1:(i)*var)};
            end
            [p,~] = size(Xint);
            wb = waitbarn(0,'Interval model...');
            if dif>=50
                Xint(p+1,1) = {X(:,(nu+1:end))};
                Xsel = cell(p+1,1);
                custvec = cel(p+1,1);
                varsI = cell(p+1,1);
                rmsecvI = cell(p+1,1);
                for i = 1:p+1
                    waitbarn(i/(p+1),wb);
                    switch calhops
                        case {'y','Y','yes','Yes','YES'}
                            switch vector
                                case {'reg','sqr','vip','nas','wght','REG','SQR','VIP','NAS','WGHT'}
                                    hopsl = chops(Xint{i,1},y,hmod,20,vector,options);
                                    resultfeedI = feedinter(Xint{i,1},y,[hmod hopsl],vector,[pe nloop],options);
                                    Xsel (i,1) = {resultfeedI.Xs};
                                    varsI(i,1) = {resultfeedI.vs};
                                    rmsecvI(i) = resultfeedI.rmsecv;
                                case {'regcor','regsqr','regvip','regnas','regurxy','regwght','regcov','corsqr','corvip','cornas','corurxy','corwght','corcov','sqrvip','sqrnas','sqrurxy','sqrwght','sqrcov','vipnas','vipurxy','vipwght','vipcov','nasurxy','naswght','nascov','urxywght','urxycov','wghtcov','prodall'}
                                    hopsl = chopscomb(Xint{i,1},y,hmod,vector,options);
                                    resultfeedI = feedinter(Xint{i,1},y,[hmod hopsl],vector,[pe nloop],options);
                                    Xsel (i,1) = {resultfeedI.Xs};
                                    varsI(i,1) = {resultfeedI.vs};
                                    rmsecvI(i) = resultfeedI.rmsecv;
                                case {'main','MAIN','inter','INTER','all','ALL'}
                                    hopsl = chops(Xint{i,1},y,hmod,20,vector,options);
                                    resultfeedI = feedinter(Xint{i,1},y,[hmod hopsl.hops],vector,[pe nloop],options);
                                    Xsel (i,1) = {resultfeedI.Xs};
                                    varsI(i,1) = {resultfeedI.vs};
                                    rmsecvI(i) = resultfeedI.rmsecv;
                                case {'cor','urxy','cov','COR','URXY','COV'}
                                    hopsl = hmodl;
                                    resultfeedI = feedinter(Xint{i,1},y,[hmod hopsl],vector,[pe nloop],options);
                                    Xsel (i,1) = {resultfeedI.Xs};
                                    varsI(i,1) = {resultfeedI.vs};
                                    rmsecvI(i) = resultfeedI.rmsecv;
                            end
                        case {'n','N','no','No','NO'}
                            switch vector
                                case {'reg','cor','sqr','vip','nas','urxy','wght','cov','REG','COR','SQR','VIP','NAS','URXY','WGHT','COV','regcor','regsqr','regvip','regnas','regurxy','regwght','regcov','corsqr','corvip','cornas','corurxy','corwght','corcov','sqrvip','sqrnas','sqrurxy','sqrwght','sqrcov','vipnas','vipurxy','vipwght','vipcov','nasurxy','naswght','nascov','urxywght','urxycov','wghtcov','prodall'}
                                    resultfeedI = feedinter(Xint{i,1},y,[hmod hops],vector,[pe nloop],options);
                                    Xsel (i,1) = {resultfeedI.Xs};
                                    varsI(i,1) = {resultfeedI.vs};
                                    rmsecvI(i) = resultfeedI.rmsecv;
                                case {'main','MAIN','inter','INTER','all','ALL'}
                                    resultfeedI = feedinter(Xint{i,1},y,[hmod hops],vector,[pe nloop],options);
                                    Xsel (i,1) = {resultfeedI.Xs};
                                    varsI(i,1) = {resultfeedI.vs};
                                    rmsecvI(i) = resultfeedI.rmsecv;
                                case {'custom','Custom','CUSTOM'}
                                    resultfeedI = customvec(Xint{i,1},y,hmod,'feed',[pe nloop],options,interv{i,1});
                                    Xsel (i,1) = {resultfeedI.Xs};
                                    custvec(i,1) = {resultfeedI.informative_vector};
                                    varsI(i,1) = {resultfeedI.vs};
                                    rmsecvI(i) = resultfeedI.rmsecv;
                            end
                    end
                end
            else
                Xint(p,1) = {X(:,((p-1)*var)+1:end)};
                Xsel = cell(p,1);
                custvec = cell(p,1);
                varsI = cell(p,1);
                rmsecvI = zeros(1,p);
                for i = 1:p
                    waitbarn(i/p,wb);
                    switch calhops
                        case {'y','Y','yes','Yes','YES'}
                            switch vector
                                case {'reg','sqr','vip','nas','wght','REG','SQR','VIP','NAS','WGHT'}
                                    hopsl = chops(Xint{i,1},y,hmod,20,vector,options);
                                    resultfeedI = feedinter(Xint{i,1},y,[hmod hopsl],vector,[pe nloop],options);
                                    Xsel (i,1) = {resultfeedI.Xs};
                                    varsI(i,1) = {resultfeedI.vs};
                                    rmsecvI(i) = resultfeedI.rmsecv;
                                case {'regcor','regsqr','regvip','regnas','regurxy','regwght','regcov','corsqr','corvip','cornas','corurxy','corwght','corcov','sqrvip','sqrnas','sqrurxy','sqrwght','sqrcov','vipnas','vipurxy','vipwght','vipcov','nasurxy','naswght','nascov','urxywght','urxycov','wghtcov','prodall'}
                                    hopsl = chopscomb(Xint{i,1},y,hmod,vector,options);
                                    resultfeedI = feedinter(Xint{i,1},y,[hmod hopsl.hops],vector,[pe nloop],options);
                                    Xsel (i,1) = {resultfeedI.Xs};
                                    varsI(i,1) = {resultfeedI.vs};
                                    rmsecvI(i) = resultfeedI.rmsecv;
                                case {'main','MAIN','inter','INTER','all','ALL'}
                                    hopsl = chops (Xint{i,1},y,hmod,20,vector,options);
                                    resultfeedI = feedinter(Xint{i,1},y,[hmod hopsl.hops],vector,[pe nloop],options);
                                    Xsel (i,1) = {resultfeedI.Xs};
                                    varsI(i,1) = {resultfeedI.vs};
                                    rmsecvI(i) = resultfeedI.rmsecv;
                                case {'cor','urxy','cov','COR','URXY','COV'}
                                    hopsl = hmod;
                                    resultfeedI = feedinter(Xint{i,1},y,[hmod hopsl],vector,[pe nloop],options);
                                    Xsel (i,1) = {resultfeedI.Xs};
                                    varsI(i,1) = {resultfeedI.vs};
                                    rmsecvI(i) = resultfeedI.rmsecv;
                            end
                        case {'n','N','no','No','NO'}
                            switch vector
                                case {'reg','cor','sqr','vip','nas','urxy','wght','cov','REG','COR','SQR','VIP','NAS','urxy','WGHT','COV','regcor','regsqr','regvip','regnas','regurxy','regwght','regcov','corsqr','corvip','cornas','corurxy','corwght','corcov','sqrvip','sqrnas','sqrurxy','sqrwght','sqrcov','vipnas','vipurxy','vipwght','vipcov','nasurxy','naswght','nascov','urxywght','urxycov','wghtcov','prodall'}
                                    resultfeedI = feedinter(Xint{i,1},y,[hmod hops],vector,[pe nloop],options);
                                    Xsel (i,1) = {resultfeedI.Xs};
                                    varsI(i,1) = {resultfeedI.vs};
                                    rmsecvI(i) = resultfeedI.rmsecv;
                                case {'main','MAIN','inter','INTER','all','ALL'}
                                    resultfeedI = feedinter(Xint{i,1},y,[hmod hops],vector,[pe nloop],options);
                                    Xsel (i,1) = {resultfeedI.Xs};
                                    varsI(i,1) = {resultfeedI.vs};
                                    rmsecvI(i) = resultfeedI.rmsecv;
                                case {'custom','Custom','CUSTOM'}
                                    resultfeedI = customvec(Xint{i,1},y,hmod,'feed',[pe nloop],options,interv{i,1});
                                    Xsel (i,1) = {resultfeedI.Xs};
                                    custvec(i,1) = {resultfeedI.informative_vector};
                                    varsI(i,1) = {resultfeedI.vs};
                                    rmsecvI(i) = resultfeedI.rmsecv;
                            end
                    end
                end
            end
            close(wb);
            [t,~] = size(Xsel);
            
            % Matrix
            Xsli = Xsel{1,1};
            for i=2:t
                Xsli = [Xsli,Xsel{i,1}];
            end
            
            % Custom vector
            switch vector
                case {'custom','Custom','CUSTOM'}
                    vecsli = custvec{1,1};
                    for i=2:t
                        vecsli = [vecsli;custvec{i,1}];
                    end
            end
            
            % Selected  variables
            Vsc = cell(t,1);
            for i=1:t
                Vsc(i,1) = {varsI{i,1}+(i-1)*var};
            end
            
            Vsli = Vsc{1,1};
            for i=2:t
                Vsli = [Vsli,Vsc{i,1}];
            end
            
            % Ordered Intervals
            switch calhops
                case {'y','Y','yes','Yes','YES'}
                    switch vector
                        case {'reg','sqr','vip','nas','wght','REG','SQR','VIP','NAS','WGHT'}
                            hopsl = chops(Xsli,y,hmod,20,vector,options);
                            resultfeedOI = feedinter(Xsli,y,[hmod hopsl],vector,[pe nloop],options);
                            varsOI = resultfeedOI.vs;
                            vecOI = resultfeedOI.vec;
                            resultOI = resultfeedOI;
                        case {'regcor','regsqr','regvip','regnas','regurxy','regwght','regcov','corsqr','corvip','cornas','corurxy','corwght','corcov','sqrvip','sqrnas','sqrurxy','sqrwght','sqrcov','vipnas','vipurxy','vipwght','vipcov','nasurxy','naswght','nascov','urxywght','urxycov','wghtcov','prodall'}
                            hopsl = chopscomb(Xsli,y,hmod,vector,options);
                            resultfeedOI = feedinter(Xsli,y,[hmod hopsl.hops],vector,[pe nloop],options);
                            varsOI = resultfeedOI.vs;
                            vecOI = resultfeedOI.vec;
                            resultOI = resultfeedOI;
                        case {'main','MAIN','inter','INTER','all','ALL'}
                            hopsl = chops(Xsli,y,hmod,20,vector,options);
                            resultfeedOI = feedinter(Xsli,y,[hmod hopsl.hops],vector,[pe nloop],options);
                            varsOI = resultfeedOI.vs;
                            vecOI = resultfeedOI.vec;
                            resultOI = resultfeedOI;
                        case {'cor','urxy','cov','COR','URXY','COV'}
                            hopsl = hmod;
                            resultfeedOI = feedinter(Xsli,y,[hmod hopsl],vector,[pe nloop],options);
                            varsOI = resultfeedOI.vs;
                            vecOI = resultfeedOI.vec;
                            resultOI = resultfeedOI;
                    end
                case {'n','N','no','No','NO'}
                    switch vector
                        case {'reg','cor','sqr','vip','nas','urxy','wght','cov','REG','COR','SQR','VIP','NAS','URXY','WGHT','COV','regcor','regsqr','regvip','regnas','regurxy','regwght','regcov','corsqr','corvip','cornas','corurxy','corwght','corcov','sqrvip','sqrnas','sqrurxy','sqrwght','sqrcov','vipnas','vipurxy','vipwght','vipcov','nasurxy','naswght','nascov','urxywght','urxycov','wghtcov','prodall'}
                            resultfeedOI = feedinter(Xsli,y,[hmod hops],vector,[pe nloop],options);
                            varsOI = resultfeedOI.vs;
                            vecOI = resultfeedOI.vec;
                            resultOI = resultfeedOI;
                        case {'main','MAIN','inter','INTER','all','ALL'}
                            resultfeedOI = feedinter(Xsli,y,[hmod hops],vector,[pe nloop],options);
                            varsOI = resultfeedOI.vs;
                            vecOI = resultfeedOI.vec;
                            resultOI = resultfeedOI;
                        case {'custom','Custom','CUSTOM'}
                            resultfeedOI = customvec(Xsli,y,hmod,'feed',[pe nloop],options,vecsli);
                            varsOI = resultfeedOI.vs;
                            vecOI = 'custom';
                            resultOI = resultfeedOI;
                    end
            end
            vsOI = Vsli(varsOI);
            bestOI = resultOI;
            
            infvect = bestOI.informative_vector;
            Inter.feedOPS.infofeed = bestOI.info;
            
            resulttab = [full;bestOI.resultSEL];
            rnames_result = {'Full'};
            [r,~] = size(bestOI.resultSEL);
            for i=1:r
                rnames_result{i+1,1} = strcat('Loop',num2str(i));
            end
            colnames_result = {'Subset' 'nVars' 'hmod' 'RMSECV' 'Rcv'};
            fprintf(' - - - - - - - - - - Feed iOPS - - - - - - - - - - - -');
            Quality_parameters = array2table(resulttab,'RowNames',rnames_result,'VariableNames',colnames_result)
            colnames_MP = {'Sample' 'y' 'ycv' 'Residual' 'Rel_error'};
            MP_results = array2table(bestOI.tabMP,'VariableNames',colnames_MP)
            
            % Saving results
            Inter.feedOPS.feed = resultOI;
            Inter.feedOPS.rmsecv_intervals = rmsecvI;
            Inter.feedOPS.hmod = bestOI.nlvs(1);
            Inter.feedOPS.hops = bestOI.nlvs(2:end);
            Inter.feedOPS.vs =  vsOI;
            switch vector
                case {'main','MAIN','inter','INTER','all','ALL'}
                    Inter.feedOPS.best_vec = vecOI;
                    Inter.feedOPS.informative_vector = infvect;
                case {'custom','Custom','CUSTOM'}
                    Inter.feedOPS.informative_vector = infvect';
                otherwise
                    Inter.feedOPS.informative_vector = infvect;
            end
            Inter.feedOPS.Quality_parameters = Quality_parameters;
            Inter.feedOPS.MP_results = MP_results;
            
            interfig(X,Inter);
            
        case {'auto','Auto','AUTO','autoops','autoOPS','AUTOOPS'}
            display('AutoOPS will be apllied in each interval.');
            Inter.autoOPS.type = 'auto';

            % Intervals models
            Xint = {};
            Xint(1,1) = {X(:,1:var)};
            for i = 2:ri
                Xint(i,1) = {X(:,((i-1)*var)+1:(i)*var)};
            end
            [p,~] = size(Xint);
            wb = waitbarn(0,'Interval model...');
            if dif>=50
                Xint(p+1,1) = {X(:,(nu+1:end))};
                Xsel = cell(p+1,1);
                custvec = cell(p+1,1);
                varsI = cell(p+1,1);
                rmsecvI = zeros(1,p+1);
                for i = 1:p+1
                    waitbarn(i/(p+1),wb);
                    switch calhops
                        case {'y','Y','yes','Yes','YES'}
                            switch vector
                                case {'reg','sqr','vip','nas','wght','REG','SQR','VIP','NAS','WGHT'}
                                    hopsl = chops(Xint{i,1},y,hmod,20,vector,options);
                                    resultautoI = autofeed(Xint{i,1},y,hmod,hopsl,vector,options);
                                    Xsel (i,1) = {resultautoI.Xs};
                                    varsI(i,1) = {resultautoI.varsels};
                                    rmsecvI(i) = resultautoI.resultSEL(1,4);
                                case {'regcor','regsqr','regvip','regnas','regurxy','regwght','regcov','corsqr','corvip','cornas','corurxy','corwght','corcov','sqrvip','sqrnas','sqrurxy','sqrwght','sqrcov','vipnas','vipurxy','vipwght','vipcov','nasurxy','naswght','nascov','urxywght','urxycov','wghtcov','prodall'}
                                    hopsl = chopscomb(Xint{i,1},y,hmod,vector,options);
                                    resultautoI = autofeed(Xint{i,1},y,hmod,hopsl.hops,vector,options);
                                    Xsel (i,1) = {resultautoI.Xs};
                                    varsI(i,1) = {resultautoI.varsels};
                                    rmsecvI(i) = resultautoI.resultSEL(1,4);
                                case {'main','MAIN','inter','INTER','all','ALL'}
                                    hopsl = chops(Xint{i,1},y,hmod,20,vector,options);
                                    resultautoI = autofeed(Xint{i,1},y,hmod,hopsl.hops,vector,options);
                                    Xsel (i,1) = {resultautoI.Xs};
                                    varsI(i,1) = {resultautoI.varsels};
                                    rmsecvI(i) = resultautoI.resultSEL(1,4);
                                case {'cor','urxy','cov','COR','URXY','COV'}
                                    hopsl = hmod;
                                    resultautoI = autofeed(Xint{i,1},y,hmod,hopsl,vector,options);
                                    Xsel (i,1) = {resultautoI.Xs};
                                    varsI(i,1) = {resultautoI.varsels};
                                    rmsecvI(i) = resultautoI.resultSEL(1,4);
                            end
                        case {'n','N','no','No','NO'}
                            switch vector
                                case {'reg','cor','sqr','vip','nas','urxy','wght','cov','REG','COR','SQR','VIP','NAS','URXY','WGHT','COV','regcor','regsqr','regvip','regnas','regurxy','regwght','regcov','corsqr','corvip','cornas','corurxy','corwght','corcov','sqrvip','sqrnas','sqrurxy','sqrwght','sqrcov','vipnas','vipurxy','vipwght','vipcov','nasurxy','naswght','nascov','urxywght','urxycov','wghtcov','prodall'}
                                    resultautoI = autofeed(Xint{i,1},y,hmod,hops,vector,options);
                                    Xsel (i,1) = {resultautoI.Xs};
                                    varsI(i,1) = {resultautoI.varsels};
                                    rmsecvI(i) = resultautoI.resultSEL(1,4);
                                case {'main','MAIN','inter','INTER','all','ALL'}
                                    resultautoI = autofeed(Xint{i,1},y,hmod,hops,vector,options);
                                    Xsel (i,1) = {resultautoI.Xs};
                                    varsI(i,1) = {resultautoI.varsels};
                                    rmsecvI(i) = resultautoI.resultSEL(1,4);
                                case {'custom','Custom','CUSTOM'}
                                    resultautoI = customvec(Xint{i,1},y,hmod,'auto',[pe nloop],options,interv{i,1});
                                    Xsel (i,1) = {resultautoI.Xs};
                                    custvec(i,1) = {resultautoI.informative_vector};
                                    varsI(i,1) = {resultautoI.vs};
                                    rmsecvI(i) = resultautoI.rmsecv;
                            end
                    end
                    
                end
            else
                Xint(p,1) = {X(:,((p-1)*var)+1:end)};
                Xsel = cell(p,1);
                custvec = cell(p,1);
                varsI = cell(p,1);
                rmsecvI = zeros(1,p);
                for i = 1:p
                    waitbarn(i/p,wb);
                    switch calhops
                        case {'y','Y','yes','Yes','YES'}
                            switch vector
                                case {'reg','sqr','vip','nas','wght','REG','SQR','VIP','NAS','WGHT'}
                                    hopsl = chops(Xint{i,1},y,hmod,20,vector,options);
                                    resultautoI = autofeed(Xint{i,1},y,hmod,hopsl,vector,options);
                                    Xsel (i,1) = {resultautoI.Xs};
                                    varsI(i,1) = {resultautoI.varsels};
                                    rmsecvI(i) = resultautoI.resultSEL(1,4);
                                case {'regcor','regsqr','regvip','regnas','regurxy','regwght','regcov','corsqr','corvip','cornas','corurxy','corwght','corcov','sqrvip','sqrnas','sqrurxy','sqrwght','sqrcov','vipnas','vipurxy','vipwght','vipcov','nasurxy','naswght','nascov','urxywght','urxycov','wghtcov','prodall'}
                                    hopsl = chopscomb(Xint{i,1},y,hmod,vector,options);
                                    resultautoI = autofeed(Xint{i,1},y,hmod,hopsl.hops,vector,options);
                                    Xsel (i,1) = {resultautoI.Xs};
                                    varsI(i,1) = {resultautoI.varsels};
                                    rmsecvI(i) = resultautoI.resultSEL(1,4);
                                case {'main','MAIN','inter','INTER','all','ALL'}
                                    hopsl = chops(Xint{i,1},y,hmod,20,vector,options);
                                    resultautoI = autofeed(Xint{i,1},y,hmod,hopsl.hops,vector,options);
                                    Xsel (i,1) = {resultautoI.Xs};
                                    varsI(i,1) = {resultautoI.varsels};
                                    rmsecvI(i) = resultautoI.resultSEL(1,4);
                                case {'cor','urxy','cov','COR','URXY','COV'}
                                    hopsl = hmod;
                                    resultautoI = autofeed(Xint{i,1},y,hmod,hopsl,vector,options);
                                    Xsel (i,1) = {resultautoI.Xs};
                                    varsI(i,1) = {resultautoI.varsels};
                                    rmsecvI(i) = resultautoI.resultSEL(1,4);
                            end
                        case {'n','N','no','No','NO'}
                            switch vector
                                case {'reg','cor','sqr','vip','nas','urxy','wght','cov','REG','COR','SQR','VIP','NAS','URXY','WGHT','COV','regcor','regsqr','regvip','regnas','regurxy','regwght','regcov','corsqr','corvip','cornas','corurxy','corwght','corcov','sqrvip','sqrnas','sqrurxy','sqrwght','sqrcov','vipnas','vipurxy','vipwght','vipcov','nasurxy','naswght','nascov','urxywght','urxycov','wghtcov','prodall'}
                                    resultautoI = autofeed(Xint{i,1},y,hmod,hops,vector,options);
                                    Xsel (i,1) = {resultautoI.Xs};
                                    varsI(i,1) = {resultautoI.varsels};
                                    rmsecvI(i) = resultautoI.resultSEL(1,4);
                                case {'main','MAIN','inter','INTER','all','ALL'}
                                    resultautoI = autofeed(Xint{i,1},y,hmod,hops,vector,options);
                                    Xsel (i,1) = {resultautoI.Xs};
                                    varsI(i,1) = {resultautoI.varsels};
                                    rmsecvI(i) = resultautoI.resultSEL(1,4);
                                case {'custom','Custom','CUSTOM'}
                                    resultautoI = customvec(Xint{i,1},y,hmod,'auto',[pe nloop],options,interv{i,1});
                                    Xsel (i,1) = {resultautoI.Xs};
                                    custvec(i,1) = {resultautoI.informative_vector};
                                    varsI(i,1) = {resultautoI.vs};
                                    rmsecvI(i) = resultautoI.rmsecv;
                            end
                    end
                    
                end
            end
            close(wb);
            
            [t,~] = size(Xsel);
            
            % Matrix
            Xsli = Xsel{1,1};
            for i=2:t
                Xsli = [Xsli,Xsel{i,1}];
            end
            
            % Custom vector
            switch vector
                case {'custom','Custom','CUSTOM'}
                    vecsli = custvec{1,1};
                    for i=2:t
                        vecsli = [vecsli;custvec{i,1}];
                    end
            end
            
            % Selected Variables
            Vsc = cell(t,1);
            for i=1:t
                Vsc(i,1) = {varsI{i,1}+(i-1)*var};
            end
            
            Vsli = Vsc{1,1};
            for i=2:t
                Vsli = [Vsli,Vsc{i,1}];
            end
            
            % Ordered Intervals
            switch calhops
                case {'y','Y','yes','Yes','YES'}
                    switch vector
                        case {'reg','sqr','vip','nas','wght','REG','SQR','VIP','NAS','WGHT'}
                            hopsl = chops(Xsli,y,hmod,20,vector,options);
                            resultautoOI = autofeed(Xsli,y,hmod,hopsl,vector,options);
                            varsOI = resultautoOI.varsels;
                            vecOI = resultautoOI.vec;
                            resultOI = resultautoOI;
                        case {'regcor','regsqr','regvip','regnas','regurxy','regwght','regcov','corsqr','corvip','cornas','corurxy','corwght','corcov','sqrvip','sqrnas','sqrurxy','sqrwght','sqrcov','vipnas','vipurxy','vipwght','vipcov','nasurxy','naswght','nascov','urxywght','urxycov','wghtcov','prodall'}
                            hopsl = chopscomb(Xsli,y,hmod,vector,options);
                            resultautoOI = autofeed(Xsli,y,hmod,hopsl.hops,vector,options);
                            varsOI = resultautoOI.varsels;
                            vecOI = resultautoOI.vec;
                            resultOI = resultautoOI;
                        case {'main','MAIN','inter','INTER','all','ALL'}
                            hopsl = chops(Xsli,y,hmod,20,vector,options);
                            resultautoOI = autofeed(Xsli,y,hmod,hopsl.hops,vector,options);
                            varsOI = resultautoOI.varsels;
                            vecOI = resultautoOI.vec;
                            resultOI = resultautoOI;
                        case {'cor','urxy','cov','COR','URXY','COV'}
                            hopsl = hmod;
                            resultautoOI = autofeed(Xsli,y,hmod,hopsl,vector,options);
                            varsOI = resultautoOI.varsels;
                            vecOI = resultautoOI.vec;
                            resultOI = resultautoOI;
                    end
                case {'n','N','no','No','NO'}
                    switch vector
                        case {'reg','cor','sqr','vip','nas','urxy','wght','cov','REG','COR','SQR','VIP','NAS','URXY','WGHT','COV','regcor','regsqr','regvip','regnas','regurxy','regwght','regcov','corsqr','corvip','cornas','corurxy','corwght','corcov','sqrvip','sqrnas','sqrurxy','sqrwght','sqrcov','vipnas','vipurxy','vipwght','vipcov','nasurxy','naswght','nascov','urxywght','urxycov','wghtcov','prodall'}
                            resultautoOI = autofeed(Xsli,y,hmod,hops,vector,options);
                            varsOI = resultautoOI.varsels;
                            vecOI = resultautoOI.vec;
                            resultOI = resultautoOI;
                        case {'main','MAIN','inter','INTER','all','ALL'}
                            resultautoOI = autofeed(Xsli,y,hmod,hops,vector,options);
                            varsOI = resultautoOI.varsels;
                            vecOI = resultautoOI.vec;
                            resultOI = resultautoOI;
                        case {'custom','Custom','CUSTOM'}
                            resultautoOI = customvec(Xsli,y,hmod,'auto',[pe nloop],options,vecsli);
                            varsOI = resultautoOI.vs;
                            vecOI = 'custom';
                            resultOI = resultautoOI;
                    end
            end
            
            vsOI = Vsli(varsOI);
            bestOI = resultOI;
            
            infvect = bestOI.informative_vector;
            
            resulttab = [full;bestOI.resultSEL];
            rnames_result = {'Full';'Sel'};
            colnames_result = {'Subset' 'nVars' 'hmod' 'RMSECV' 'Rcv'};
            fprintf(' - - - - - - - - - - Auto iOPS - - - - - - - - - - - -');
            Quality_parameters = array2table(resulttab,'RowNames',rnames_result,'VariableNames',colnames_result)
            colnames_MP = {'Sample' 'y' 'ycv' 'Residual' 'Rel_error'};
            MP_results = array2table(bestOI.tabMP,'VariableNames',colnames_MP)
            
            % Saving results
            Inter.autoOPS.auto = resultOI;
            Inter.autoOPS.rmsecv_intervals = rmsecvI;
            Inter.autoOPS.hmod = bestOI.nlvs(1);
            Inter.autoOPS.hops = bestOI.nlvs(2:end);
            Inter.autoOPS.vs =  vsOI;
            switch vector
                case {'main','MAIN','inter','INTER','all','ALL'}
                    Inter.autoOPS.best_vec = vecOI;
                    Inter.autoOPS.informative_vector = infvect;
                case {'custom','Custom','CUSTOM'}
                    Inter.autoOPS.informative_vector = infvect';
                otherwise
                    Inter.autoOPS.informative_vector = infvect;
            end
            Inter.autoOPS.Quality_parameters = Quality_parameters;
            Inter.autoOPS.MP_results = MP_results;
            
            interfig(X,Inter);
    end
end

end

function selOPS(X,vs,w)
%---------------------------------------------------------------------------
%  SELOPS algorithm shows the selected variables.
%---------------------------------------------------------------------------
%   Input:
%
%      X (m,n): Independent variable.
%
%     vs (m,1): Selected variables.
%
%      w (1,n): (OPTIONAL) Vector with the original variables.
%               For exemple: a vector of wavelengths.
%
%   I/O:  selOPS(X,vs);
%
%         selOPS(X,vs,w);
%
% See also: NLVOPS, NEWOPS, CREATEOPTIONS
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%
if nargin < 3
    plot(X');
    axis tight
    vline(vs,'r:')
    xlabel ('Variables')
    ylabel ('Intensity')
    
elseif nargin == 3
    w2 = w(vs);
    plot(w,X);
    axis tight
    vline(w2,'r:')
    xlabel ('Variables')
    ylabel ('Intensity')
end

end

function resultauto = autofeed(X,y,hmod,hops,vector,options)

% Secondary Function OPS: AutoOPS for running with feedOPS and iOPS!
%
% I/O: resultauto = autofeed(X,y,hmod,hops,vector,options);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%

[~,n] = size(X);
if hmod > n
    hmod = n;
end

switch vector
    case {'reg','cor','sqr','vip','nas','urxy','wght','cov','REG','COR','SQR','VIP','NAS','URXY','WGHT','COV'}
        if hops > n
            hops = n;
        end
        resultvec = vecops(X,y,hmod,hops,vector,options);
        resultauto.result.graph = resultvec.graph;
        resultauto.result.vect = resultvec.vect;
        vect = resultvec.vect';
        set = resultvec.sel(1);
        varsels = resultvec.varsel{set};
        sel = resultvec.sel;
        b = {vector};
    case {'regcor','regsqr','regvip','regnas','regurxy','regwght','regcov','corsqr','corvip','cornas','corurxy','corwght','corcov','sqrvip','sqrnas','sqrurxy','sqrwght','sqrcov','vipnas','vipurxy','vipwght','vipcov','nasurxy','naswght','nascov','urxywght','urxycov','wghtcov'}
        if hops > n
            hops = n;
        end
        resultvec = vecops(X,y,hmod,hops,vector,options);
        resultauto.result.graph = resultvec.(char(vector)).graph;
        resultauto.result.vect = resultvec.vect;
        vect = resultvec.vect';
        sel = resultvec.sel;
        varsels = resultvec.(char(vector)).setvar{1,sel(1)};
        b = {vector};
    case {'prodall'}
        f = hops > n;
        hops(f) = n;
        resultvec = vecops(X,y,hmod,hops,vector,options);
        resultauto.result.graph = resultvec.(char(vector)).graph;
        resultauto.result.vect = resultvec.vect;
        vect = resultvec.vect';
        sel = resultvec.sel;
        varsels = resultvec.(char(vector)).setvar{1,sel(1)};
        b = {vector};
    case {'main','MAIN'}
        f = hops > n;
        hops(f) = n;
        resultvec = vecops(X,y,hmod,hops,vector,options);
        resultauto.result.sel = resultvec.sel;
        resultauto.result.vectors = resultvec.vectors;
        prod = (resultvec.sel(:,4))./resultvec.sel(:,5);
        [~,i] = sort(prod,1,'ascend');
        vec = i(1);
        vect = resultvec.Xvet(i(1),:);
        set = resultvec.sel(vec);
        b = resultvec.vectors(vec);
        [~,varsels,sel] = setmain(vec,set,resultvec);
        resultauto.result.graph = resultvec.(char(b)).graph;
        resultauto.result.vect = vect;
    case {'inter','INTER'}
        f = hops > n;
        hops(f) = n;
        resultvec = vecops(X,y,hmod,hops,vector,options);
        resultauto.result.sel = resultvec.sel;
        resultauto.result.vectors = resultvec.vectors;
        prod = (resultvec.sel(:,4))./resultvec.sel(:,5);
        [~,i] = sort(prod,1,'ascend');
        vec = i(1);
        vect = resultvec.Xvet(i(1),:);
        set = resultvec.sel(vec);
        b = resultvec.vectors(vec);
        [~,varsels,sel] = setinter(vec,set,resultvec);
        resultauto.result.graph = resultvec.(char(b)).graph;
        resultauto.result.vect = vect;
    case {'all','ALL'}
        f = hops > n;
        hops(f) = n;
        resultvec = vecops(X,y,hmod,hops,vector,options);
        resultauto.result.sel = resultvec.sel;
        resultauto.result.vectors = resultvec.vectors;
        prod = (resultvec.sel(:,4))./resultvec.sel(:,5);
        [~,i] = sort(prod,1,'ascend');
        vec = i(1);
        vect = resultvec.Xvet(i(1),:);
        set = resultvec.sel(vec);
        b = resultvec.vectors(vec);
        [~,varsels,sel] = setall(vec,set,resultvec);
        resultauto.result.graph = resultvec.(char(b)).graph;
        resultauto.result.vect = vect;
end

cval = options.crossvalidation;
splits = options.split;
prep = options.preprocessing;

Xs = X(:,varsels);
[~,ns] = size(Xs);
if hmod > ns
    hmod = ns;
end
PLSmodelS = cval_pls(Xs,y,hmod,prep,cval,splits);
ycv = PLSmodelS.Ypred(:,:,hmod);
err = (((y-ycv)./y)*100);
res = y-ycv;
nu = (1:length(y))';
T = [nu y ycv res err];

resultauto.nlvs = [hmod hops];
resultauto.varsels = varsels;
resultauto.vec = char(b);
resultauto.informative_vector = vect(varsels);
resultauto.Xs = Xs;
resultauto.resultSEL = sel;
resultauto.tabMP = T;


end

function NhOPS = chops(X,y,nlv,fn,vector,options)

% Secondary Function OPS
%
% I/O: NhOPS = chops(X,y,nlv,fn,vector,options);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%

in=nlv;

[~,n] = size(X);
if fn > n
    fn = n;
end

switch vector
    case {'reg','sqr','vip','nas','wght','REG','SQR','VIP','NAS','WGHT'}
        for i = in:fn
            resultvec = vecops(X,y,nlv,i,vector,options);
            res = resultvec.sel;
            Rp(i,:) = res;
        end
        Rp = Rp(in:fn,:);
        nin = in:fn;
        [~,b] = min(Rp(:,4));
        NhOPS = nin(b);
    case{'main','MAIN','inter','INTER','all','ALL'}
        for i = in:fn
            resultreg = vecops(X,y,nlv,i,'reg',options);
            resreg = resultreg.sel;
            Rpreg(i,:) = resreg;
            resultsqr = vecops(X,y,nlv,i,'sqr',options);
            ressqr = resultsqr.sel;
            Rpsqr(i,:) = ressqr;
            resultvip = vecops(X,y,nlv,i,'vip',options);
            resvip = resultvip.sel;
            Rpvip(i,:) = resvip;
            resultnas = vecops(X,y,nlv,i,'nas',options);
            resnas = resultnas.sel;
            Rpnas(i,:) = resnas;
            resultwght = vecops(X,y,nlv,i,'wght',options);
            reswght = resultwght.sel;
            Rpwght(i,:) = reswght;
        end
        Rpreg = Rpreg(in:fn,:);
        [~,b1] = min(Rpreg(:,4));
        Rpsqr = Rpsqr(in:fn,:);
        [~,b2] = min(Rpsqr(:,4));
        Rpvip = Rpvip(in:fn,:);
        [~,b3] = min(Rpvip(:,4));
        Rpnas = Rpnas(in:fn,:);
        [~,b4] = min(Rpnas(:,4));
        Rpwght = Rpwght(in:fn,:);
        [~,b5] = min(Rpwght(:,4));
        
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [nin(b1) nin(b2) nin(b3) nin(b4) nin(b5)];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
end

end

function NhOPS = chopscomb(X,y,nlv,vector,options)

% Secondary Function OPS
%
% I/O: NhOPS = chops(X,y,nlv,fn,vector,options);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%

fn = options.fn;
in=nlv;
[~,n] = size(X);
if fn > n
    fn = n;
end
lvs = in:fn;
steps = length(lvs);

switch vector
    case {'regcor','regurxy','regcov'}
        Rp = zeros(steps,1);
        for i = 1:steps
            resultvec = vecops(X,y,nlv,lvs(i),'reg',options);
            res = resultvec.sel(4);
            Rp(i,:) = res;
        end
        nin = in:fn;
        [~,b] = min(Rp);
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        NhOPS.vec = vec;
        NhOPS.hops = [nin(b) 0 0 0 0];
    case {'corsqr','sqrurxy','sqrcov'}
        Rp = zeros(steps,1);
        for i = 1:steps
            resultvec = vecops(X,y,nlv,lvs(i),'sqr',options);
            res = resultvec.sel(4);
            Rp(i,:) = res;
        end
        nin = in:fn;
        [~,b] = min(Rp);
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        NhOPS.vec = vec;
        NhOPS.hops = [0 nin(b) 0 0 0];
    case {'corvip','vipurxy','vipcov'}
        Rp = zeros(steps,1);
        for i = 1:steps
            resultvec = vecops(X,y,nlv,lvs(i),'vip',options);
            res = resultvec.sel(4);
            Rp(i,:) = res;
        end
        nin = in:fn;
        [~,b] = min(Rp);
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        NhOPS.vec = vec;
        NhOPS.hops = [0 0 nin(b) 0 0];
    case {'cornas','nasurxy','nascov'}
        Rp = zeros(steps,1);
        for i = 1:steps
            resultvec = vecops(X,y,nlv,lvs(i),'nas',options);
            res = resultvec.sel(4);
            Rp(i,:) = res;
        end
        nin = in:fn;
        [~,b] = min(Rp);
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        NhOPS.vec = vec;
        NhOPS.hops = [0 0 0 nin(b) 0];
    case {'corwght','urxywght','wghtcov'}
        Rp = zeros(steps,1);
        for i = 1:steps
            resultvec = vecops(X,y,nlv,lvs(i),'wght',options);
            res = resultvec.sel(4);
            Rp(i,:) = res;
        end
        nin = in:fn;
        [~,b] = min(Rp);
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        NhOPS.vec = vec;
        NhOPS.hops = [0 0 0 0 nin(b)];        
    case{'prodall'}
        Rpreg = zeros(steps,1);
        Rpsqr = zeros(steps,1);
        Rpvip = zeros(steps,1);
        Rpnas = zeros(steps,1);
        Rpwght = zeros(steps,1);
        for i = 1:steps
            resultreg = vecops(X,y,nlv,lvs(i),'reg',options);
            resreg = resultreg.sel(4);
            Rpreg(i,:) = resreg;
            resultsqr = vecops(X,y,nlv,lvs(i),'sqr',options);
            ressqr = resultsqr.sel(4);
            Rpsqr(i,:) = ressqr;
            resultvip = vecops(X,y,nlv,lvs(i),'vip',options);
            resvip = resultvip.sel(4);
            Rpvip(i,:) = resvip;
            resultnas = vecops(X,y,nlv,lvs(i),'nas',options);
            resnas = resultnas.sel(4);
            Rpnas(i,:) = resnas;
            resultwght = vecops(X,y,nlv,lvs(i),'wght',options);
            reswght = resultwght.sel(4);
            Rpwght(i,:) = reswght;
        end
        [~,b1] = min(Rpreg);
        [~,b2] = min(Rpsqr);
        [~,b3] = min(Rpvip);
        [~,b4] = min(Rpnas);
        [~,b5] = min(Rpwght);
        
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [nin(b1) nin(b2) nin(b3) nin(b4) nin(b5)];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
    case {'corurxy','corcov','urxycov'}
        NhOPS.hops = [0 0 0 0 0];
    case {'regsqr'}
        Rpreg = zeros(steps,1);
        Rpsqr = zeros(steps,1);
        for i = 1:steps
            resultreg = vecops(X,y,nlv,lvs(i),'reg',options);
            resreg = resultreg.sel(4);
            Rpreg(i,:) = resreg;
            resultsqr = vecops(X,y,nlv,lvs(i),'sqr',options);
            ressqr = resultsqr.sel(4);
            Rpsqr(i,:) = ressqr;
        end
        [~,b1] = min(Rpreg);
        [~,b2] = min(Rpsqr);
        
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [nin(b1) nin(b2) 0 0 0];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
    case {'regvip'}
        Rpreg = zeros(steps,1);
        Rpvip = zeros(steps,1);
        for i = 1:steps
            resultreg = vecops(X,y,nlv,lvs(i),'reg',options);
            resreg = resultreg.sel(4);
            Rpreg(i,:) = resreg;
            resultvip = vecops(X,y,nlv,lvs(i),'vip',options);
            resvip = resultvip.sel(4);
            Rpvip(i,:) = resvip;
            
        end
        [~,b1] = min(Rpreg);
        [~,b3] = min(Rpvip);
        
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [nin(b1) 0 nin(b3) 0 0];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
    case {'regnas'}
        Rpreg = zeros(steps,1);
        Rpnas = zeros(steps,1);
        for i = 1:steps
            resultreg = vecops(X,y,nlv,lvs(i),'reg',options);
            resreg = resultreg.sel(4);
            Rpreg(i,:) = resreg;
            resultnas = vecops(X,y,nlv,lvs(i),'nas',options);
            resnas = resultnas.sel(4);
            Rpnas(i,:) = resnas;
        end
        [~,b1] = min(Rpreg);
        [~,b4] = min(Rpnas);
              
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [nin(b1) 0 0 nin(b4) 0];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
    case {'regwght'}
        Rpreg = zeros(steps,1);
        Rpwght = zeros(steps,1);
        for i = 1:steps
            resultreg = vecops(X,y,nlv,lvs(i),'reg',options);
            resreg = resultreg.sel(4);
            Rpreg(i,:) = resreg;
            resultwght = vecops(X,y,nlv,lvs(i),'wght',options);
            reswght = resultwght.sel(4);
            Rpwght(i,:) = reswght;
        end
        [~,b1] = min(Rpreg);
        [~,b5] = min(Rpwght);
        
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [nin(b1) 0 0 0 nin(b5)];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
    case {'sqrvip'}
        Rpsqr = zeros(steps,1);
        Rpvip = zeros(steps,1);
        for i = 1:steps
            resultsqr = vecops(X,y,nlv,lvs(i),'sqr',options);
            ressqr = resultsqr.sel(4);
            Rpsqr(i,:) = ressqr;
            resultvip = vecops(X,y,nlv,lvs(i),'vip',options);
            resvip = resultvip.sel(4);
            Rpvip(i,:) = resvip;
        end
        [~,b2] = min(Rpsqr);
        [~,b3] = min(Rpvip);
              
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [0 nin(b2) nin(b3) 0 0];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
    case{'sqrnas'}
        Rpsqr = zeros(steps,1);
        Rpnas = zeros(steps,1);
        for i = 1:steps
            resultsqr = vecops(X,y,nlv,lvs(i),'sqr',options);
            ressqr = resultsqr.sel(4);
            Rpsqr(i,:) = ressqr;
            resultnas = vecops(X,y,nlv,lvs(i),'nas',options);
            resnas = resultnas.sel(4);
            Rpnas(i,:) = resnas;
        end
        [~,b2] = min(Rpsqr);
        [~,b4] = min(Rpnas);
        
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [0 nin(b2) 0 nin(b4) 0];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
    case {'sqrwght'}
        Rpsqr = zeros(steps,1);
        Rpwght = zeros(steps,1);
        for i = 1:steps
            resultsqr = vecops(X,y,nlv,lvs(i),'sqr',options);
            ressqr = resultsqr.sel(4);
            Rpsqr(i,:) = ressqr;
            resultwght = vecops(X,y,nlv,lvs(i),'wght',options);
            reswght = resultwght.sel(4);
            Rpwght(i,:) = reswght;
        end
        [~,b2] = min(Rpsqr);
        [~,b5] = min(Rpwght);
        
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [0 nin(b2) 0 0 nin(b5)];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
    case {'vipnas'}
        Rpvip = zeros(steps,1);
        Rpnas = zeros(steps,1);
        for i = 1:steps
            resultvip = vecops(X,y,nlv,lvs(i),'vip',options);
            resvip = resultvip.sel(4);
            Rpvip(i,:) = resvip;
            resultnas = vecops(X,y,nlv,lvs(i),'nas',options);
            resnas = resultnas.sel(4);
            Rpnas(i,:) = resnas;
        end
        [~,b3] = min(Rpvip);
        [~,b4] = min(Rpnas);
        
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [0 0  nin(b3) nin(b4) 0];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
    case {'vipwght'}
        Rpvip = zeros(steps,1);
        Rpwght = zeros(steps,1);
        for i = 1:steps
            resultvip = vecops(X,y,nlv,lvs(i),'vip',options);
            resvip = resultvip.sel(4);
            Rpvip(i,:) = resvip;
            resultwght = vecops(X,y,nlv,lvs(i),'wght',options);
            reswght = resultwght.sel(4);
            Rpwght(i,:) = reswght;
        end
        [~,b3] = min(Rpvip);
        [~,b5] = min(Rpwght);
        
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [0 0 nin(b3) 0 nin(b5)];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
    case {'naswght'}
        Rpnas = zeros(steps,1);
        Rpwght = zeros(steps,1);
        for i = 1:steps
            resultnas = vecops(X,y,nlv,lvs(i),'nas',options);
            resnas = resultnas.sel(4);
            Rpnas(i,:) = resnas;
            resultwght = vecops(X,y,nlv,lvs(i),'wght',options);
            reswght = resultwght.sel(4);
            Rpwght(i,:) = reswght;
        end
        [~,b4] = min(Rpnas);
        [~,b5] = min(Rpwght);
        
        nin = in:fn;
        vec = {'reg' 'sqr' 'vip' 'nas' 'wght'};
        hops = [0 0 0 nin(b4) nin(b5)];
        
        NhOPS.vec = vec;
        NhOPS.hops = hops;
end

end

function resultvec = customvec(X,y,hmod,type,criteria,options,custvector)

% Secondary Function OPS: iOPS for custom vector!
%
% I/O: resultfeed = feedinter(X,y,nlvs,vector,criteria,options);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%

[~,n1] = size(X);
if hmod > n1
    hmod = n1;
end
custvector = reshape(custvector,1,length(custvector));
switch type
    case {'feed'}
        
        p = criteria(1);
        nloop = criteria(2);
        
        cval = options.crossvalidation;
        splits = options.split;
        prep = options.preprocessing;
        
        a = 1;
        PLSmodel = cval_pls(X,y,hmod,prep,cval,splits);
        e(a) = PLSmodel.RMSECV(1,hmod(1));
        
        % First selection
        [Xor,ind] = ordx(X,custvector);
        [self,resultcustomf] = perfops(Xor,y,ind,options,hmod);
        vsel = resultcustomf.setvar{1,self(1)};
        Xa = X(:,vsel);
        PLSmodelS = cval_pls(Xa,y,hmod,prep,cval,splits);
        ycv = PLSmodelS.Ypred(:,:,hmod);
        err = (((y-ycv)./y)*100);
        res = y-ycv;
        nu = (1:length(y))';
        T = {};
        T{1,1} = [nu y ycv res err];
        a = 2;
        e(a) = self(4);
        select =  self;
        SELresult (1,:) = select;
        VS = {};
        VS(1,1) = {vsel};
        custvector = custvector(vsel);
        InfVec = {};
        InfVec(1,1) = {custvector'};
        per = 1;
        counter = 2;
        while per > p
            [Xor,ind] = ordx(Xa,custvector);
            [sels,resultcustoms] = perfops(Xor,y,ind,options,hmod);
            vsel = resultcustoms.setvar{1,sels(1)};
            Xa = X(:,vsel);
            [~,n2] = size(Xa);
            if hmod > n2
                hmod = n2;
            end
            PLSmodelS = cval_pls(Xa,y,hmod,prep,cval,splits);
            ycv = PLSmodelS.Ypred(:,:,hmod);
            err = (((y-ycv)./y)*100);
            res = y-ycv;
            nu = (1:length(y))';
            T{counter,1} = [nu y ycv res err];
            a = a+1;
            rmsecv = sels(4);
            e(a) = rmsecv;
            select =  sels;
            SELresult (counter,:) = select;
            vsn = vsel;
            vso = VS{1,(counter-1)};
            vs = vso(vsn);
            VS(1,counter) = {vs};
            custvector = custvector(vsel);
            InfVec(1,counter) = {custvector'};
            ne = length(e);
            de = e(ne) - e(ne-1);
            per = abs(de)/e(ne);
            counter = counter + 1;            
            if (de >= 0) || ((counter -1) == nloop)
                break
            end
        end
        
        resultvec.vectors = InfVec;
        
        [r,~] = size(SELresult);
        if (r == nloop) || (r == 1)
            resultvec.bestloop = r;
            resultvec.info = strcat(['The best variables were selected in loop ',num2str(r),'.']);
            vsopt = VS{1,r};
            infvecopt = InfVec{1,r};
            Topt = T{r,1};
            hmod = SELresult(r,3);
            f = e(r+1);
        else
            resultvec.bestloop = r;
            resultvec.info = strcat(['The best variables were selected in loop ',num2str(r-1),'.']);
            vsopt = VS{1,r-1};
            infvecopt = InfVec{1,r-1};
            Topt = T{r-1,:};
            hmod = SELresult(r,3);
            f = e(r);
        end
        Xs = X(:,vsopt);
        
        resultvec.nlvs = hmod;
        resultvec.Xs = Xs;
        resultvec.vs = vsopt;
        resultvec.vec = 'custom';
        resultvec.informative_vector = infvecopt;
        resultvec.resultSEL = SELresult;
        resultvec.tabMP = Topt;
        resultvec.rmsecv = f;
    case {'auto'}
        [Xor,ind] = ordx(X,custvector);
        [selt,resultcustom] = perfops(Xor,y,ind,options,hmod);
        vect = custvector';
        varsels = resultcustom.setvar{1,selt(1)};
        
        resultvec.result.graph = resultcustom.graph;
        resultvec.result.vect = vect;
        cval = options.crossvalidation;
        splits = options.split;
        prep = options.preprocessing;
        
        Xs = X(:,varsels);
        [~,ns] = size(Xs);
        if hmod > ns
            hmod = ns;
        end
        
        PLSmodelS = cval_pls(Xs,y,hmod,prep,cval,splits);
        ycv = PLSmodelS.Ypred(:,:,hmod);
        err = (((y-ycv)./y)*100);
        res = y-ycv;
        nu = (1:length(y))';
        T = [nu y ycv res err];
        
        resultvec.nlvs = hmod;
        resultvec.Xs = Xs;
        resultvec.vs = varsels;
        resultvec.vec = 'custom';
        resultvec.informative_vector = vect(varsels);
        resultvec.resultSEL = selt;
        resultvec.tabMP = T;
        resultvec.rmsecv = selt(1,4);
end

end

function resultfeed = feedinter(X,y,nlvs,vector,criteria,options)

% Secondary Function OPS: FeedOPS for running with iOPS!
%
% I/O: resultfeed = feedinter(X,y,nlvs,vector,criteria,options);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%

hmod = nlvs(1);
hops = nlvs(2:end);
p = criteria(1);
nloop = criteria(2);

cval = options.crossvalidation;
splits = options.split;
prep = options.preprocessing;

a = 1;
PLSmodel = cval_pls(X,y,hmod,prep,cval,splits);
e(a) = PLSmodel.RMSECV(1,hmod(1));

% First selection
resultf = autofeed(X,y,hmod,hops,vector,options);
result = {};
result(1,1) = {resultf};
Xa = resultf.Xs;
a = 2;
e(a) = resultf.resultSEL(4);
select =  resultf.resultSEL;
SELresult = [];
SELresult (1,:) = select;
vs = resultf.varsels;
VS = {};
VS(1,1) = {vs};
Vectors = {};
Vectors(1,1) = {resultf.vec};
InfVec = {};
InfVec(1,1) = {resultf.informative_vector};
per = 1;
counter = 2;
while per > p
    [~,n] = size(Xa);
    f = hops > n;
    hops(f) = n;
    resultauto = autofeed(Xa,y,hmod,hops,vector,options);
    Xa = resultauto.Xs;
    a = a+1;
    result(1,counter) = {resultauto};
    rmsecv = resultauto.resultSEL(4);
    e(a) = rmsecv;
    select =  resultauto.resultSEL;
    SELresult (counter,:) = select;
    vsn = resultauto.varsels;
    vso = VS{1,(counter-1)};
    vs = vso(vsn);
    VS(1,counter) = {vs};
    Vectors(1,counter) = {resultauto.vec};
    InfVec(1,counter) = {resultauto.informative_vector};
    ne = length(e);
    de = e(ne) - e(ne-1);
    per = abs(de)/e(ne);
    counter = counter + 1;
    if (de >= 0) || ((counter -1) == nloop)
        break
    end
end

resultfeed.vectors = Vectors;

[r,~] = size(SELresult);
if (r == nloop) || (r == 1)
    resultfeed.bestloop = r;
    resultfeed.info = strcat(['The best variables were selected in loop ',num2str(r),'.']);
    vsopt = VS{1,r};
    infvecopt = InfVec{1,r};
    resultopt = result{1,r};
    hmod = SELresult(r,3);
    Vectors = Vectors{1,1:r};
    f = e(r+1);
else
    
    resultfeed.bestloop = r-1;
    resultfeed.info = strcat(['The best variables were selected in loop ',num2str(r-1),'.']);
    vsopt = VS{1,r-1};
    infvecopt = InfVec{1,r-1};
    resultopt = result{1,r-1};
    hmod = SELresult(r,3);
    Vectors = Vectors{1,1:r-1};
    f = e(r);
end
Xs = X(:,vsopt);
T = resultopt.tabMP;

switch vector
    case {'main','MAIN','inter','INTER','all','ALL'}
        resultfeed.nlvs = [hmod hops];
    otherwise
        resultfeed.nlvs = hmod;
end

resultfeed.vec = Vectors;
resultfeed.Xs = Xs;
resultfeed.vs = vsopt;
resultfeed.informative_vector = infvecopt;
resultfeed.resultSEL = SELresult;
resultfeed.tabMP = T;
resultfeed.rmsecv = f;

end

function vecall = allvec(X,y,hmod,hops,options)

% Vector Function OPS
%
% I/O: vecall = allvec(X,y,hmod,hops,options);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%

%Regression (REG) vector
reg = vecreg(X,y,hops(1));
[Xor,ind] = ordx(X,reg);
[selreg,resultreg] = perfops(Xor,y,ind,options,hmod);
vecall.reg = resultreg;

% Correlation (COR) vector
cor = veccorr(X,y);
[Xor,ind] = ordx(X,cor);
[selcor,resultcor] = perfops(Xor,y,ind,options,hmod);
vecall.cor = resultcor;

% Residual (SQR) vector
sqr = vecsqr(X,y,hops(2));
[Xor,ind] = ordx(X,sqr);
[selsqr,resultsqr] = perfops(Xor,y,ind,options,hmod);
vecall.sqr = resultsqr;

% Variable influence on projection (VIP) vector
vip = vecvip(X,y,hops(3));
[Xor,ind] = ordx(X,vip);
[selvip,resultvip] = perfops(Xor,y,ind,options,hmod);
vecall.vip = resultvip;

% Net analyte signal (NAS) vector
nas = vecnas(X,y,hops(4));
[Xor,ind] = ordx(X,nas);
[selnas, resultnas] = perfops(Xor,y,ind,options,hmod);
vecall.nas = resultnas;

% Error in y (URXY) vector
urxy = vecerry (X,y);
[Xor,ind] = ordx(X,urxy);
[selurxy,resulturxy] = perfops(Xor,y,ind,options,hmod);
vecall.urxy = resulturxy;

% Weight (WGHT) vector
wght = vecweight (X,y,hops(5));
[Xor,ind] = ordx(X,wght);
[selwght,resultwght] = perfops(Xor,y,ind,options,hmod);
vecall.wght = resultwght;

% Covariance procedures (COV) vector
cov = veccov (X,y);
[Xor,ind] = ordx(X,cov);
[selcov,resultcov] = perfops(Xor,y,ind,options,hmod);
vecall.cov = resultcov;

% Product REG x COR
regcor = reg.*cor;
[Xor,ind] = ordx(X,regcor);
[sel1,resultregcor] = perfops(Xor,y,ind,options,hmod);
vecall.regcor = resultregcor;

% Product REG x SQR
regsqr = reg.*sqr;
[Xor,ind] = ordx(X,regsqr);
[sel2,resultregsqr] = perfops(Xor,y,ind,options,hmod);
vecall.regsqr = resultregsqr;

% Product REG x VIP
regvip = reg.*vip;
[Xor,ind] = ordx(X,regvip);
[sel3,resultregvip] = perfops(Xor,y,ind,options,hmod);
vecall.regvip = resultregvip;

% Product REG x NAS
regnas = reg.*nas;
[Xor,ind] = ordx(X,regnas);
[sel4,resultregnas] = perfops(Xor,y,ind,options,hmod);
vecall.regnas = resultregnas;

% Product REG x URXY
regurxy = reg.*urxy;
[Xor,ind] = ordx(X,regurxy);
[sel5,resultregurxy] = perfops(Xor,y,ind,options,hmod);
vecall.regurxy = resultregurxy;

% Product REG x WGHT
regwght = reg.*wght;
[Xor,ind] = ordx(X,regwght);
[sel6,resultregwght] = perfops(Xor,y,ind,options,hmod);
vecall.regwght = resultregwght;

% Product REG x COV
regcov = reg.*cov;
[Xor,ind] = ordx(X,regcov);
[sel7,resultregcov] = perfops(Xor,y,ind,options,hmod);
vecall.regcov = resultregcov;

% Product COR x SQR
corsqr = reg.*sqr;
[Xor,ind] = ordx(X,corsqr);
[sel8,resultcorsqr] = perfops(Xor,y,ind,options,hmod);
vecall.corsqr = resultcorsqr;

% Product COR x VIP
corvip = cor.*vip;
[Xor,ind] = ordx(X,corvip);
[sel9,resultcorvip] = perfops(Xor,y,ind,options,hmod);
vecall.corvip = resultcorvip;

% Product COR x NAS
cornas = cor.*nas;
[Xor,ind] = ordx(X,cornas);
[sel10,resultcornas] = perfops(Xor,y,ind,options,hmod);
vecall.cornas = resultcornas;

% Product COR x URXY
corurxy = cor.*urxy;
[Xor,ind] = ordx(X,corurxy);
[sel11,resultcorurxy] = perfops(Xor,y,ind,options,hmod);
vecall.corurxy = resultcorurxy;

% Product COR x WGHT
corwght = cor.*wght;
[Xor,ind] = ordx(X,corwght);
[sel12,resultcorwght] = perfops(Xor,y,ind,options,hmod);
vecall.corwght = resultcorwght;

% Product COR x COV
corcov = cor.*cov;
[Xor,ind] = ordx(X,corcov);
[sel13,resultcorcov] = perfops(Xor,y,ind,options,hmod);
vecall.corcov = resultcorcov;

% Product SQR x VIP
sqrvip = sqr.*vip;
[Xor,ind] = ordx(X,sqrvip);
[sel14,resultsqrvip] = perfops(Xor,y,ind,options,hmod);
vecall.sqrvip = resultsqrvip;

% Product SQR x NAS
sqrnas = sqr.*nas;
[Xor,ind] = ordx(X,sqrnas);
[sel15,resultsqrnas] = perfops(Xor,y,ind,options,hmod);
vecall.sqrnas = resultsqrnas;

% Product SQR x URXY
sqrurxy = sqr.*urxy;
[Xor,ind] = ordx(X,sqrurxy);
[sel16,resultsqrurxy] = perfops(Xor,y,ind,options,hmod);
vecall.sqrurxy = resultsqrurxy;

% Product SQR x WGHT
sqrwght = sqr.*wght;
[Xor,ind] = ordx(X,sqrwght);
[sel17,resultsqrwght] = perfops(Xor,y,ind,options,hmod);
vecall.sqrwght = resultsqrwght;

% Product SQR x COV
sqrcov = sqr.*cov;
[Xor,ind] = ordx(X,sqrcov);
[sel18,resultsqrcov] = perfops(Xor,y,ind,options,hmod);
vecall.sqrcov = resultsqrcov;

% Product VIP x NAS
vipnas = vip.*nas;
[Xor,ind] = ordx(X,vipnas);
[sel19,resultvipnas] = perfops(Xor,y,ind,options,hmod);
vecall.vipnas = resultvipnas;

% Product VIP x URXY
vipurxy = vip.*urxy;
[Xor,ind] = ordx(X,vipurxy);
[sel20,resultvipurxy] = perfops(Xor,y,ind,options,hmod);
vecall.vipurxy = resultvipurxy;

% Product VIP x WGHT
vipwght = vip.*wght;
[Xor,ind] = ordx(X,vipwght);
[sel21,resultvipwght] = perfops(Xor,y,ind,options,hmod);
vecall.vipwght = resultvipwght;

% Product VIP x COV
vipcov = vip.*cov;
[Xor,ind] = ordx(X,vipcov);
[sel22,resultvipcov] = perfops(Xor,y,ind,options,hmod);
vecall.vipcov = resultvipcov;

% Product NAS x URXY
nasurxy = nas.*urxy;
[Xor,ind] = ordx(X,nasurxy);
[sel23,resultnasurxy] = perfops(Xor,y,ind,options,hmod);
vecall.nasurxy = resultnasurxy;

% Product NAS x WGHT
naswght = nas.*wght;
[Xor,ind] = ordx(X,naswght);
[sel24,resultnaswght] = perfops(Xor,y,ind,options,hmod);
vecall.naswght = resultnaswght;

% Product NAS x COV
nascov = nas.*cov;
[Xor,ind] = ordx(X,nascov);
[sel25,resultnascov] = perfops(Xor,y,ind,options,hmod);
vecall.nascov = resultnascov;

% Product URXY x WGHT
urxywght = urxy.*wght;
[Xor,ind] = ordx(X,urxywght);
[sel26,resulturxywght] = perfops(Xor,y,ind,options,hmod);
vecall.urxywght = resulturxywght;

% Product URXY x COV
urxycov = urxy.*cov;
[Xor,ind] = ordx(X,urxycov);
[sel27,resulturxycov] = perfops(Xor,y,ind,options,hmod);
vecall.urxycov = resulturxycov;

% Product WGHT x COV
wghtcov = wght.*cov;
[Xor,ind] = ordx(X,wghtcov);
[sel28,resultwghtcov] = perfops(Xor,y,ind,options,hmod);
vecall.wghtcov = resultwghtcov;

% Product of all vectors
prodall = reg.*cor.*sqr.*vip.*nas.*urxy.*wght.*cov;
[Xor,ind] = ordx(X,prodall);
[sel29,resultprodall] = perfops(Xor,y,ind,options,hmod);
vecall.prodall = resultprodall;

sel = [selreg;selcor;selsqr;selvip;selnas;selurxy;selwght;selcov;sel1;sel2;sel3;sel4;sel5;sel6;sel7;sel8;sel9;sel10;sel11;sel12;sel13;sel14;sel15;sel16;sel17;sel18;sel19;sel20;sel21;sel22;sel23;sel24;sel25;sel26;sel27;sel28;sel29];
vect = {'reg';'cor';'sqr';'vip';'nas';'urxy';'wght';'cov';'regcor';'regsqr';'regvip';'regnas';'regurxy';'regwght';'regcov';'corsqr';'corvip';'cornas';'corurxy';'corwght';'corcov';'sqrvip';'sqrnas';'sqrurxy';'sqrwght';'sqrcov';'vipnas';'vipurxy';'vipwght';'vipcov';'nasurxy';'naswght';'nascov';'urxywght';'urxycov';'wghtcov';'prodall'};
Xvet = [reg';cor';sqr';vip';nas';urxy';wght';cov';regcor';regsqr';regvip';regnas';regurxy';regwght';regcov';corsqr';corvip';cornas';corurxy';corwght';corcov';sqrvip';sqrnas';sqrurxy';sqrwght';sqrcov';vipnas';vipurxy';vipwght';vipcov';nasurxy';naswght';nascov';urxywght';urxycov';wghtcov';prodall'];
vecall.sel = sel;
vecall.namesel =  {'Subset' 'nVars' 'nlv' 'RMSECV' 'Rcv'};
vecall.vectors = vect;
vecall.Xvet = Xvet;

end

function vecinter = combvec(X,y,hmod,hops,options,vec)

% Vector Function OPS
%
% I/O: [vect,sel,vecinter] = combvec(X,y,hmod,hops,options);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%

switch vec
    case 1
        %Regression (REG) vector
        reg = vecreg(X,y,hops(1));
        % Correlation (COR) vector
        cor = veccorr(X,y);
        % Product REG x COR
        regcor = reg.*cor;
        [Xor,ind] = ordx(X,regcor);
        [sel1,resultregcor] = perfops(Xor,y,ind,options,hmod);
        vecinter.regcor = resultregcor;
        vecinter.sel = sel1;
        vecinter.vect = regcor';
    case 2
        %Regression (REG) vector
        reg = vecreg(X,y,hops(1));
        % Residual (SQR) vector
        sqr = vecsqr(X,y,hops(2));
        % Product REG x SQR
        regsqr = reg.*sqr;
        [Xor,ind] = ordx(X,regsqr);
        [sel2,resultregsqr] = perfops(Xor,y,ind,options,hmod);
        vecinter.regsqr = resultregsqr;
        vecinter.sel = sel2;
        vecinter.vect = regsqr';
    case 3
        %Regression (REG) vector
        reg = vecreg(X,y,hops(1));
        % Variable influence on projection (VIP) vector
        vip = vecvip(X,y,hops(3));
        % Product REG x VIP
        regvip = reg.*vip;
        [Xor,ind] = ordx(X,regvip);
        [sel3,resultregvip] = perfops(Xor,y,ind,options,hmod);
        vecinter.regvip = resultregvip;
        vecinter.sel = sel3;
        vecinter.vect = regvip';
    case 4
        %Regression (REG) vector
        reg = vecreg(X,y,hops(1));
        % Net analyte signal (NAS) vector
        nas = vecnas(X,y,hops(4));
        % Product REG x NAS
        regnas = reg.*nas;
        [Xor,ind] = ordx(X,regnas);
        [sel4,resultregnas] = perfops(Xor,y,ind,options,hmod);
        vecinter.regnas = resultregnas;
        vecinter.sel = sel4;
        vecinter.vect = regnas';
    case 5
        %Regression (REG) vector
        reg = vecreg(X,y,hops(1));
        % Error in y (URXY) vector
        urxy = vecerry (X,y);
        % Product REG x URXY
        regurxy = reg.*urxy;
        [Xor,ind] = ordx(X,regurxy);
        [sel5,resultregurxy] = perfops(Xor,y,ind,options,hmod);
        vecinter.regurxy = resultregurxy;
        vecinter.sel = sel5;
        vecinter.vect = regurxy';
    case 6
        %Regression (REG) vector
        reg = vecreg(X,y,hops(1));
        % Weight (WGHT) vector
        wght = vecweight (X,y,hops(5));
        % Product REG x WGHT
        regwght = reg.*wght;
        [Xor,ind] = ordx(X,regwght);
        [sel6,resultregwght] = perfops(Xor,y,ind,options,hmod);
        vecinter.regwght = resultregwght;
        vecinter.sel = sel6;
        vecinter.vect = regwght';
    case 7
        %Regression (REG) vector
        reg = vecreg(X,y,hops(1));
        % Covariance procedures (COV) vector
        cov = veccov (X,y);
        % Product REG x COV
        regcov = reg.*cov;
        [Xor,ind] = ordx(X,regcov);
        [sel7,resultregcov] = perfops(Xor,y,ind,options,hmod);
        vecinter.regcov = resultregcov;
        vecinter.sel = sel7;
        vecinter.vect = regcov';
    case 8
        % Correlation (COR) vector
        cor = veccorr(X,y);
        % Residual (SQR) vector
        sqr = vecsqr(X,y,hops(2));
        % Product COR x SQR
        corsqr = cor.*sqr;
        [Xor,ind] = ordx(X,corsqr);
        [sel8,resultcorsqr] = perfops(Xor,y,ind,options,hmod);
        vecinter.corsqr = resultcorsqr;
        vecinter.sel = sel8;
        vecinter.vect = corsqr';
    case 9
        % Correlation (COR) vector
        cor = veccorr(X,y);
        % Variable influence on projection (VIP) vector
        vip = vecvip(X,y,hops(3));
        % Product COR x VIP
        corvip = cor.*vip;
        [Xor,ind] = ordx(X,corvip);
        [sel9,resultcorvip] = perfops(Xor,y,ind,options,hmod);
        vecinter.corvip = resultcorvip;
        vecinter.sel = sel9;
        vecinter.vect = corvip';
    case 10
        % Correlation (COR) vector
        cor = veccorr(X,y);
        % Net analyte signal (NAS) vector
        nas = vecnas(X,y,hops(4));
        % Product COR x NAS
        cornas = cor.*nas;
        [Xor,ind] = ordx(X,cornas);
        [sel10,resultcornas] = perfops(Xor,y,ind,options,hmod);
        vecinter.cornas = resultcornas;
        vecinter.sel = sel10;
        vecinter.vect = cornas';
    case 11
        % Correlation (COR) vector
        cor = veccorr(X,y);
        % Error in y (URXY) vector
        urxy = vecerry (X,y);
        % Product COR x urxy
        corurxy = cor.*urxy;
        [Xor,ind] = ordx(X,corurxy);
        [sel11,resultcorurxy] = perfops(Xor,y,ind,options,hmod);
        vecinter.corurxy = resultcorurxy;
        vecinter.sel = sel11;
        vecinter.vect = corurxy';
    case 12
        % Correlation (COR) vector
        cor = veccorr(X,y);
        % Weight (WGHT) vector
        wght = vecweight (X,y,hops(5));
        % Product COR x WGHT
        corwght = cor.*wght;
        [Xor,ind] = ordx(X,corwght);
        [sel12,resultcorwght] = perfops(Xor,y,ind,options,hmod);
        vecinter.corwght = resultcorwght;
        vecinter.sel = sel12;
        vecinter.vect = corwght';
    case 13
        % Correlation (COR) vector
        cor = veccorr(X,y);
        % Covariance procedures (COV) vector
        cov = veccov (X,y);
        % Product COR x COV
        corcov = cor.*cov;
        [Xor,ind] = ordx(X,corcov);
        [sel13,resultcorcov] = perfops(Xor,y,ind,options,hmod);
        vecinter.corcov = resultcorcov;
        vecinter.sel = sel13;
        vecinter.vect = corcov';
    case 14
        % Residual (SQR) vector
        sqr = vecsqr(X,y,hops(2));
        % Variable influence on projection (VIP) vector
        vip = vecvip(X,y,hops(3));
        % Product SQR x VIP
        sqrvip = sqr.*vip;
        [Xor,ind] = ordx(X,sqrvip);
        [sel14,resultsqrvip] = perfops(Xor,y,ind,options,hmod);
        vecinter.sqrvip = resultsqrvip;
        vecinter.sel = sel14;
        vecinter.vect = sqrvip';
    case 15
        % Residual (SQR) vector
        sqr = vecsqr(X,y,hops(2));
        % Net analyte signal (NAS) vector
        nas = vecnas(X,y,hops(4));
        % Product SQR x NAS
        sqrnas = sqr.*nas;
        [Xor,ind] = ordx(X,sqrnas);
        [sel15,resultsqrnas] = perfops(Xor,y,ind,options,hmod);
        vecinter.sqrnas = resultsqrnas;
        vecinter.sel = sel15;
        vecinter.vect = sqrnas';
    case 16
        % Residual (SQR) vector
        sqr = vecsqr(X,y,hops(2));
        % Error in y (URXY) vector
        urxy = vecerry (X,y);
        % Product SQR x URXY
        sqrurxy = sqr.*urxy;
        [Xor,ind] = ordx(X,sqrurxy);
        [sel16,resultsqrurxy] = perfops(Xor,y,ind,options,hmod);
        vecinter.sqrurxy = resultsqrurxy;
        vecinter.sel = sel16;
        vecinter.vect = sqrurxy';
    case 17
        % Residual (SQR) vector
        sqr = vecsqr(X,y,hops(2));
        % Weight (WGHT) vector
        wght = vecweight (X,y,hops(5));
        % Product SQR x WGHT
        sqrwght = sqr.*wght;
        [Xor,ind] = ordx(X,sqrwght);
        [sel17,resultsqrwght] = perfops(Xor,y,ind,options,hmod);
        vecinter.sqrwght = resultsqrwght;
        vecinter.sel = sel17;
        vecinter.vect = sqrwght';
    case 18
        % Residual (SQR) vector
        sqr = vecsqr(X,y,hops(2));
        % Covariance procedures (COV) vector
        cov = veccov (X,y);
        % Product SQR x COV
        sqrcov = sqr.*cov;
        [Xor,ind] = ordx(X,sqrcov);
        [sel18,resultsqrcov] = perfops(Xor,y,ind,options,hmod);
        vecinter.sqrcov = resultsqrcov;
        vecinter.sel = sel18;
        vecinter.vect = sqrcov';
    case 19
        % Variable influence on projection (VIP) vector
        vip = vecvip(X,y,hops(3));
        % Net analyte signal (NAS) vector
        nas = vecnas(X,y,hops(4));
        % Product VIP x NAS
        vipnas = vip.*nas;
        [Xor,ind] = ordx(X,vipnas);
        [sel19,resultvipnas] = perfops(Xor,y,ind,options,hmod);
        vecinter.vipnas = resultvipnas;
        vecinter.sel = sel19;
        vecinter.vect = vipnas';
    case 20
        % Variable influence on projection (VIP) vector
        vip = vecvip(X,y,hops(3));
        % Error in y (URXY) vector
        urxy = vecerry (X,y);
        % Product VIP x URXY
        vipurxy = vip.*urxy;
        [Xor,ind] = ordx(X,vipurxy);
        [sel20,resultvipurxy] = perfops(Xor,y,ind,options,hmod);
        vecinter.vipurxy = resultvipurxy;
        vecinter.sel = sel20;
        vecinter.vect = vipurxy';
    case 21
        % Variable influence on projection (VIP) vector
        vip = vecvip(X,y,hops(3));
        % Weight (WGHT) vector
        wght = vecweight (X,y,hops(5));
        % Product VIP x WGHT
        vipwght = vip.*wght;
        [Xor,ind] = ordx(X,vipwght);
        [sel21,resultvipwght] = perfops(Xor,y,ind,options,hmod);
        vecinter.vipwght = resultvipwght;
        vecinter.sel = sel21;
        vecinter.vect = vipwght';
    case 22
        % Variable influence on projection (VIP) vector
        vip = vecvip(X,y,hops(3));
        % Covariance procedures (COV) vector
        cov = veccov (X,y);
        % Product VIP x COV
        vipcov = vip.*cov;
        [Xor,ind] = ordx(X,vipcov);
        [sel22,resultvipcov] = perfops(Xor,y,ind,options,hmod);
        vecinter.vipcov = resultvipcov;
        vecinter.sel = sel22;
        vecinter.vect = vipcov';
    case 23
        % Net analyte signal (NAS) vector
        nas = vecnas(X,y,hops(4));
        % Error in y (URXY) vector
        urxy = vecerry (X,y);
        % Product NAS x URXY
        nasurxy = nas.*urxy;
        [Xor,ind] = ordx(X,nasurxy);
        [sel23,resultnasurxy] = perfops(Xor,y,ind,options,hmod);
        vecinter.nasurxy = resultnasurxy;
        vecinter.sel = sel23;
        vecinter.vect = nasurxy';
    case 24
        % Net analyte signal (NAS) vector
        nas = vecnas(X,y,hops(4));
        % Weight (WGHT) vector
        wght = vecweight (X,y,hops(5));
        % Product NAS x WGHT
        naswght = nas.*wght;
        [Xor,ind] = ordx(X,naswght);
        [sel24,resultnaswght] = perfops(Xor,y,ind,options,hmod);
        vecinter.naswght = resultnaswght;
        vecinter.sel = sel24;
        vecinter.vect = naswght';
    case 25
        % Net analyte signal (NAS) vector
        nas = vecnas(X,y,hops(4));
        % Covariance procedures (COV) vector
        cov = veccov (X,y);
        % Product NAS x COV
        nascov = nas.*cov;
        [Xor,ind] = ordx(X,nascov);
        [sel25,resultnascov] = perfops(Xor,y,ind,options,hmod);
        vecinter.nascov = resultnascov;
        vecinter.sel = sel25;
        vecinter.vect = nascov';
    case 26
        % Error in y (URXY) vector
        urxy = vecerry (X,y);
        % Weight (WGHT) vector
        wght = vecweight (X,y,hops(5));
        % Product URXY x WGHT
        urxywght = urxy.*wght;
        [Xor,ind] = ordx(X,urxywght);
        [sel26,resulturxywght] = perfops(Xor,y,ind,options,hmod);
        vecinter.urxywght = resulturxywght;
        vecinter.sel = sel26;
        vecinter.vect = urxywght';
    case 27
        % Error in y (URXY) vector
        urxy = vecerry (X,y);
        % Covariance procedures (COV) vector
        cov = veccov (X,y);
        % Product URXY x COV
        urxycov = urxy.*cov;
        [Xor,ind] = ordx(X,urxycov);
        [sel27,resulturxycov] = perfops(Xor,y,ind,options,hmod);
        vecinter.urxycov = resulturxycov;
        vecinter.sel = sel27;
        vecinter.vect = urxycov';
    case 28
        % Weight (WGHT) vector
        wght = vecweight (X,y,hops(5));
        % Covariance procedures (COV) vector
        cov = veccov (X,y);
        % Product WGHT x COV
        wghtcov = wght.*cov;
        [Xor,ind] = ordx(X,wghtcov);
        [sel28,resultwghtcov] = perfops(Xor,y,ind,options,hmod);
        vecinter.wghtcov = resultwghtcov;
        vecinter.sel = sel28;
        vecinter.vect = wghtcov';
    case 29
        %Regression (REG) vector
        reg = vecreg(X,y,hops(1));
        % Correlation (COR) vector
        cor = veccorr(X,y);
        % Residual (SQR) vector
        sqr = vecsqr(X,y,hops(2));
        % Variable influence on projection (VIP) vector
        vip = vecvip(X,y,hops(3));
        % Net analyte signal (NAS) vector
        nas = vecnas(X,y,hops(4));
        % Error in y (URXY) vector
        urxy = vecerry (X,y);
        % Weight (WGHT) vector
        wght = vecweight (X,y,hops(5));
        % Covariance procedures (COV) vector
        cov = veccov (X,y);
        % Product of all vectors
        prodall = reg.*cor.*sqr.*vip.*nas.*urxy.*wght.*cov;
        [Xor,ind] = ordx(X,prodall);
        [sel29,resultprodall] = perfops(Xor,y,ind,options,hmod);
        vecinter.prodall = resultprodall;
        vecinter.sel = sel29;
        vecinter.vect = prodall';
end

end

function vecinter = intervec(X,y,hmod,hops,options)

% Vector Function OPS
%
% I/O: [Xvet,sel,vecinter] = intervec(X,y,hmod,hops,options);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%

%Regression (REG) vector
reg = vecreg(X,y,hops(1));

% Correlation (COR) vector
cor = veccorr(X,y);

% Residual (SQR) vector
sqr = vecsqr(X,y,hops(2));

% Variable influence on projection (VIP) vector
vip = vecvip(X,y,hops(3));

% Net analyte signal (NAS) vector
nas = vecnas(X,y,hops(4));

% Error in y (URXY) vector
urxy = vecerry (X,y);

% Weight (WGHT) vector
wght = vecweight (X,y,hops(5));

% Covariance procedures (COV) vector
cov = veccov (X,y);

% Product REG x COR
regcor = reg.*cor;
[Xor,ind] = ordx(X,regcor);
[sel1,resultregcor] = perfops(Xor,y,ind,options,hmod);
vecinter.regcor = resultregcor;

% Product REG x SQR
regsqr = reg.*sqr;
[Xor,ind] = ordx(X,regsqr);
[sel2,resultregsqr] = perfops(Xor,y,ind,options,hmod);
vecinter.regsqr = resultregsqr;

% Product REG x VIP
regvip = reg.*vip;
[Xor,ind] = ordx(X,regvip);
[sel3,resultregvip] = perfops(Xor,y,ind,options,hmod);
vecinter.regvip = resultregvip;

% Product REG x NAS
regnas = reg.*nas;
[Xor,ind] = ordx(X,regnas);
[sel4,resultregnas] = perfops(Xor,y,ind,options,hmod);
vecinter.regnas = resultregnas;

% Product REG x URXY
regurxy = reg.*urxy;
[Xor,ind] = ordx(X,regurxy);
[sel5,resultregurxy] = perfops(Xor,y,ind,options,hmod);
vecinter.regurxy = resultregurxy;

% Product REG x WGHT
regwght = reg.*wght;
[Xor,ind] = ordx(X,regwght);
[sel6,resultregwght] = perfops(Xor,y,ind,options,hmod);
vecinter.regwght = resultregwght;

% Product REG x COV
regcov = reg.*cov;
[Xor,ind] = ordx(X,regcov);
[sel7,resultregcov] = perfops(Xor,y,ind,options,hmod);
vecinter.regcov = resultregcov;

% Product COR x SQR
corsqr = reg.*sqr;
[Xor,ind] = ordx(X,corsqr);
[sel8,resultcorsqr] = perfops(Xor,y,ind,options,hmod);
vecinter.corsqr = resultcorsqr;

% Product COR x VIP
corvip = cor.*vip;
[Xor,ind] = ordx(X,corvip);
[sel9,resultcorvip] = perfops(Xor,y,ind,options,hmod);
vecinter.corvip = resultcorvip;

% Product COR x NAS
cornas = cor.*nas;
[Xor,ind] = ordx(X,cornas);
[sel10,resultcornas] = perfops(Xor,y,ind,options,hmod);
vecinter.cornas = resultcornas;

% Product COR x URXY
corurxy = cor.*urxy;
[Xor,ind] = ordx(X,corurxy);
[sel11,resultcorurxy] = perfops(Xor,y,ind,options,hmod);
vecinter.corurxy = resultcorurxy;

% Product COR x WGHT
corwght = cor.*wght;
[Xor,ind] = ordx(X,corwght);
[sel12,resultcorwght] = perfops(Xor,y,ind,options,hmod);
vecinter.corwght = resultcorwght;

% Product COR x COV
corcov = cor.*cov;
[Xor,ind] = ordx(X,corcov);
[sel13,resultcorcov] = perfops(Xor,y,ind,options,hmod);
vecinter.corcov = resultcorcov;

% Product SQR x VIP
sqrvip = sqr.*vip;
[Xor,ind] = ordx(X,sqrvip);
[sel14,resultsqrvip] = perfops(Xor,y,ind,options,hmod);
vecinter.sqrvip = resultsqrvip;

% Product SQR x NAS
sqrnas = sqr.*nas;
[Xor,ind] = ordx(X,sqrnas);
[sel15,resultsqrnas] = perfops(Xor,y,ind,options,hmod);
vecinter.sqrnas = resultsqrnas;

% Product SQR x URXY
sqrurxy = sqr.*urxy;
[Xor,ind] = ordx(X,sqrurxy);
[sel16,resultsqrurxy] = perfops(Xor,y,ind,options,hmod);
vecinter.sqrurxy = resultsqrurxy;

% Product SQR x WGHT
sqrwght = sqr.*wght;
[Xor,ind] = ordx(X,sqrwght);
[sel17,resultsqrwght] = perfops(Xor,y,ind,options,hmod);
vecinter.sqrwght = resultsqrwght;

% Product SQR x COV
sqrcov = sqr.*cov;
[Xor,ind] = ordx(X,sqrcov);
[sel18,resultsqrcov] = perfops(Xor,y,ind,options,hmod);
vecinter.sqrcov = resultsqrcov;

% Product VIP x NAS
vipnas = vip.*nas;
[Xor,ind] = ordx(X,vipnas);
[sel19,resultvipnas] = perfops(Xor,y,ind,options,hmod);
vecinter.vipnas = resultvipnas;

% Product VIP x URXY
vipurxy = vip.*urxy;
[Xor,ind] = ordx(X,vipurxy);
[sel20,resultvipurxy] = perfops(Xor,y,ind,options,hmod);
vecinter.vipurxy = resultvipurxy;

% Product VIP x WGHT
vipwght = vip.*wght;
[Xor,ind] = ordx(X,vipwght);
[sel21,resultvipwght] = perfops(Xor,y,ind,options,hmod);
vecinter.vipwght = resultvipwght;

% Product VIP x COV
vipcov = vip.*cov;
[Xor,ind] = ordx(X,vipcov);
[sel22,resultvipcov] = perfops(Xor,y,ind,options,hmod);
vecinter.vipcov = resultvipcov;

% Product NAS x URXY
nasurxy = nas.*urxy;
[Xor,ind] = ordx(X,nasurxy);
[sel23,resultnasurxy] = perfops(Xor,y,ind,options,hmod);
vecinter.nasurxy = resultnasurxy;

% Product NAS x WGHT
naswght = nas.*wght;
[Xor,ind] = ordx(X,naswght);
[sel24,resultnaswght] = perfops(Xor,y,ind,options,hmod);
vecinter.naswght = resultnaswght;

% Product NAS x COV
nascov = nas.*cov;
[Xor,ind] = ordx(X,nascov);
[sel25,resultnascov] = perfops(Xor,y,ind,options,hmod);
vecinter.nascov = resultnascov;

% Product URXY x WGHT
urxywght = urxy.*wght;
[Xor,ind] = ordx(X,urxywght);
[sel26,resulturxywght] = perfops(Xor,y,ind,options,hmod);
vecinter.urxywght = resulturxywght;

% Product URXY x COV
urxycov = urxy.*cov;
[Xor,ind] = ordx(X,urxycov);
[sel27,resulturxycov] = perfops(Xor,y,ind,options,hmod);
vecinter.urxycov = resulturxycov;

% Product WGHT x COV
wghtcov = wght.*cov;
[Xor,ind] = ordx(X,wghtcov);
[sel28,resultwghtcov] = perfops(Xor,y,ind,options,hmod);
vecinter.wghtcov = resultwghtcov;

% Product of all vectors
prodall = reg.*cor.*sqr.*vip.*nas.*urxy.*wght.*cov;
[Xor,ind] = ordx(X,prodall);
[sel29,resultprodall] = perfops(Xor,y,ind,options,hmod);
vecinter.prodall = resultprodall;

sel = [sel1;sel2;sel3;sel4;sel5;sel6;sel7;sel8;sel9;sel10;sel11;sel12;sel13;sel14;sel15;sel16;sel17;sel18;sel19;sel20;sel21;sel22;sel23;sel24;sel25;sel26;sel27;sel28;sel29];
vect = {'regcor';'regsqr';'regvip';'regnas';'regurxy';'regwght';'regcov';'corsqr';'corvip';'cornas';'corurxy';'corwght';'corcov';'sqrvip';'sqrnas';'sqrurxy';'sqrwght';'sqrcov';'vipnas';'vipurxy';'vipwght';'vipcov';'nasurxy';'naswght';'nascov';'urxywght';'urxycov';'wghtcov';'prodall'};
Xvet = [regcor';regsqr';regvip';regnas';regurxy';regwght';regcov';corsqr';corvip';cornas';corurxy';corwght';corcov';sqrvip';sqrnas';sqrurxy';sqrwght';sqrcov';vipnas';vipurxy';vipwght';vipcov';nasurxy';naswght';nascov';urxywght';urxycov';wghtcov';prodall'];

vecinter.sel = sel;
vecinter.namesel =  {'Subset' 'nVars' 'nlv' 'RMSECV' 'Rcv'};
vecinter.vectors = vect;
vecinter.Xvet = Xvet;

end

function vecmain = mainvec(X,y,hmod,hops,options)

% Vector Function OPS
%
% I/O: [Xvet,sel,vecmain] = mainvec(X,y,hmod,hops,options);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%

%Regression (REG) vector
reg = vecreg(X,y,hops(1));
[Xor,ind] = ordx(X,reg);
[selreg,resultreg] = perfops(Xor,y,ind,options,hmod);
vecmain.reg = resultreg;

% Correlation (COR) vector
cor = veccorr(X,y);
[Xor,ind] = ordx(X,cor);
[selcor,resultcor] = perfops(Xor,y,ind,options,hmod);
vecmain.cor = resultcor;

% Residual (SQR) vector
sqr = vecsqr(X,y,hops(2));
[Xor,ind] = ordx(X,sqr);
[selsqr,resultsqr] = perfops(Xor,y,ind,options,hmod);
vecmain.sqr = resultsqr;

% Variable influence on projection (VIP) vector
vip = vecvip(X,y,hops(3));
[Xor,ind] = ordx(X,vip);
[selvip,resultvip] = perfops(Xor,y,ind,options,hmod);
vecmain.vip = resultvip;

% Net analyte signal (NAS) vector
nas = vecnas(X,y,hops(4));
[Xor,ind] = ordx(X,nas);
[selnas, resultnas] = perfops(Xor,y,ind,options,hmod);
vecmain.nas = resultnas;

% Error in y (URXY) vector
urxy = vecerry (X,y);
[Xor,ind] = ordx(X,urxy);
[selurxy,resulturxy] = perfops(Xor,y,ind,options,hmod);
vecmain.urxy = resulturxy;

% Weight (WGHT) vector
wght = vecweight (X,y,hops(5));
[Xor,ind] = ordx(X,wght);
[selwght,resultwght] = perfops(Xor,y,ind,options,hmod);
vecmain.wght = resultwght;

% Covariance procedures (COV) vector
cov = veccov (X,y);
[Xor,ind] = ordx(X,cov);
[selcov,resultcov] = perfops(Xor,y,ind,options,hmod);
vecmain.cov = resultcov;

sel = [selreg;selcor;selsqr;selvip;selnas;selurxy;selwght;selcov];
vect = {'reg';'cor';'sqr';'vip';'nas';'urxy';'wght';'cov'};
Xvet = [reg'; cor'; sqr'; vip'; nas'; urxy'; wght'; cov'];

vecmain.sel = sel;
vecmain.namesel =  {'Subset' 'nVars' 'nlv' 'RMSECV' 'Rcv'};
vecmain.vectors = vect;
vecmain.Xvet = Xvet;

end

function [result,varsels,sel] = setall(vec,set,resultvec)

% Vector Function OPS
%
% I/O: [result,varsels,sel] = setall(vec,set,resultvec);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%

switch vec
    case 1
        varsels = resultvec.reg.setvar{set};
        sel = resultvec.sel(1,:);
        result = resultvec.reg;
    case 2
        varsels = resultvec.cor.setvar{set};
        sel = resultvec.sel(2,:);
        result = resultvec.cor;
    case 3
        varsels = resultvec.sqr.setvar{set};
        sel = resultvec.sel(3,:);   
        result = resultvec.sqr;
    case 4
        varsels = resultvec.vip.setvar{set};
        sel = resultvec.sel(4,:);
        result = resultvec.vip;
    case 5
        varsels = resultvec.nas.setvar{set};
        sel = resultvec.sel(5,:);
        result = resultvec.nas;
    case 6
        varsels = resultvec.urxy.setvar{set};
        sel = resultvec.sel(6,:);
        result = resultvec.urxy;
    case 7
        varsels = resultvec.wght.setvar{set};
        sel = resultvec.sel(7,:);
        result = resultvec.wght;
    case 8
        varsels = resultvec.cov.setvar{set};
        sel = resultvec.sel(8,:);
        result = resultvec.cov;
    case 9
        varsels = resultvec.regcor.setvar{set};
        sel = resultvec.sel(9,:);
        result = resultvec.regcor;
    case 10
        varsels = resultvec.regsqr.setvar{set};
        sel = resultvec.sel(10,:);
        result = resultvec.regsqr;
    case 11
        varsels = resultvec.regvip.setvar{set};
        sel = resultvec.sel(11,:);   
        result = resultvec.regvip;
    case 12
        varsels = resultvec.regnas.setvar{set};
        sel = resultvec.sel(12,:);
        result = resultvec.regnas;
    case 13
        varsels = resultvec.regurxy.setvar{set};
        sel = resultvec.sel(13,:);
        result = resultvec.regurxy;
    case 14
        varsels = resultvec.regwght.setvar{set};
        sel = resultvec.sel(14,:);
        result = resultvec.regwght;
    case 15
        varsels = resultvec.regcov.setvar{set};
        sel = resultvec.sel(15,:);
        result = resultvec.regcov;
    case 16
        varsels = resultvec.corsqr.setvar{set};
        sel = resultvec.sel(16,:);
        result = resultvec.corsqr;
    case 17
        varsels = resultvec.corvip.setvar{set};
        sel = resultvec.sel(17,:);
        result = resultvec.corvip;
    case 18
        varsels = resultvec.cornas.setvar{set};
        sel = resultvec.sel(18,:);
        result = resultvec.cornas;
    case 19
        varsels = resultvec.corurxy.setvar{set};
        sel = resultvec.sel(19,:);  
        result = resultvec.corurxy;
    case 20
        varsels = resultvec.corwght.setvar{set};
        sel = resultvec.sel(20,:);
        result = resultvec.corwght;
    case 21
        varsels = resultvec.corcov.setvar{set};
        sel = resultvec.sel(21,:);
        result = resultvec.corcov;
    case 22
        varsels = resultvec.sqrvip.setvar{set};
        sel = resultvec.sel(22,:);
        result = resultvec.sqrvip;
    case 23
        varsels = resultvec.sqrnas.setvar{set};
        sel = resultvec.sel(23,:);
        result = resultvec.sqrnas;
    case 24
        varsels = resultvec.sqrurxy.setvar{set};
        sel = resultvec.sel(24,:);
        result = resultvec.sqrurxy;
    case 25
        varsels = resultvec.sqrwght.setvar{set};
        sel = resultvec.sel(25,:);
        result = resultvec.sqrwght;
    case 26
        varsels = resultvec.sqrcov.setvar{set};
        sel = resultvec.sel(26,:);
        result = resultvec.sqrcov;
    case 27
        varsels = resultvec.vipnas.setvar{set};
        sel = resultvec.sel(27,:);   
        result = resultvec.vipnas;
    case 28
        varsels = resultvec.vipurxy.setvar{set};
        sel = resultvec.sel(28,:);
        result = resultvec.vipurxy;
    case 29
        varsels = resultvec.vipwght.setvar{set};
        sel = resultvec.sel(29,:);
        result = resultvec.vipwght;
    case 30
        varsels = resultvec.vipcov.setvar{set};
        sel = resultvec.sel(30,:);
        result = resultvec.vipcov;
    case 31
        varsels = resultvec.nasurxy.setvar{set};
        sel = resultvec.sel(31,:);
        result = resultvec.nasurxy;
    case 32
        varsels = resultvec.naswght.setvar{set};
        sel = resultvec.sel(32,:);
        result = resultvec.naswght;
    case 33
        varsels = resultvec.nascov.setvar{set};
        sel = resultvec.sel(33,:);
        result = resultvec.nascov;
    case 34
        varsels = resultvec.urxywght.setvar{set};
        sel = resultvec.sel(34,:);
        result = resultvec.urxywght;
    case 35
        varsels = resultvec.urxycov.setvar{set};
        sel = resultvec.sel(35,:); 
        result = resultvec.urxycov;
    case 36
        varsels = resultvec.wghtcov.setvar{set};
        sel = resultvec.sel(36,:);
        result = resultvec.wghtcov;
    case 37
        varsels = resultvec.prodall.setvar{set};
        sel = resultvec.sel(37,:);
        result = resultvec.prodall;
end

end

function [result,varsels,sel] = setinter(vec,set,resultvec)

% Vector Function OPS
%
% I/O: [result,varsels,sel] = setinter(vec,set,resultvec);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%

switch vec
    case 1
        varsels = resultvec.regcor.setvar{set};
        sel = resultvec.sel(1,:);
        result = resultvec.regcor;
    case 2
        varsels = resultvec.regsqr.setvar{set};
        sel = resultvec.sel(2,:);
        result = resultvec.regsqr;
    case 3
        varsels = resultvec.regvip.setvar{set};
        sel = resultvec.sel(3,:); 
        result = resultvec.regvip;
    case 4
        varsels = resultvec.regnas.setvar{set};
        sel = resultvec.sel(4,:);
        result = resultvec.regnas;
    case 5
        varsels = resultvec.regurxy.setvar{set};
        sel = resultvec.sel(5,:);
        result = resultvec.regurxy;
    case 6
        varsels = resultvec.regwght.setvar{set};
        sel = resultvec.sel(6,:);
        result = resultvec.regwght;
    case 7
        varsels = resultvec.regcov.setvar{set};
        sel = resultvec.sel(7,:);
        result = resultvec.regcov;
    case 8
        varsels = resultvec.corsqr.setvar{set};
        sel = resultvec.sel(8,:);
        result = resultvec.corsqr;
    case 9
        varsels = resultvec.corvip.setvar{set};
        sel = resultvec.sel(9,:);
        result = resultvec.corvip;
    case 10
        varsels = resultvec.cornas.setvar{set};
        sel = resultvec.sel(10,:);
        result = resultvec.cornas;
    case 11
        varsels = resultvec.corurxy.setvar{set};
        sel = resultvec.sel(11,:); 
        result = resultvec.corurxy;
    case 12
        varsels = resultvec.corwght.setvar{set};
        sel = resultvec.sel(12,:);
        result = resultvec.corwght;
    case 13
        varsels = resultvec.corcov.setvar{set};
        sel = resultvec.sel(13,:);
        result = resultvec.corcov;
    case 14
        varsels = resultvec.sqrvip.setvar{set};
        sel = resultvec.sel(14,:);
        result = resultvec.sqrvip;
    case 15
        varsels = resultvec.sqrnas.setvar{set};
        sel = resultvec.sel(15,:);
        result = resultvec.sqrnas;
    case 16
        varsels = resultvec.sqrurxy.setvar{set};
        sel = resultvec.sel(16,:);
        result = resultvec.sqrurxy;
    case 17
        varsels = resultvec.sqrwght.setvar{set};
        sel = resultvec.sel(17,:);
        result = resultvec.sqrwght;
    case 18
        varsels = resultvec.sqrcov.setvar{set};
        sel = resultvec.sel(18,:);
        result = resultvec.sqrcov;
    case 19
        varsels = resultvec.vipnas.setvar{set};
        sel = resultvec.sel(19,:);
        result = resultvec.vipnas;
    case 20
        varsels = resultvec.vipurxy.setvar{set};
        sel = resultvec.sel(20,:);
        result = resultvec.vipurxy;
    case 21
        varsels = resultvec.vipwght.setvar{set};
        sel = resultvec.sel(21,:);
        result = resultvec.vipwght;
    case 22
        varsels = resultvec.vipcov.setvar{set};
        sel = resultvec.sel(22,:);
        result = resultvec.vipcov;
    case 23
        varsels = resultvec.nasurxy.setvar{set};
        sel = resultvec.sel(23,:);
        result = resultvec.nasurxy;
    case 24
        varsels = resultvec.naswght.setvar{set};
        sel = resultvec.sel(24,:);
        result = resultvec.naswght;
    case 25
        varsels = resultvec.nascov.setvar{set};
        sel = resultvec.sel(25,:);
        result = resultvec.nascov;
    case 26
        varsels = resultvec.urxywght.setvar{set};
        sel = resultvec.sel(26,:);
        result = resultvec.urxywght;
    case 27
        varsels = resultvec.urxycov.setvar{set};
        sel = resultvec.sel(27,:);
        result = resultvec.urxycov;
    case 28
        varsels = resultvec.wghtcov.setvar{set};
        sel = resultvec.sel(28,:);
        result = resultvec.wghtcov;
    case 29
        varsels = resultvec.prodall.setvar{set};
        sel = resultvec.sel(29,:);
        result = resultvec.prodall;
 end


end

function [result,varsels,sel] = setmain(vec,set,resultvec)

% Vector Function OPS
%
% I/O: [result,varsels,sel] = setmain(vec,set,resultvec);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%
switch vec
    case 1
        varsels = resultvec.reg.setvar{set};
        sel = resultvec.sel(1,:);
        result = resultvec.reg;
    case 2
        varsels = resultvec.cor.setvar{set};
        sel = resultvec.sel(2,:);
        result = resultvec.cor;
    case 3
        varsels = resultvec.sqr.setvar{set};
        sel = resultvec.sel(3,:);
        result = resultvec.sqr;
    case 4
        varsels = resultvec.vip.setvar{set};
        sel = resultvec.sel(4,:);
        result = resultvec.vip;
    case 5
        varsels = resultvec.nas.setvar{set};
        sel = resultvec.sel(5,:);
        result = resultvec.nas;
    case 6
        varsels = resultvec.urxy.setvar{set};
        sel = resultvec.sel(6,:);
        result = resultvec.urxy;
    case 7
        varsels = resultvec.wght.setvar{set};
        sel = resultvec.sel(7,:);
        result = resultvec.wght;
    case 8
        varsels = resultvec.cov.setvar{set};
        sel = resultvec.sel(8,:);
        result = resultvec.cov;
end

end

function corr = veccorr(X,y)

% Vector Function OPS
%
% I/O: corr = veccorr(X,y);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%
c = corrcoef([X,y]);
[~,n] = size(X);
corr = abs(c(1:n,n+1));

end

function cov = veccov (X,y)

% Vector Function OPS
%
% I/O: cov = veccov (X,y);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%

cp = diag(X'*y*y'*X);
cov = abs(cp)/max(abs(cp));

end

function urxy = vecerry (X,y)

% Vector Function OPS
%
% I/O: urxy = vecerry (X,y);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%
[m,n] = size(X);
for i = 1:n
    x = [ones(m,1) X(:,i)];
    b = x\y;
    ye = x*b;
    ey = y-ye;
    e(:,i) = b(2)/(ey'*ey);
end
ey = abs(e');
urxy = ey./max(ey);

end

function nas = vecnas(X,y,hops)

% Vector Function OPS
%
% I/O: nas = vecnas(X,y,hops);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%
B = plsbdg(X,y,hops);
b  = abs(B(:,hops));
yn = X*b;
[m,~] = size(X);
for i = 1:m
    vnas(i,:)  = (yn(i)*pinv(b'))';
end
snas = sum(vnas)';
nas = abs(snas)/max(abs(snas));

end

function reg = vecreg(X,y,hops)

% Vector Function OPS
%
% I/O: reg = vecreg(X,y,hops);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%

B = plsbdg(X,y,hops);
b  = abs(B(:,hops));
reg = b/max(b);

end

function sqr = vecsqr(X,y,hops)

% Vector Function OPS
%
% I/O: sqr = vecsqr(X,y,hops);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%
[~,t,s,p] = plsbdg(X,y,hops);
E = X - (t(:,1:hops)*s(1:hops,1:hops)*p(:,1:hops)');
[~,n] = size(X);
sqr = zeros(n,1);
for i = 1:n
    sqr(i,1) = abs(1/(E(:,i)'*E(:,i)));
end

j = isinf(sqr);
sqr(j==1)=0;

sqr = sqr/max(sqr);


end

function vip = vecvip(X,y,hops)

% Vector Function OPS
%
% I/O: vip = vectvip(X,y,hops);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%

[T,w,b,~,P,~] = nipops(X,y,hops);
b = b(:,hops);
l = length(b);
t = sum(T.^2,1);
w_norm = (w*diag(1./sqrt(sum(w.^2,1))));
bp  = (w*inv(P'*w))\b;
SS = bp.^2.*t';
v = l*w_norm.^2*SS./sum(SS);
vip = v/max(v);

end

function wght = vecweight(X,y,hops)

% Vector Function OPS
%
% I/O: wght = vecweight(X,y,hops);
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039

%%

[~,w] = nipops(X,y,hops);
wg = zeros(size(w,1),1);
for i = 1:size(w,1)
    wg(i,1) = norm(w(i,:));
end
wght = abs(wg)/max(abs(wg));

end


function validation

DAT=datenum(date);

DATref=datenum('31-Dez-2030');

if  DATref-DAT<0
    
    error('ALGORITHM EXPIRED. PLEASE CONTACT jussararoque@gmail.com.');
    
end

end

function hhh=vline(x,in1,in2)
% function h=vline(x, linetype, label)
% 
% Draws a vertical line on the current axes at the location specified by 'x'.  Optional arguments are
% 'linetype' (default is 'r:') and 'label', which applies a text label to the graph near the line.  The
% label appears in the same color as the line.
%
% The line is held on the current axes, and after plotting the line, the function returns the axes to
% its prior hold state.
%
% The HandleVisibility property of the line object is set to "off", so not only does it not appear on
% legends, but it is not findable by using findobj.  Specifying an output argument causes the function to
% return a handle to the line, so it can be manipulated or deleted.  Also, the HandleVisibility can be 
% overridden by setting the root's ShowHiddenHandles property to on.
%
% h = vline(42,'g','The Answer')
%
% returns a handle to a green vertical line on the current axes at x=42, and creates a text object on
% the current axes, close to the line, which reads "The Answer".
%
% vline also supports vector inputs to draw multiple lines at once.  For example,
%
% vline([4 8 12],{'g','r','b'},{'l1','lab2','LABELC'})
%
% draws three lines with the appropriate labels and colors.
% 
% By Brandon Kuczenski for Kensington Labs.
% brandon_kuczenski@kensingtonlabs.com
% 8 November 2001

if length(x)>1  % vector input
    for I=1:length(x)
        switch nargin
        case 1
            linetype='r:';
            label='';
        case 2
            if ~iscell(in1)
                in1={in1};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            label='';
        case 3
            if ~iscell(in1)
                in1={in1};
            end
            if ~iscell(in2)
                in2={in2};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            if I>length(in2)
                label=in2{end};
            else
                label=in2{I};
            end
        end
        h(I)=vline(x(I),linetype,label);
    end
else
    switch nargin
    case 1
        linetype='r:';
        label='';
    case 2
        linetype=in1;
        label='';
    case 3
        linetype=in1;
        label=in2;
    end

    
    
    
    g=ishold(gca);
    hold on

    y=get(gca,'ylim');
    h=plot([x x],y,linetype);
    if length(label)
        xx=get(gca,'xlim');
        xrange=xx(2)-xx(1);
        xunit=(x-xx(1))/xrange;
        if xunit<0.8
            text(x+0.01*xrange,y(1)+0.1*(y(2)-y(1)),label,'color',get(h,'color'))
        else
            text(x-.05*xrange,y(1)+0.1*(y(2)-y(1)),label,'color',get(h,'color'))
        end
    end     

    if g==0
    hold off
    end
    set(h,'tag','vline','handlevisibility','off')
end % else

if nargout
    hhh=h;
end
end

function hhh=hline(y,in1,in2)
% function h=hline(y, linetype, label)
% 
% Draws a horizontal line on the current axes at the location specified by 'y'.  Optional arguments are
% 'linetype' (default is 'r:') and 'label', which applies a text label to the graph near the line.  The
% label appears in the same color as the line.
%
% The line is held on the current axes, and after plotting the line, the function returns the axes to
% its prior hold state.
%
% The HandleVisibility property of the line object is set to "off", so not only does it not appear on
% legends, but it is not findable by using findobj.  Specifying an output argument causes the function to
% return a handle to the line, so it can be manipulated or deleted.  Also, the HandleVisibility can be 
% overridden by setting the root's ShowHiddenHandles property to on.
%
% h = hline(42,'g','The Answer')
%
% returns a handle to a green horizontal line on the current axes at y=42, and creates a text object on
% the current axes, close to the line, which reads "The Answer".
%
% hline also supports vector inputs to draw multiple lines at once.  For example,
%
% hline([4 8 12],{'g','r','b'},{'l1','lab2','LABELC'})
%
% draws three lines with the appropriate labels and colors.
% 
% By Brandon Kuczenski for Kensington Labs.
% brandon_kuczenski@kensingtonlabs.com
% 8 November 2001

if length(y)>1  % vector input
    for I=1:length(y)
        switch nargin
        case 1
            linetype='r:';
            label='';
        case 2
            if ~iscell(in1)
                in1={in1};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            label='';
        case 3
            if ~iscell(in1)
                in1={in1};
            end
            if ~iscell(in2)
                in2={in2};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            if I>length(in2)
                label=in2{end};
            else
                label=in2{I};
            end
        end
        h(I)=hline(y(I),linetype,label);
    end
else
    switch nargin
    case 1
        linetype='r:';
        label='';
    case 2
        linetype=in1;
        label='';
    case 3
        linetype=in1;
        label=in2;
    end

    
    
    
    g=ishold(gca);
    hold on

    x=get(gca,'xlim');
    h=plot(x,[y y],linetype);
    if ~isempty(label)
        yy=get(gca,'ylim');
        yrange=yy(2)-yy(1);
        yunit=(y-yy(1))/yrange;
        if yunit<0.2
            text(x(1)+0.02*(x(2)-x(1)),y+0.02*yrange,label,'color',get(h,'color'))
        else
            text(x(1)+0.02*(x(2)-x(1)),y-0.02*yrange,label,'color',get(h,'color'))
        end
    end

    if g==0
    hold off
    end
    set(h,'tag','hline','handlevisibility','off') % this last part is so that it doesn't show up on legends
end % else

if nargout
    hhh=h;
end
end

function h = refplot(p)
%REFCURVE Add a reference curve (polynomial) to a plot.
%   REFCURVE(P) adds a line with the given polynomial coefficients 
%   to the current figure.
%   H = REFCURVE(P) returns the handle to the line object
%   in H.
%   REFCURVE with no arguments plots the function Y = 0.
%   Example: 
%       y = p(1)*x^d + p(2)*x^(d-1) + ... + p(d)*x + p(d+1)
%       Shows a polynomial of degree, d.
%       Note that p(1) goes with the highest order term.
%
%   See also POLYFIT, POLYVAL.   

%   B.A. Jones 2-3-95
%   Copyright 1993-2004 The MathWorks, Inc. 


if nargin == 0
   p = zeros(2,1);
end

xlimits = get(gca,'Xlim');
ylimits = get(gca,'Ylim');

np = get(gcf,'NextPlot');
set(gcf,'NextPlot','add');

xdat = linspace(xlimits(1),xlimits(2),100);
ydat = polyval(p,xdat);
maxy = max(ydat);
miny = min(ydat);

if maxy > ylimits(2)
  if miny < ylimits(1)
     set(gca,'YLim',[miny maxy]);
  else
     set(gca,'YLim',[ylimits(1) maxy]);
  end
else
  if miny < ylimits(1)
     set(gca,'YLim',[miny ylimits(2)]);
  end
end

if nargout == 1
   h = line(xdat,ydat);
   set(h,'LineStyle','-','Color','r');
else
   hh = line(xdat,ydat);
   set(hh,'LineStyle','-','Color','r');
end

set(gcf,'NextPlot',np);
end