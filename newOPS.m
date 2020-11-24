%---------------------------------------------------------------------------
%  The NEWOPS algorithm performs the variable selection using the new OPS
%  approaches: autoOPS, feedOPS and iOPS.
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
%              Vector with 1 or 5 elements. [reg sqr vip nas wght] are the
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
%           'iOPS': run iOPS;
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
%      Window: {10}  number of initial variables.
%      Increment: {5} increment of variables to be added over the window.
%      Percentage: {100} percentage of the variables to be considered in
%                  the selection.
%      cv: [ 'full' | 'cblocks' | {'vblinds'} | 'random' ]
%          A cell containing standard cross-validation (CV).
%          cross-validation options:
%          'full': leave-one-out;
%          'cblocks': contiguous block;
%          'vblinds': venetian blinds;
%          'random': random subset.
%      Split: {10} number of subsets for CV;
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
%      Var: Number of variables included in each iOPS interval.
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
%      OPS : Struct with the results of each approach: best informative
%      vector, selected variables, quality parameters (RMSECV and Rcv), and
%      predicted values of Ycv.
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
% See also: NLVOPS, SELOPS, CREATEOPTIONS, OPSPLOTS
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli, and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039