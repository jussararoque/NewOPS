%---------------------------------------------------------------------------
%  NLVOPS algorithm calculates the hmod (the number of latent
%  variables to build PLS models) and the hops (number of latent
%  variables used by OPS algorithm) at the same time.
%---------------------------------------------------------------------------
% Input:
%
%    X (m,n): Independent variable.
%
%    y (m,1): Dependent variable.
%
%       hmod: (OPTIONAL) Number of latent variables of model. If the hmod
%       is already calculated, please enter with it.
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
%          'random': random subset).
%      Split: {10} number of subsets for CV;
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
%
%   Output:
%
%      hmod: number of latent variables of model.
%      hops: number of latent variables used by OPS algorithm.
%       vec: hops identification.
%
%   I/O:  [hmod,hops] = nlvOPS(X,y)
%
%         [hmod,hops] = nlvOPS(X,y,hmod)
%
%         [hmod,hops] = nlvOPS(X,y,[],vector)
%
%         [hmod,hops] = nlvOPS(X,y,hmod,vector)
%
%         [hmod,hops] = nlvOPS(X,y,[],vector,options)
%
%         [hmod,hops] = nlvOPS(X,y,hmod,vector,options)
%
% See also: NEWOPS, SELOPS, CREATEOPTIONS, OPSPLOTS
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli, and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039