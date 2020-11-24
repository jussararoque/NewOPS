%---------------------------------------------------------------------------
%  CREATEOPTIONS gerenates a options struct to be used in OPS algorithms.
%---------------------------------------------------------------------------
% Input:
%      win: Number of initial variables (window).
%
%      incr: Increment of variables to be added over the window.
%
%      perc: Percentage of the variables to be considered in the selection.
%
%      cv: Cross-validation (CV) options:
%          'full': leave-one-out;
%          'cblocks': contiguous block;
%          'vblinds': venetian blinds;
%          'random': random subset.
%
%      split: Number of subsets for CV;
%
%      prep_X: preprocessing for X.
%      prep_y: preprocessing for y.
%          Preprocessing options (for X and y blocks):
%           'mean': mean center;
%           'auto': autoscale;
%           'none': no preprocessing.
%
%      criteria: Vector with 2 inputs [p nloop] used by feedOPS.
%             p: percentage difference between two consecutive RMSECVS.
%             nloop: number of loops.
%
%      var: Number of variables included in each iOPS interval.
%
%      type: The option to apply autoOPS or feedOPS in each iOPS interval.
%           'auto': autoOPS;
%           'feed': feedOPS;

%      calhops: 'no' or 'yes'. The option to perform (or not) the
%            calculation of hOPS for each interval.
%
%   Output:
%
%      options: Struct with the settings to perform OPS.
%
%      I/O:  options = createoptions(win,incr,perc,cv,split,prep_X,prep_y);
%
%            options = createoptions(win,incr,perc,cv,split,prep_X,prep_y,criteria);
%
%            options = createoptions(win,incr,perc,cv,split,prep_X,prep_y,criteria,var,type,calhops);
%
% See also: NLVOPS, NEWOPS, SELOPS, OPSPLOTS
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli, and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039