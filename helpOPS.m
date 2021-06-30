%---------------------------------------------------------------------------
% NLVOPS: calculates the hmod (the number of latent
% variables to build PLS models) and the hops (number of latent
% variables used by OPS algorithm) at the same time.
%
% NEWOPS: performs the variable selection using the new OPS
% approaches: autoOPS, feedOPS and iOPS.
%
% SELOPS: shows the selected variables.
%
% CREATEOPTIONS: creates the options struct.
%
% OPSPLOTS: generates figure plots of OPS using OPS struct output of NEWOPS.
%---------------------------------------------------------------------------
%
% See also: NLVOPS, NEWOPS, SELOPS, CREATEOPTIONS, OPSPLOTS
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 30/06/2021
%
% J.V. Roque, W. Cardoso, L.A. Peternelli, and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039