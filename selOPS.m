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
% See also: NLVOPS, NEWOPS, CREATEOPTIONS, OPSPLOTS
%
% Copyright: Jussara V. Roque, 2019.
% Checked by JVR: 28/05/2019
%
% J.V. Roque, W. Cardoso, L.A. Peternelli and R.F. Teófilo.
% Comprehensive new approaches for variable selection using ordered
% predictors selection. Analytica Chimica Acta (2019).
% https://doi.org/10.1016/j.aca.2019.05.039