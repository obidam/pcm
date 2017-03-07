function [post, a] = gmmLINpost(mix, x, W)
%GMMPOST Computes the class posterior probabilities of a Gaussian mixture model.
%
%	Description
%	This function computes the posteriors POST (i.e. the probability of
%	each component conditioned on the data P(J|X)) for a Gaussian mixture
%	model.   The data structure MIX defines the mixture model, while the
%	matrix X contains the data vectors.  Each row of X represents a
%	single vector.
%
%	See also
%	GMM, GMMACTIV, GMMPROB
%
%	Copyright (c) Ian T Nabney (1996-2001)

% Check that inputs are consistent
errstring = consist(mix, 'gmm', x);
if ~isempty(errstring)
  error(errstring);
end

ndata = size(x, 1);
a     = zeros(ndata, mix.ncentres); 
normal = (2*pi)^(mix.nin/2);
for j = 1:mix.ncentres
  diffs = x - (mix.centres(j, :) * W(j).w')';
  %% a(:, j) = exp(-0.5*sum(diffs.*diffs,2) / mix.covars(j)) ./ (normal*sqrt(mix.covars(j)));%% modif 27/11/09
  C     = inv(diag( mix.covars(j,:)));%% modif 27/11/09
  a(:, j) = exp(-0.5*(diffs' * C *diffs))./(normal*sqrt(det(C))); %% modif 27/11/09
end


post = (ones(ndata, 1)*mix.priors).*a;
s = sum(post, 2);
if any(s==0)
   warning('Some zero posterior probabilities')
   % Set any zeros to one before dividing
   zero_rows = find(s==0);
   s = s + (s==0);
   post(zero_rows, :) = 1/mix.ncentres;
end
post = post./(s*ones(1, mix.ncentres));
