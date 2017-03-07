function [mix, options, errlog] = gmmem(mix, x, options,dim_1,dim_2)
%GMMEM	EM algorithm for Gaussian mixture model.
%
%	Description
%	[MIX, OPTIONS, ERRLOG] = GMMEM(MIX, X, OPTIONS) uses the Expectation
%	Maximization algorithm of Dempster et al. to estimate the parameters
%	of a Gaussian mixture model defined by a data structure MIX. The
%	matrix X represents the data whose expectation is maximized, with
%	each row corresponding to a vector.    The optional parameters have
%	the following interpretations.
%
%	OPTIONS(1) is set to 1 to display error values; also logs error
%	values in the return argument ERRLOG. If OPTIONS(1) is set to 0, then
%	only warning messages are displayed.  If OPTIONS(1) is -1, then
%	nothing is displayed.
%
%	OPTIONS(3) is a measure of the absolute precision required of the
%	error function at the solution. If the change in log likelihood
%	between two steps of the EM algorithm is less than this value, then
%	the function terminates.
%
%	OPTIONS(5) is set to 1 if a covariance matrix is reset to its
%	original value when any of its singular values are too small (less
%	than MIN_COVAR which has the value eps).   With the default value of
%	0 no action is taken.
%
%	OPTIONS(14) is the maximum number of iterations; default 100.
%
%	The optional return value OPTIONS contains the final error value
%	(i.e. data log likelihood) in OPTIONS(8).
%
%	See also
%	GMM, GMMINIT
%

%	Copyright (c) Ian T Nabney (1996-2001)

% Check that inputs are consistent
errstring = consist(mix, 'gmm', x);
if ~isempty(errstring)
  error(errstring);
end

[ndata, xdim] = size(x);

% Sort out the options
if (options(14))
  niters = options(14);
else
  niters = 100;
end

display = options(1);
store = 0;
if (nargout > 2)
  store = 1;	% Store the error values to return them
  errlog = zeros(1, niters).*NaN;
end
test = 0;
if options(3) > 0.0
  test = 1;	% Test log likelihood for termination
end

check_covars = 0;
if options(5) >= 1
  if display >= 0
    disp('check_covars is on');
  end
  check_covars = 1;	% Ensure that covariances don't collapse
  MIN_COVAR = eps;	% Minimum singular value of covariance matrix
  init_covars = mix.covars;
end

if options(19) == 1
	ffland;iw=2;jw=1;ipl=0;
	% Compute raw data pdf to fit:
	tx = linspace(min(range(x)),max(range(x)),100);
	Opdf = hist(x,tx);
	Opdf = Opdf/sum(Opdf*diff(tx(1:2))); % Scale obs pdf to make it unit integral
	tx_init = 0;
end% if 

% Main loop of algorithm
for n = 1:niters

	% Calculate posteriors based on old parameters
	[post, act] = gmmpost(mix, x);

	% Debug plot
	if options(19) == 1 & rem(n,250) == 0

		ipl=0; %clf;
		if exist('tx_ph','var'),delete(tx_ph),end% if 
		if exist('tx_ph2','var'),delete(tx_ph2),end% if 
		clear tx_pdf tx_act tx_pdfc leg		
		tx_act = gmmactiv(mix,tx')';   % GMM activations
		tx_pdf = tx_act' * (mix.priors)';  % GMM pdf
		for ic = 1 : mix.ncentres
			tx_pdfc(ic,:) = mix.priors(ic)*tx_act(ic,:);	
			leg(ic) = {sprintf('\\lambda_%i=%0.2f',ic,mix.priors(ic))};
%			n(ic) = length(find(olabels==ic)); 
		end% for ic
		if ~exist('pdfc500','var') & n == 500
			pdfc500 = tx_pdfc;
		else
			pdfc500 = zeros(size(tx_pdfc)).*NaN;
		end% if 

		ipl=ipl+1;subp(ipl)=subplot(iw,jw,ipl);hold on
		tx_ph  = plot(tx,tx_pdfc,'color','b','linewidth',2);
		tx_ph2 = plot(tx,tx_pdf,'color','k','linewidth',2);	
		plot(tx,pdfc500,'b--');
		
		if ~tx_init
			plot(tx,Opdf,'color','r','linewidth',2);	
			set(gca,'xlim',range(tx));
			grid on, box on
			la(1)=xlabel('data');la(2)=ylabel('probability density function');
		end% if 
		%legend(tx_ph,leg,'location','northwestoutside');	
		set([la gca],'fontsize',12)
		title(sprintf('Error %0.4f\niter = %i / %i\nGMM details for K=%i',errlog(n-1),n,niters,mix.ncentres),'fontweight','bold','fontsize',14)

		ipl=ipl+1;subp(ipl)=subplot(iw,jw,ipl);
		plot(errlog,'linewidth',2);
		xlim([1 niters]);
		grid on,box on
		title(sprintf('Error value as negative log likelihood of data'),'fontweight','bold','fontsize',14)

		drawnow
		tx_init = 1;
	end% if 


	% Calculate error value if needed
	if (display | store | test)
		prob = act*(mix.priors)';
		% Error value is negative log likelihood of data
		e = - sum(log(prob));
		if store
			errlog(n) = e;
		end% if 
		if display > 0
			fprintf(1, 'Cycle %4d  Error %11.6f\n', n, e);
		end% if 
		if test
			if (n > 1 & abs(e - eold) < options(3))
				options(8) = e;
				return;
			else
				eold = e;
			end% if 
		end% if 
	end% if 
  
	% Adjust the new estimates for the parameters
	post   = (post + 1e-10) ./ (sum(post,2)*ones(1,mix.ncentres));
	new_pr = sum(post, 1);
	new_c = post' * x;

	% Now move new estimates to old parameter vectors
	mix.priors = new_pr ./ ndata;

	mix.centres = new_c ./ (new_pr' * ones(1, mix.nin));

	switch mix.covar_type
	
		case 'spherical'
			n2 = dist2(x, mix.centres);
			for j = 1:mix.ncentres
				v(j) = (post(:,j)'*n2(:,j));
			end
			mix.covars = ((v./new_pr))./mix.nin;
			if check_covars
				% Ensure that no covariance is too small
				for j = 1:mix.ncentres
					if mix.covars(j) < MIN_COVAR
						mix.covars(j) = init_covars(j);
					end
				end
			end

		case 'diag'
			for j = 1:mix.ncentres
				diffs = x - (ones(ndata, 1) * mix.centres(j,:));
				mix.covars(j,:) = sum((diffs.*diffs).*(post(:,j)*ones(1,mix.nin)), 1)./new_pr(j);
			end
			if check_covars
				% Ensure that no covariance is too small
				for j = 1:mix.ncentres
					if min(mix.covars(j,:)) < MIN_COVAR
						mix.covars(j,:) = init_covars(j,:);
					end
				end
			end

		case 'mix' % This case is similar to 'diag' but has to be used when more than two non-standardized variables are used
			for j = 1:mix.ncentres
				diffs = x - (ones(ndata, 1) * mix.centres(j,:));
				mix.covars(j,:) = sum((diffs.*diffs).*(post(:,j)*ones(1,mix.nin)), 1)./new_pr(j);
			end
			if check_covars
				% Ensure that no covariance is too small
				for j = 1:mix.ncentres
					if min(mix.covars(j,:)) < MIN_COVAR
						mix.covars(j,:) = init_covars(j,:);
					end
				end
			end
			for j = 1:mix.ncentres
				mix.covars(j,1:dim_1) = ones(1,dim_1)*(sum(mix.covars(j,1:dim_1))/dim_1);
				mix.covars(j,dim_1+1:dim_1+dim_2)=ones(1,dim_2)*(sum(mix.covars(j,dim_1+1:dim_1+dim_2))/dim_2);
			end% for j
	
		case 'full'
			for j = 1:mix.ncentres
				diffs = x - (ones(ndata, 1) * mix.centres(j,:));
				diffs = diffs.*(sqrt(post(:,j))*ones(1, mix.nin));
				mix.covars(:,:,j) = (diffs'*diffs)/new_pr(j);
			end
			if check_covars
				% Ensure that no covariance is too small
				for j = 1:mix.ncentres
					if min(svd(mix.covars(:,:,j))) < MIN_COVAR
						mix.covars(:,:,j) = init_covars(:,:,j);
					end
				end
			end

		case 'ppca'
			for j = 1:mix.ncentres
				diffs = x - (ones(ndata, 1) * mix.centres(j,:));
				diffs = diffs.*(sqrt(post(:,j))*ones(1, mix.nin));
				[tempcovars, tempU, templambda] = ...
					ppca((diffs'*diffs)/new_pr(j), mix.ppca_dim);
				if length(templambda) ~= mix.ppca_dim
					error('Unable to extract enough components');
				else 
					mix.covars(j) = tempcovars;
					mix.U(:, :, j) = tempU;
					mix.lambda(j, :) = templambda;
				end
			end
			if check_covars
				if mix.covars(j) < MIN_COVAR
					mix.covars(j) = init_covars(j);
				end
			end

		otherwise
			error(['Unknown covariance type ', mix.covar_type]);               
	end% switch 


end% for n

options(8) = -sum(log(gmmprob(mix, x)));
if (display >= 0)
	disp(maxitmess);
end
  