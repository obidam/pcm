function fh=confratefig(y, t,h)
%CONFFIG Display a confusion matrix.
%
%	Description
%	CONFFIG(Y, T) displays the confusion matrix  and classification
%	performance for the predictions mat{y} compared with the targets T.
%	The data is assumed to be in a 1-of-N encoding, unless there is just
%	one column, when it is assumed to be a 2 class problem with a 0-1
%	encoding.  Each row of Y and T corresponds to a single example.
%
%	In the confusion matrix, the rows represent the true classes and the
%	columns the predicted classes.
%
%	FH = CONFFIG(Y, T) also returns the figure handle FH which  can be
%	used, for instance, to delete the figure when it is no longer needed.
%
%	See also
%	CONFMAT, DEMTRAIN
%

%	Copyright (c) Ian T Nabney (1996-2001)

[C, rate] = confmat(y, t);

% classif rates
Cr = C ./ (sum(C')' * ones(1,size(C,2)) );

% calcul de l'erreur moyenne
x = [1:size(C,1)]'*ones(1,size(C,1));x=(x-x').^2;
err = sqrt(mean(sum((x.*Cr)')));

% calssif rate [+1,-1]
x = [1:size(C,1)]'*ones(1,size(C,1));ind=find(abs(x-x')<=1);
x = zeros(size(C,1)); x(ind) = 1;
C2 = sum((Cr.*x)');
x  = diag(ones(size(C,1),1));
C1 = mean(sum((Cr.*x)'));

% $$$ fh = figure('Name', 'Confusion matrix', ...
% $$$   'NumberTitle', 'off');

% figure
figure(h);
set(h,'Name', 'Confusion matrix');
set(h,'NumberTitle', 'off')
plotmat(0.1*floor(1000*[Cr C2']), 'k', 'k', 14);
tit = ['R = ' num2str(rate(1)) '% - RC =',num2str(100*C1),'% -  mse = ',num2str(err),' - R[-1,+1] =',num2str(mean(100*C2(:))),'%'];
title(tit, 'FontSize', 14);

fh = h;