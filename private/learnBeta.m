function beta = learnBeta( X, Y )
% betas = learnBeta( X, Y )
% X should be N-by-D, Y should be N-by-NumProblems.
% Solves subproblems independently, optionally in parallel.

D = size(X,2); P = size(Y,2); beta = zeros(D,P);
par = ~isempty(gcp('nocreate'));
if(par), parfor i=1:P, beta(:,i) = solve(X,Y(:,i)); end
else for i=1:P, beta(:,i) = solve(X,Y(:,i)); end; end

function beta = solve( X, y )
% Solve subproblem k using lsqcurvefit with analytic derivatives
[N,D] = size(X); asserteq( size(y), [N 1] );
opts = optimoptions('lsqcurvefit','Jacobian','on','Display','off');
beta = lsqcurvefit( @myfun, ones(D,1), X, double(y), [], [], opts );

function [F,J] = myfun(beta,X)
% size( beta )   = [D 1] 
% size(  X   )   = [N D]
% size( X*beta ) = [N 1]
[N,D] = size(X);
F = 1-exp(-X*beta);
F = double(F);
if nargout > 1
  J = repmat(1-F,[1 D]) .* X;
  J = double(J);
end
