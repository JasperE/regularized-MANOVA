function [eigval, eigvec,stdval] = invwb_prod(x,hyp,lambda_cor,lambda_var,w)
% Computation of eigenvectors and eigenvalues of the generalized eigenvalue problem Ba=bWa,
% where W is a shrinkage estimator of the within-group covariance matrix. The shrinkage constant
% lambda is estimated using the Ledoit-Wolf theorem.
%
% To do: add correct weighting for residuals

[n,~] = size(x);
alpha = -1; % Estimate W^alpha;

if nargin < 5 || isempty(w)
    w = ones(n,1)*(1/n);% Can be adjusted later if I want to put in weighting for samples.
end

if isempty(lambda_var) % Estimate correlation shrinkage
    lambda_var = estimate_lambda_var(x,w); 
elseif lambda_var > 1
    lambda_var = 1;
elseif lambda_var < 0
    lambda_var = 0;
end

if isempty(lambda_cor) % Estimate correlation shrinkage
    lambda = estimate_lambda(x,w);
else 
    lambda = lambda_cor;
    if lambda > 1
        lambda = 1;
    elseif lambda < 0 
        lambda = 0;
    end
end

% Shrinkage scale factors
[vs, ~] = svar(x,lambda_var, w);
sc = sqrt(vs*(n-1));
stdval = ones(n,1)*sc; % Matrix for scaling of mean vectors 

% Scale hypothesis matrix (hyp) by target variances - include special when
% lambda is 0
if lambda == 0
    hyps = hyp./stdval;
else
    hyps = hyp./(stdval./sqrt(lambda^alpha)); %Scale by means chosen target variances as well as correlation shrinkage parameter
end

% Autoscale residual data (x)
mx = mean(x);
sx = x - ones(n, 1)*mx;
sd = ones(n, 1)*std(sx);
xs = sx./sd; 

% Dimension reduction
sum_data = hyps+xs;
[~,svdxss,svdxsv] = svdecon(sum_data,1e-7);
m = length(diag(svdxss)); % Rank of data (should be n-1)

% Between matrix
Fb = (hyps*svdxsv)'*(hyps*svdxsv);

% Inverse within matrix
h1 = 1/(1-sum(w.*w)); % Bias correction factor; For w = 1/n this equals the usual h1 = n/(n-1)
xsscr = xs*svdxsv;
C = xsscr'*diag(w)*xsscr;
C = (1-lambda)*h1*C;
C = (C+C')/2;
if lambda == 0 % Use eigenvalue decomposition computing the matrix power
    Fw = mpower(C,alpha);
	WB = Fw*Fb;
else
    Fw = eye(m)-mpower(C/lambda+eye(m),alpha);
	WB = Fb - Fw*Fb;
end

[a, b] = eig(WB);
[eigvec, eigval] = cdf2rdf(a,b);
eigvec = svdxsv*eigvec;
eigval = diag(eigval);
[~, b] = sort(eigval,'descend');
eigval = eigval(b); eigvec = eigvec(:,b);
[zero_pos,~] = find(eigval == 0);
eigval(zero_pos) = []; eigvec(:,zero_pos) = [];
