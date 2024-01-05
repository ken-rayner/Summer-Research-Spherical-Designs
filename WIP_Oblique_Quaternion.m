function M = obliquequaternionfactory(n, m, transposed)
% Returns a manifold struct defining quaternion matrices w/ unit-norm columns.
%
% function M = obliquequaternionfactory(n, m)
% function M = obliquequaternionfactory(n, m, transposed)
%
% Oblique manifold: deals with quaternion matrices of size n x m such that
% each column has unit 4-norm, i.e., is a point on the unit sphere in Q^n.
% The geometry is a product geometry of m unit spheres in Q^n. For the
% metric, Q^n is treated as R^(4n), so that the real part and imaginary
% parts are treated separately as 4n real coordinates. As such, the complex
% oblique manifold is a Riemannian submanifold of (R^4)^(n x m), with the
% usual metric <u, v> = real(u'*v).
%
% If transposed is set to true (it is false by default), then the matrices
% are transposed: a point Y on the manifold is a matrix of size m x n and
% each row has unit 4-norm. It is the same geometry, just a different
% representation.
%
% In transposed form, a point Y is such that Y*Y' is a Hermitian, positive
% semidefinite matrix of size m and of rank at most n, such that all the
% diagonal entries are equal to 1.
%
%% REMOVED NOT RELEVANT TO QUATERNIONS%%
%
% See also: spherecomplexfactory complexcirclefactory obliquefactory

% This file is adapted from Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Sep. 3, 2014.
% Contributors:
% Change log:
%
%   Oct. 21, 2016 (NB)
%       Formatted for inclusion in Manopt release.
%
%   July 20, 2017 (NB)
%       Distance function is now accurate for close-by points. See notes
%       inside the spherefactory file for details. Also improves distances
%       computation as part of the log function.
%
%   May 28, 2023 (NB)
%       Fixed bug in M.log in case 'transposed' is true (bug reported by
%       Lingping Kong).
%
%%   January 5, 2024 (NB)
%       Started adapting source for the case of the quaternions


    if ~exist('transposed', 'var') || isempty(transposed)
        transposed = false;
    end

    if transposed
        trnsp = @(X) X.';
    else
        trnsp = @(X) X;
    end

    M.name = @() sprintf('Quaternion oblique manifold COB(%d, %d)', n, m);

    M.dim = @() (2*n-1)*m; %EDIT

    M.inner = @(x, d1, d2) real(d1(:)'*d2(:));%EDIT

    M.norm = @(x, d) norm(d(:));%EDIT

    M.dist = @(x, y) norm(real(2*asin(.5*sqrt(sum(trnsp(abs(x - y).^2), 1)))));%EDIT

    M.typicaldist = @() pi*sqrt(m);%EDIT

    M.proj = @(X, U) trnsp(projection(trnsp(X), trnsp(U)));%EDIT

    M.tangent = M.proj;%EDIT

    % For Riemannian submanifolds, converting a Euclidean gradient into a
    % Riemannian gradient amounts to an orthogonal projection.
    M.egrad2rgrad = M.proj;%EDIT

    M.ehess2rhess = @ehess2rhess;%EDIT
    function rhess = ehess2rhess(X, egrad, ehess, U)%EDIT
        X = trnsp(X);%EDIT
        egrad = trnsp(egrad);%EDIT
        ehess = trnsp(ehess);%EDIT
        U = trnsp(U);%EDIT

        PXehess = projection(X, ehess);%EDIT
        inners = sum(real(conj(X).*egrad), 1);%EDIT
        rhess = PXehess - bsxfun(@times, U, inners);%EDIT

        rhess = trnsp(rhess);%EDIT
    end

    M.exp = @exponential;%EDIT
    % Exponential on the complex oblique manifold
    function y = exponential(x, d, t)%EDIT
        x = trnsp(x);%EDIT
        d = trnsp(d);%EDIT

        if nargin == 2%EDIT
            % t = 1;%EDIT
            td = d;%EDIT
        else
            td = t*d;%EDIT
        end

        nrm_td = sqrt(sum(real(td).^2 + imag(td).^2, 1));%EDIT

        y = bsxfun(@times, x, cos(nrm_td)) + ...%EDIT
            bsxfun(@times, td, sinxoverx(nrm_td));%EDIT

        y = trnsp(y);%EDIT
    end

    M.log = @logarithm;%EDIT
    function v = logarithm(x1, x2)%EDIT
        x1 = trnsp(x1);%EDIT
        x2 = trnsp(x2);%EDIT

        v = projection(x1, x2 - x1);%EDIT
        dists = real(2*asin(.5*sqrt(sum(abs(x2 - x1).^2, 1))));%EDIT
        norms = sqrt(sum(real(v).^2 + imag(v).^2, 1));%EDIT
        factors = dists./norms;%EDIT
        % For very close points, dists is almost equal to norms, but
        % because they are both almost zero, the division above can return
        % NaN's. To avoid that, we force those ratios to 1.
        factors(dists <= 1e-10) = 1;%EDIT
        v = bsxfun(@times, v, factors);%EDIT

        v = trnsp(v);%EDIT
    end

    M.retr = @retraction;%EDIT
    % Retraction on the oblique manifold
    function y = retraction(x, d, t)%EDIT
        x = trnsp(x);%EDIT
        d = trnsp(d);%EDIT

        if nargin < 3%EDIT
            td = d;%EDIT
        else
            td = t*d;%EDIT
        end

        y = normalize_columns(x + td);%EDIT

        y = trnsp(y);%EDIT
    end

    M.hash = @(x) ['z' hashmd5([real(x(:)) ; imag(x(:))])];%EDIT

    M.rand = @() trnsp(random(n, m));%EDIT

    M.randvec = @(x) trnsp(randomvec(n, m, trnsp(x)));%EDIT

    M.lincomb = @matrixlincomb;%EDIT

    M.zerovec = @(x) trnsp(zeros(n, m));%EDIT

    M.transp = @(x1, x2, d) M.proj(x2, d);%EDIT

    M.pairmean = @pairmean;%EDIT
    function y = pairmean(x1, x2)%EDIT
        y = trnsp(x1+x2);%EDIT
        y = normalize_columns(y);%EDIT
        y = trnsp(y);%EDIT
    end

    % vec returns a vector representation of an input tangent vector which
    % is represented as a matrix. mat returns the original matrix
    % representation of the input vector representation of a tangent
    % vector. vec and mat are thus inverse of each other. They are
    % furthermore isometries between a subspace of R^2nm and the tangent
    % space at x.
    vect = @(X) X(:);%EDIT
    M.vec = @(x, u_mat) [vect(real(trnsp(u_mat))) ; ...%EDIT
                         vect(imag(trnsp(u_mat)))];%EDIT
    M.mat = @(x, u_vec)    trnsp(reshape(u_vec(1:(n*m)),     [n, m])) + ...%EDIT
                        1i*trnsp(reshape(u_vec((n*m+1):end), [n, m]));%EDIT
    M.vecmatareisometries = @() true;%EDIT

end

% Given a matrix X, returns the same matrix but with each column scaled so
% that they have unit 2-norm.
function X = normalize_columns(X)%EDIT
    norms = sqrt(sum(real(X).^2 + imag(X).^2, 1));%EDIT
    X = bsxfun(@times, X, 1./norms);%EDIT
end

% Orthogonal projection of the ambient vector H onto the tangent space at X
function PXH = projection(X, H)%EDIT

    % Compute the inner product between each vector H(:, i) with its root
    % point X(:, i), that is, real(X(:, i)' * H(:, i)).
    % Returns a row vector.
    inners = real(sum(conj(X).*H, 1));%EDIT

    % Subtract from H the components of the H(:, i)'s that are parallel to
    % the root points X(:, i).
    PXH = H - bsxfun(@times, X, inners);%EDIT

end

% Uniform random sampling on the sphere.
function x = random(n, m)%EDIT

    x = normalize_columns(randn(n, m) + 1i*randn(n, m));%EDIT

end

% Random normalized tangent vector at x.
function d = randomvec(n, m, x)%EDIT

    d = randn(n, m) + 1i*randn(n, m);%EDIT
    d = projection(x, d);%EDIT
    d = d / norm(d(:));%EDIT

end