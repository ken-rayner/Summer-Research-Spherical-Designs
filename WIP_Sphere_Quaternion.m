function M = spherequaternionfactory(n, m)
% Returns a manifold struct to optimize over unit-norm quaternion matrices.
%
% function M = spherequaternionfactory(n)
% function M = spherequaternionfactory(n, m)
%
% Manifold of n-by-m quaternion matrices of unit Frobenius norm.
% By default, m = 1, which corresponds to the unit sphere in Q^n. The
% metric is such that the sphere is a Riemannian submanifold of the space
% of 4nx4m real matrices with the usual trace inner product, i.e., the
% usual metric.
% 
% See also: spherefactory

% This file is adapted from Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 
%
%   Sep. 4, 2014 (NB):
%       Added ehess2rhess.
%
%   April 7, 2015 (NB):
%       Added vec/mat pair (for use with hessianspectrum, for example).
%
%   April 13, 2015 (NB):
%       Added logarithm
%
%   Oct. 8, 2016 (NB)
%       Code for exponential was simplified to only treat the zero vector
%       as a particular case.
%
%   Oct. 22, 2016 (NB)
%       Distance function dist now significantly more accurate for points
%       within 1e-7 and less from each other.
%
%%   January 5, 2024 (NB)
%       Started adapting source for the case of the quaternions

    
    if ~exist('m', 'var')
        m = 1;
    end

    if m == 1
        M.name = @() sprintf('Quaternion sphere S^%d', n-1);
    else
        M.name = @() sprintf('Unit F-norm %dx%d Quaternion matrices', n, m);
    end
    
    M.dim = @() 4*(n*m)-1;
    
    M.inner = @(x, d1, d2) parts(d1(:)'*d2(:));
    
    M.norm = @(x, d) norm(d, 'fro');%EDIT
    
    M.dist = @(x, y) real(2*asin(.5*norm(x - y, 'fro')));%EDIT
    
    M.typicaldist = @() pi;%EDIT
    
    M.proj = @(x, d) reshape(d(:) - x(:)*(real(x(:)'*d(:))), n, m);%EDIT
    
    % For Riemannian submanifolds, converting a Euclidean gradient into a
    % Riemannian gradient amounts to an orthogonal projection.
    M.egrad2rgrad = M.proj;%EDIT
    
    M.ehess2rhess = @ehess2rhess;%EDIT
    function rhess = ehess2rhess(x, egrad, ehess, u)
        rhess = M.proj(x, ehess) - real((x(:)'*egrad(:)))*u;%EDIT
    end
    
    M.tangent = M.proj;%EDIT
    
    M.exp = @exponential;%EDIT
    
    M.retr = @retraction;%EDIT

    M.log = @logarithm;%EDIT
    function v = logarithm(x1, x2)
        v = M.proj(x1, x2 - x1);%EDIT
        di = M.dist(x1, x2);%EDIT
        % If the two points are "far apart", correct the norm.
        if di > 1e-6
            nv = norm(v, 'fro');%EDIT
            v = v * (di / nv);%EDIT
        end
    end
    
    M.hash = @(x) ['z' hashmd5([real(x(:)) ; imag(x(:))])];%EDIT
    
    M.rand = @() random(n, m);%EDIT
    
    M.randvec = @(x) randomvec(n, m, x);%EDIT
    
    M.lincomb = @matrixlincomb;%EDIT
    
    M.zerovec = @(x) zeros(n, m);%EDIT
    
    M.transp = @(x1, x2, d) M.proj(x2, d);%EDIT
    
    M.pairmean = @pairmean;%EDIT
    function y = pairmean(x1, x2)
        y = x1+x2;%EDIT
        y = y / norm(y, 'fro');%EDIT
    end

    mn = m*n;%EDIT
    M.vec = @(x, u_mat) [real(u_mat(:)) ; imag(u_mat(:))];%EDIT
    M.mat = @(x, u_vec) reshape(u_vec(1:mn), m, n) + 1i*reshape(u_vec((mn+1):end), m, n);%EDIT
    M.vecmatareisometries = @() true;%EDIT

end

% Exponential on the sphere
function y = exponential(x, d, t)%EDIT

    if nargin == 2%EDIT
        % t = 1;
        td = d;%EDIT
    else
        td = t*d;%EDIT
    end
    
    nrm_td = norm(td, 'fro');%EDIT
    
    if nrm_td > 0%EDIT
        y = x*cos(nrm_td) + td*(sin(nrm_td)/nrm_td);%EDIT
    else
        y = x;%EDIT
    end

end

% Retraction on the sphere
function y = retraction(x, d, t)%EDIT

    if nargin == 2%EDIT
        t = 1;%EDIT
    end
    
    y = x+t*d;%EDIT
    y = y/norm(y, 'fro');%EDIT

end

% Uniform random sampling on the sphere.
function x = random(n, m)%EDIT

    x = randn(n, m) + 1i*randn(n, m);%EDIT
    x = x/norm(x, 'fro');%EDIT

end

% Random normalized tangent vector at x.
function d = randomvec(n, m, x)%EDIT

    d = randn(n, m) + 1i*randn(n, m);%EDIT
    d = reshape(d(:) - x(:)*(real(x(:)'*d(:))), n, m);%EDIT
    d = d / norm(d, 'fro');%EDIT

end