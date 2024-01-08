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
%
%%   January 7, 2024 (NB)
%       Completed initial adaptation to quaternions
    
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
    
    M.norm = @(x, d) norm(norm(d(:)));
    
    M.dist = @(x, y) norm(norm(x(:)-y(:)));%Need to verify this

    M.typicaldist = @() pi;
    
    M.proj = @(x, d) reshape(d(:) - x(:)*(real(x(:)'*d(:))), n, m);%VERIFY
    
    % For Riemannian submanifolds, converting a Euclidean gradient into a
    % Riemannian gradient amounts to an orthogonal projection.
    M.egrad2rgrad = M.proj;
    
    M.ehess2rhess = @ehess2rhess;%VERIFY
    function rhess = ehess2rhess(x, egrad, ehess, u)
        rhess = M.proj(x, ehess) - real((x(:)'*egrad(:)))*u;%VERIFY
    end
    
    M.tangent = M.proj;
    
    M.exp = @exponential;
    
    M.retr = @retraction;

    M.log = @logarithm;
    function v = logarithm(x1, x2)
        v = M.proj(x1, x2 - x1);
        di = M.dist(x1, x2);
        % If the two points are "far apart", correct the norm.
        if di > 1e-6
            nv = norm(norm(v));
            v = v * (di / nv);
        end
    end
    
    M.hash = @qhash;
    function mhash = qhash(x);
        [x_r, x_i, x_j, x_k] = parts(x);
        mhash = ['q' hashmd5([x_r(:); x_i(:); x_j(:); x_k(:)])];
    end
    
    M.rand = @() random(n, m);
    
    M.randvec = @(x) randomvec(n, m, x);
    
    M.lincomb = @matrixlincomb;
    
    M.zerovec = @(x) zeros(n, m);
    
    M.transp = @(x1, x2, d) M.proj(x2, d);
    
    M.pairmean = @pairmean;
    function y = pairmean(x1, x2)
        y = x1+x2;
        y = y / norm(norm(y));
    end

    mn = m*n;

    M.vec = @mvec;%Should verify that this is correct
    function u = mvec(x, u_mat)
        [u_r, u_i, u_j, u_k] = parts(u_mat);
        u = [u_r(:); u_i(:); u_j(:); u_k(:)];
    end

    M.mat = @mmat;%Should verify that this is correct
    function X = mmat(x, u_vec)
        X = quaternion(reshape(u_vec(1:mn),[m(:)', n(:)']),reshape(u_vec(mn+1:2*mn),[m(:)', n(:)']),reshape(u_vec(2*mn+1:3*mn),[m(:)', n(:)']),reshape(u_vec(3*mn+1:end),[m(:)', n(:)']));
    end

    M.vecmatareisometries = @() true;

end

% Exponential on the sphere
function y = exponential(x, d, t)

    if nargin == 2
        % t = 1;
        td = d;
    else
        td = t*d;
    end
    
    nrm_td = norm(norm(td));
    
    if nrm_td > 0
        y = x*cos(nrm_td) + td*(sin(nrm_td)/nrm_td); %Need to verify if this is correct
    else
        y = x;%Need to verify if this is correct
    end

end

% Retraction on the sphere
function y = retraction(x, d, t)%EDIT

    if nargin == 2%EDIT
        t = 1;%EDIT
    end
    
    y = x+t*d;%EDIT
    y = y/norm(norm(y));

end

% Uniform random sampling on the sphere.
function x = random(n, m)%Verify that this is appropriate

    x = quaternion(randn([m(:)', n(:)']),randn([m(:)', n(:)']),randn([m(:)', n(:)']),randn([m(:)', n(:)']));
    x = x/norm(norm(x));

end

% Random normalized tangent vector at x.
function d = randomvec(n, m, x)

    d = quaternion(randn([m(:)', n(:)']),randn([m(:)', n(:)']),randn([m(:)', n(:)']),randn([m(:)', n(:)']));
    d = reshape(d(:) - x(:)*(real(x(:)'*d(:))), n, m);%EDIT
    d = d / norm(norm(d));

end