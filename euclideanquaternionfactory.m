function M = euclideanquaternionfactory(m, n)
% Returns a manifold struct to optimize over quaternion matrices.
%
% function M = euclideanquaternionfactory(m)
% function M = euclideanquaternionfactory(m, n)
% function M = euclideanquaternionfactory([n1, n2, ...])
%
% Returns M, a structure describing the vector space of quaternion matrices,
% as a manifold for Manopt.
%
% The quaternion plane is here viewed as R^4. The inner product between two
% m-by-n matrices A and B is given by: real(trace(A'*B)). This choice
% guides the proper definition of gradient and Hessian for this geometry.
% This is not the classical Euclidean inner product for complex matrices;
% it is a real inner product.
%
% See also: euclideanfactory

% This file is adapted from Manopt: www.manopt.org.
% Original author: Nicolas Boumal, April 7, 2015.
% Contributors: 
% Change log: 
%
%   Jan. 25, 2017 (NB):
%       Added functionality to handle multidimensional arrays.
%
%%   January 5, 2024 (NB)
%       Started adapting source for the case of the quaternions
%%   January 6, 2024 (NB)
%       Finished adapting source for the case of the quaternions

    % The size can be defined using both m and n, or simply with m.
    % If m is a scalar, then n is implicitly 1.
    % This mimics the use of built-in Matlab functions such as zeros(...).
    if ~exist('n', 'var') || isempty(n)
        if numel(m) == 1
            n = 1;
        else
            n = [];
        end
    end
    
    dimensions_vec = [m(:)', n(:)']; % We have a row vector.
    
    M.size = @() dimensions_vec;

    M.name = @() sprintf('Euclidean space Q^(%s)', num2str(dimensions_vec));
    
    M.dim = @() 4*prod(dimensions_vec);
    
    M.inner = @(x, d1, d2) parts(d1(:)'*d2(:));
    
    M.norm = @(x, d) norm(norm(d(:)));
    
    M.dist = @(x, y) norm(norm(x(:)-y(:)));%Need to verify this
    
    M.typicaldist = @() sqrt(prod(dimensions_vec));%No change required?
    
    M.proj = @(x, d) d;%No change required?
    
    M.egrad2rgrad = @(x, g) g;%No change required?
    
    M.ehess2rhess = @(x, eg, eh, d) eh;%No change required?
    
    M.tangent = M.proj;%No change required?
    
    M.exp = @exp;%No change required?
    function y = exp(x, d, t)
        if nargin == 3
            y = x + t*d;
        else
            y = x + d;
        end
    end
    
    M.retr = M.exp;%No change required?
    
    M.log = @(x, y) y-x;%No change required?

    M.hash = @qhash;
    function mhash = qhash(x);
        [x_r, x_i, x_j, x_k] = parts(x);
        mhash = ['q' hashmd5([x_r(:); x_i(:); x_j(:); x_k(:)])];
    end
    
    M.rand = @() quaternion(randn(dimensions_vec),randn(dimensions_vec),randn(dimensions_vec),randn(dimensions_vec))/2;
    
    M.randvec = @randvec;%Not sure if this is correct
    function u = randvec(x) %#ok<INUSD>
        u = quaternion(randn(dimensions_vec),randn(dimensions_vec),randn(dimensions_vec),randn(dimensions_vec));
        u = u / norm(norm(u(:)));
    end
    
    M.lincomb = @matrixlincomb;%No change required?
    
    M.zerovec = @(x) zeros(dimensions_vec);%No change required?
    
    M.transp = @(x1, x2, d) d;%No change required?
    
    M.pairmean = @(x1, x2) .5*(x1+x2);%No change required?
    
    sz = prod(dimensions_vec);%No change required?

    M.vec = @mvec;%Should verify that this is correct
    function u = mvec(x, u_mat)
        [u_r, u_i, u_j, u_k] = parts(u_mat);
        u = [u_r(:); u_i(:); u_j(:); u_k(:)];
    end

    M.mat = @mmat;%Should verify that this is correct
    function X = mmat(x, u_vec)
        X = quaternion(reshape(u_vec(1:sz),dimensions_vec),reshape(u_vec(sz+1:2*sz),dimensions_vec),reshape(u_vec(2*sz+1:3*sz),dimensions_vec),reshape(u_vec(3*sz+1:end),dimensions_vec));
    end

    M.vecmatareisometries = @() true;%No change required?
    M.lie_identity = @() zeros(dimensions_vec);%No change required?

end
