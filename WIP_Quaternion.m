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
    
    M.norm = @(x, d) norm(d(:), 'fro');%EDIT
    
    M.dist = @(x, y) norm(x(:)-y(:), 'fro');%EDIT
    
    M.typicaldist = @() sqrt(prod(dimensions_vec));%EDIT
    
    M.proj = @(x, d) d;%EDIT
    
    M.egrad2rgrad = @(x, g) g;%EDIT
    
    M.ehess2rhess = @(x, eg, eh, d) eh;%EDIT
    
    M.tangent = M.proj;%EDIT
    
    M.exp = @exp;%EDIT
    function y = exp(x, d, t)%EDIT
        if nargin == 3%EDIT
            y = x + t*d;%EDIT
        else
            y = x + d;%EDIT
        end
    end
    
    M.retr = M.exp;%EDIT
    
    M.log = @(x, y) y-x;%EDIT

    M.hash = @(x) ['z' hashmd5([real(x(:)) ; imag(x(:))])];%EDIT
    
    M.rand = @() (randn(dimensions_vec) + 1i*randn(dimensions_vec))/sqrt(2);%EDIT
    
    M.randvec = @randvec;%EDIT
    function u = randvec(x) %#ok<INUSD>
        u = randn(dimensions_vec) + 1i*randn(dimensions_vec);%EDIT
        u = u / norm(u(:), 'fro');%EDIT
    end
    
    M.lincomb = @matrixlincomb;%EDIT
    
    M.zerovec = @(x) zeros(dimensions_vec);%EDIT
    
    M.transp = @(x1, x2, d) d;%EDIT
    
    M.pairmean = @(x1, x2) .5*(x1+x2);%EDIT
    
    sz = prod(dimensions_vec);%EDIT
    M.vec = @(x, u_mat) [real(u_mat(:)) ; imag(u_mat(:))];%EDIT
    M.mat = @(x, u_vec) reshape(u_vec(1:sz), dimensions_vec) ...
                        + 1i*reshape(u_vec((sz+1):end), dimensions_vec);%EDIT
    M.vecmatareisometries = @() true;%EDIT
    M.lie_identity = @() zeros(dimensions_vec);%EDIT

end