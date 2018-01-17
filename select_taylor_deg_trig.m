function  [M,mv,alpha,unA] = ...
           select_taylor_deg_trig(A,b,flag,m_max,p_max,prec,shift,bal,force_estm)
%SELECT_TAYLOR_DEGREE_TRIG   Select degree of Taylor approximation.
%   [M,MV,alpha,unA] = SELECT_TAYLOR_DEGREE_TRIG(A,m_max,p_max) forms a matrix M
%   for use in determining the truncated Taylor series degree in FUNMV
%   based on parameters m_max and p_max.
%   MV is the number of matrix-vector products with A or A^* computed.

%   Reference: A. H. Al-Mohy, A New Algorithm for Computing the Action 
%   of Trigonometric and Hyperbolic Matrix Functions

%   Awad H. Al-Mohy, August 08, 2017.
if nargin < 9, force_estm = false; end
if nargin < 8 || isempty(bal), bal = false; end
if nargin < 7 || isempty(shift), shift = false; end
if nargin < 5 || isempty(p_max), p_max = 5; end
if nargin < 4 || isempty(m_max), m_max = 25; end

if p_max < 2 || m_max > 25 || m_max  < p_max*(p_max - 1)
    error('>>> Invalid p_max or m_max.')
end
n = length(A);
sigma = 1;
if ~flag, sigma = 1/2; end
if bal
    [D B] = balance(A);
    if norm(B,1) < norm(A,1), A = B; end
end
if nargin < 6 || isempty(prec), prec = class(A); end
switch prec
    case 'double'
        load theta_taylor_trig_double
    case 'single'
        load theta_taylor_trig_single
    case 'half'
        load theta_taylor_trig_half
end
if shift
    mu = trace(A)/n;
    mu = full(mu); % The code could be slower without the full!
   
    A = A - mu*speye(n);    
end
mv = 0; 
bound_hold = 0;
normA = norm(A,1); 
ell = 2;
bound = 2*ell*p_max*(p_max+3)/(m_max*size(b,2)) - 1;
c = normA^sigma;
if ~force_estm && c <= theta(m_max)*bound 
    % Base choice of m on normA, not the alpha_p.
    unA = 0; bound_hold = 1;   
    alpha = c*ones(p_max-1,1);
end

if ~force_estm && flag && ~bound_hold       
    [c,cost_d2] = normAm(A,2); 
    c = c^(1/2);
    mv = mv + cost_d2;
    if c <=  theta(m_max)*(bound - cost_d2) 
        unA = 1; bound_hold = 1;
        alpha = c*ones(p_max-1,1);
    end
end

if ~bound_hold
    eta = zeros(p_max,1); alpha = zeros(p_max-1,1);
    for p = 1:p_max   
        [c,k] = normAm(A,2*sigma*(p+1));
        c = c^(1/(2*p+2));
        mv = mv + k;
        eta(p) = c;
    end
    for p = 1:p_max-1
        alpha(p) = max(eta(p),eta(p+1));
    end
    unA = 2;
end

M = zeros(m_max,p_max-1);
for p = 2:p_max
    for m = p*(p-1)-1 : m_max
        M(m,p-1) = alpha(p-1)/theta(m);
    end
end
