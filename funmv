function [C,S,s,m,mv,mvd,unA] = funmv(t,A,b,flag,M,prec,shift,bal,full_term)
       
%FUNMV   trigonometric and hyperbolic matrix functions times vectors
%   [C,S,s,m,mv,mvd] = FUNMV(t,A,B,FLAG,[],PREC) computes f(t*A)*B or
%   f(t*SQRTM(A))B, where f is a trigonometric or hyperbolic matrix function
%   and SQRTM denotes any matrix square root, without
%   explicitly forming f(t*A) or SQRTM(A). PREC is the required accuracy,
%   'double', 'single' or 'half', and defaults to CLASS(A).
%   A total of mv products with A or A^* are used, of which mvd are
%   for norm estimation.
%   The full syntax is
%   [C,S,s,m,mv,mvd,unA] = funmv(t,A,b,flag,M,prec,shift,bal,full_term).
%   unA = 1 if the alpha_p were used instead of norm(A).
%   If repeated invocation of FUNMV is required for several values of t
%   or B, it is recommended to provide M as an external parameter as
%   M = SELECT_TAYLOR_DEGREE_TRIG(A,b,flag,m_max,p_max,prec,shift,bal,true).
%   This also allows choosing different m_max and p_max.

%   The outputs C and S depend on the input FLAG as follows:

%   FLAG = 'cos.sin'        , C = cos(tA)B  and S = sin(tA)B
%   FLAG = 'cosh.sinh'      , C = cosh(tA)B and S = sinh(tA)B
%   FLAG = 'cos.sinc'       , C = cos(tA)B  and S = sinc(tA)B
%   FLAG = 'cosh.sinch'     , C = cosh(tA)B and S = sinch(tA)B
%   FLAG = 'cos.sinc.sqrt'  , C = cos(t*sqrt(A))B and S = sinc(t*sqrt(A))B
%   FLAG = 'cosh.sinch.sqrt', C = cosh(t*sqrt(A))B and S = sinch(t*sqrt(A))B

%   The parameter SHIFT is optinal only if FLAG = 'cos.sin'.

%   Example: 
%   To evaluate a combination of the form 
%                   cos(t*sqrt(A))*y0 + t*sinc(t*sqrt(A))*y1,
%   it is C(:,1) + t*S(:,2), where [C,S] = funmv(t,A,[y0,y1],'cos.sinc.sqrt') 
%   

%   Reference: A. H. Al-Mohy, A New Algorithm for Computing the Action 
%   of Trigonometric and Hyperbolic Matrix Functions

%   Awad H. Al-Mohy, January 18, 2018.

if nargin < 4
    error('funmv:NoInput', 'Not enough input parameters to funmv.')
end

if ~ismember(flag, {'cos.sin', 'cosh.sinh', 'cos.sinc'...
                    'cosh.sinch', 'cos.sinc.sqrt', 'cosh.sinch.sqrt'})
    error('funmv:NoInput', 'Choose a correct input of the parameter FLAG.')
end 
if nargin < 7 || isempty(shift), shift = true; end
flag1 = isequal( flag, 'cos.sin');
flag2 = isequal( flag, 'cosh.sinh');
if flag1                            , flag = 1; sign = 1; end
if flag2                            , flag = 1; sign = 0; end
if isequal( flag, 'cos.sinc')       , flag = 1; sign = 1; end
if isequal( flag, 'cosh.sinch')     , flag = 1; sign = 0; end
if isequal( flag, 'cos.sinc.sqrt')  , flag = 0; sign = 1; end
if isequal( flag, 'cosh.sinch.sqrt'), flag = 0; sign = 0; end

if ~(flag1 || flag2), shift = false; end
if nargin < 9  || isempty(full_term), full_term = false; end
if nargin < 8  || isempty(bal), bal = false; end
if bal
    [D,B] = balance(A);
    if norm(B,1) < norm(A,1), A = B; b = D\b; else bal = false; end
end
n = length(A); n0 = size(b,2);
pp = n0;
if flag, pp = 2*n0; end  % the number of matrix-vector prod., A(AB).  
if shift && (flag1 || flag2)
    mu = trace(A)/n;
    mu = full(mu);     % The code could be slower without the full!
    tmu = t*mu;
    A = A - mu*speye(n);     
end
if nargin < 6 || isempty(prec), prec = class(A); end
t2 = t;
if ~flag, t2 = t^2; end
if nargin < 5 || isempty(M)
   tt = 1;
   [M,mvd,~,unA] = select_taylor_deg_trig(t2*A,b,flag,[],[],prec,false,false,false);
   mv = mvd;
else
   tt = t; mv = 0; mvd = 0; unA = 0;
end
switch prec
    case 'double', tol = 2^(-53);
    case 'single', tol = 2^(-24);
    case 'half',   tol = 2^(-10);
end
s = 1;
if t == 0
    m = 0;
else
    [m_max,p] = size(M);
     S = diag(1:m_max);
     C = ( (ceil(abs(tt)*M))'*S );
     C (C == 0) = inf;
     if p > 1
         [cost, m] = min(min(C)); % cost is the overall cost.
     else
         [cost, m] = min(C);  % when C is one column. Happens if p_max = 2.
     end
     if cost == inf; cost = 0; end
     s = max(cost/m,1);            
end
undo_inside = false; undo_outside = false; 
% undo shifting inside or outside the loop
if shift
    if flag1 && ~isreal(tmu)
        cosmu = cos(tmu/s); sinmu = sin(tmu/s);   undo_inside = true;
    elseif flag1 && isreal(tmu) && abs(tmu)> 0
        cosmu = cos(tmu); sinmu = sin(tmu);       undo_outside = true;
    elseif flag2 &&  abs(real(tmu)) > 0
        cosmu = cosh(tmu/s); sinmu = sinh(tmu/s); undo_inside = true;
    elseif flag2 && ~real(tmu) && abs(tmu)> 0
        cosmu = cosh(tmu); sinmu = sinh(tmu);     undo_outside = true;
    end
end
mods = mod(s,2);
C0 = zeros(n,n0);
if mods, C0 = b/2; end
S = C0;
C1 = b;
for i = 1:s+1
    if i == s+1
        S = 2*S; C1 = S;         
    end
    V = C1;
    if undo_inside, Z = C1; end
    b = C1;  
    c1 = norm(b,inf);
    for k = 1:m            
        even = 2*k; 
        if i <= s
            odd = even-1; q = 1/(even+1); 
        else               % when i = s+1, compute Taylor poly. of sinc
            odd = even+1; q = odd; 
        end
                           
       if flag, b = A*b; end        
        b = (A*b)*((t/s)^2/(even*odd));           
        mv = mv + pp;
        V = V + ((-1)^(k*sign))*b;
        if undo_inside, Z = Z + (((-1)^(k*sign))*q)*b; end
            
        c2 = norm(b,inf);
        if ~full_term            
            if c1 + c2 <= tol*norm(V,inf) 
                break
            end
            c1 = c2;
        end
    end    
    if undo_inside
        if i <= s
            V = V*cosmu + A*(Z*(((-1)^sign)*t*sinmu/s));
            mv = mv + n0;
        else         
            V = A*(V*(t*cosmu/s)) + Z*sinmu;
            mv = mv + n0;
        end
    end
        
    if i == 1 
        C2 = V;
    elseif i <= s 
        C2 = 2*V - C0; 
    end  
    if i <= s-1 && xor( mods,mod(i,2) ) 
        S = S + C2;   % sum of C_i for even i's if s is odd and vice versa     
    end
    C0 = C1; C1 = C2;    
end
C = C2;
if undo_inside
    S = V;
elseif flag1 || flag2
    S = A*(V*(t/s)); 
    mv = mv + n0; 
else
    S = V/s;   % sinc(tA)B (or sinch(tA)B ) if flag = 1 and
               % sinc(t*sqrt(A))B (or sinch(t*sqrt(A))B) if flag = 0
end
if undo_outside
    C = cosmu*C + (((-1)^sign)*sinmu)*S; 
    S = sinmu*C2 + cosmu*S;
end

if bal, C = D*C; S = D*S; end
