FUNMV   trigonometric and hyperbolic matrix functions times vectors

   [C,S,s,m,mv,mvd] = FUNMV(t,A,B,FLAG,[],PREC) computes f(t*A)*B or
   f(t*SQRTM(A))B, where f is a trigonometric or hyperbolic matrix function
   and SQRTM denotes any matrix square root, without
   explicitly forming f(t*A) or SQRTM(A). PREC is the required accuracy,
   'double', 'single' or 'half', and defaults to CLASS(A).
