FUNMV - computes trigonometric and hyperbolic matrix functions times vectors

   [C,S,s,m,mv,mvd] = FUNMV(t,A,B,FLAG,[],PREC) computes f(t*A)*B or
   f(t*SQRTM(A))B, where f is a trigonometric or hyperbolic matrix function
   and SQRTM denotes any matrix square root, without
   explicitly forming f(t*A) or SQRTM(A). PREC is the required accuracy,
   'double', 'single' or 'half', and defaults to CLASS(A).
   
   Details on the underlying algorithms can be found in
   
   A. H. Al-Mohy, "A New Algorithm for Computing the Actions of Trigonometric and Hyperbolic Matrix Functions", https://arxiv.org/pdf/1708.08360.pdf
