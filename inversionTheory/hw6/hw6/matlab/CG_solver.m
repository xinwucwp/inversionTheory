function dh = CG_solver(J, beta, WmtWm, WdtWd, rhs)

[~, nm] = size(J);

dh0      = zeros(nm,1);   % initial guess is a zero-model
dh       = dh0; 


r        = rhs;           % initial residual for the chosen
                          % initial zero-model
                                   
p        = r;             % initial search direction for the CG
                          % (it is not what we are solving for!)
                          
for k=0:nm                % since the maximum # of iterations
                          % is the # of unknowns

  if (norm(r/rhs) <= eps*10000) % convergence criterion, 10^(-12)
                                % for double precision
                                % eps = built-in Matlab constant
    break;
  end
                          
  saved   = r'*r;
  
  Ap      = getAp(J, beta, WmtWm, WdtWd, rhs);
  
  alphaCG = saved/(p'*Ap); 
  
  dh      = dh + alphaCG*p;
  r       = r  - alphaCG*Ap;
  
  betaCG  = r'*r/saved;
  
  p       = r + betaCG*p;
  
end

end