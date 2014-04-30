function Ap = getAp(J, beta, WmtWm, WdtWd, rhs)

step1 = beta*WmtWm*rhs;

step2 = J*rhs;

step3 = WdtWd*(step2);

step4 =    J'*(step3);

Ap    = step1 + step4;

end

