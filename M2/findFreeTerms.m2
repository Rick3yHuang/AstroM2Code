restart
R = ZZ/101[x,y,Q11,Q12,Q13,Q21,Q22,Q23,Q31,Q32,Q33];
f = (Q23^2-Q22*Q33-Q12^2+Q11*Q33)+(-Q12^2+Q11*Q22)*(y^2-x^2)+2*(Q12*Q23-Q13*Q22)*x-2*(Q12*Q13-Q11*Q23)*y;
g = 2*(-Q12^2+Q11*Q22)*x*y+2*(Q13*Q23-Q12*Q33)-2*(Q12*Q23-Q13*Q22)*y-2*(Q12*Q13-Q11*Q23)*x;
h = eliminate(y,ideal(f,g));
toString(h)
toString(sub(h,x=>0))
