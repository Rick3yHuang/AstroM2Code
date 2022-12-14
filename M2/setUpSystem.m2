
--restart
--recursionLimit = 100
--FF = ZZ/911

R = FF[-*p_(1,1),p_(1,2),p_(1,3),p_(1,4),p_(1,5),p_(2,1),p_(2,2),p_(2,3),p_(2,4),
    p_(2,5),p_(3,1),p_(3,2),p_(3,3),p_(3,4),p_(3,5),*-w_1,w_2,w_3,v_(1,1),
    v_(1,2),v_(1,3),v_(2,1),v_(2,2),v_(2,3),-*u_(1,1),u_(1,2),u_(1,3),u_(1,4),
    u_(1,5),u_(2,1),u_(2,2),u_(2,3),u_(2,4),u_(2,5),u_(3,1),u_(3,2),u_(3,3),
    u_(3,4),u_(3,5),*-o_(1,1),o_(1,2),o_(1,3),o_(1,4),o_(1,5),o_(2,1),o_(2,2),
    o_(2,3),o_(2,4),o_(2,5),c_1,c_2,c_3,c_4,c_5,c_6]

------------------------------ Small functions that will be used ---------------------------------
-- Cross Product
cross = method(Vector,Vector) := (u,v) -> (
    f := new MutableList from {3:0};
    f#0 = u_(1)*v_(2) - u_(2)*v_(1);
    f#1 = u_(2)*v_(0) - u_(0)*v_(2);
    f#2 = u_(0)*v_(1) - u_(1)*v_(0);
    out := vector(toList f);
    return out
    )

-- Pointwise Multiplication
ptwsProduct = method(Vector,Vector) := (u,v) -> (
    return diagonalMatrix(entries(u))*v
    )

-- Dot Product
dot = method(Vector,Vector) := (u,v) -> (
    return sum(entries(diagonalMatrix(entries(u))*v))
    )

------------------------------------------ Step 1 -------------------------------------------------

--- Inputs:
-- Direction matrix of observers (bearing) u(3x5)
-*u = matrix{{u_(1,1),u_(1,2),u_(1,3),u_(1,4),u_(1,5)},
    {u_(2,1),u_(2,2),u_(2,3),u_(2,4),u_(2,5)},
    {u_(3,1),u_(3,2),u_(3,3),u_(3,4),u_(3,5)}}*-
u = random(FF^3,FF^5)
u1 = vector u_{0}
-- Basis vectors w(3x1),v_1(3x1),v_2(3x1)
w = vector{w_1,w_2,w_3}
v_1 = vector{v_(1,1),v_(1,2),v_(1,3)}
v_2 = vector{v_(2,1),v_(2,2),v_(2,3)}

--- Detailed Procedure
-- Let a denotes the cross product of vector w and p, then we have the equation:
-- v = a/|a| => v_i = a_i/|a| => v_i·|a| = a_i => v_i^2·|a|^2 = a_i^2 => |a|^2·v_i^2-a_i^2 = 0
funcStep1 = method(Vector,Vector,Vector) := (w,p,v) -> (
    n := #entries(w);
    f := new MutableList from {n:0};
    a := cross(w,p);
    for i from 0 to n-1 do(
    	f#i = dot(a,a)*v_(i)^2-a_(i)^2
	);
    return f;
    )

--- Outputs: each one of f1 and f2 is a list of three polynomials
-- v_2 = (w x u1)/|w x u1|
f1 = funcStep1(w,u1,v_2)
-- v_1 = (w x v_2)/|w x v_2|
f2 = funcStep1(w,v_2,v_1)

------------------------------------------ Step 2 -------------------------------------------------

--- New Inputs:
-- Position matrix of observers p(3x5)
-*p = matrix{{p_(1,1),p_(1,2),p_(1,3),p_(1,4),p_(1,5)},
    {p_(2,1),p_(2,2),p_(2,3),p_(2,4),p_(2,5)},
    {p_(3,1),p_(3,2),p_(3,3),p_(3,4),p_(3,5)}}*-
p = random(FF^3,FF^5)
-- Observations (on the quadric) o(2x5)
o = matrix{{o_(1,1),o_(1,2),o_(1,3),o_(1,4),o_(1,5)},
    {o_(2,1),o_(2,2),o_(2,3),o_(2,4),o_(2,5)}}

--- Detailed Procedure
-- let pl denotes p^T·w, ul denotes u^T·w, and let r_i denotes the i^th obseration in the orginal
-- coordinates. Then we have:
--    r_i = p_i+(-pl_i/ul_i)·u_i and dot(v_j,r_i) = o_(j,i)
-- => dot(v_j,p_i+(-pl_i/ul_i)·u_i) = o(j,i)
-- => dot(v_j,ul_i·p_i-pl_i·u_i) = ul_i·o_(j,i)
-- => ul_i·o_(j,i)-dot(v_j,ul_i·p_i-pl_i·u_i) = 0
funcStep2 = method(Matrix,Matrix,Vector,ZZ) := (p,u,w,j) -> (
    n := #entries(transpose p);
    f := new MutableList from {n:0};
    pl := (transpose p)*w;
    ul := (transpose u)*w;
    for i from 0 to n-1 do(
    	f#i = ul_(i)*o_(j-1,i) - dot(v_(j),(ul_(i)*substitute(p_i,R)-pl_(i)*substitute(u_i,R)))
	);
    return f;
    )

--- Outputs: each one of f3 and f4 is a list of five polynomials
-- Compute the v_1-coordinates
f3 = funcStep2(p,u,w,1)
-- Compute the v_2-coordinates
f4 = funcStep2(p,u,w,2)

------------------------------------------ Step 3 -------------------------------------------------
--- New Inputs
-- Coefficient vector c(6x1)
c_6 = 1
c = vector{c_1,2*c_2,c_3,2*c_4,2*c_5,c_6}
-- tA(6x5) is the transpose of the data matrix A(5x6)
-- tA = [x^2,xy,y^2,x,y,1]
xsq = ptwsProduct(vector (transpose o)_{0}, vector (transpose o)_{0})
xy = ptwsProduct(vector (transpose o)_{0}, vector (transpose o)_{1})
ysq = ptwsProduct(vector (transpose o)_{1}, vector (transpose o)_{1})
one = vector{1,1,1,1,1}
ones = substitute(one,R)
tA = matrix{toList entries xsq,toList entries xy,toList entries ysq,
    toList entries (transpose o)_0,toList entries (transpose o)_1,
    toList entries ones}

--- Detailed Procedure
-- dot(tA_i,c)=0
funcStep3 = method(Matrix,Vector) := (tA,c) -> (
    n := #(transpose tA);
    f = new MutableList from (n:0);
    for i from 0 to n-1 do (
	f#i = dot(tA_i,c)
	);
    return f
    )

--- Output: f5 is a list of five polynomials
f5 = funcStep3(tA,c)

------------------------------------------ Step 4 -------------------------------------------------
--- New Input:
-- F_1,F_2 are two free terms when finding foci
--- Quadric matrix Q
Q = matrix {{c_(0),c_(1),c_(3)},{c_(1),c_(2),c_(4)},{c_(3),c_(4),c_(5)}}

--- Output: f6 and f7 are two polynomials
f6 = Q_(1,2)^2-Q_(1,1)*Q_(2,2)-Q_(0,2)^2+Q_(0,0)*Q_(2,2)
f7 = 2*(Q_(0,2)*Q_(1,2)-Q_(0,1)*Q_(2,2))
f8 = dot(w,w)-1
------------------------------------ Set up the system --------------------------------------------
f = flatten{toList f1,toList f2,toList f3,toList f4,toList f5,f6,f7,f8}
netList f
I = ideal(f)
end

restart
recursionLimit = 100
FF = ZZ/911
load "setUpSystem.m2"
gbTrace = 3
gbI = ideal groebnerBasis(I, Strategy=>"F4"); 
elapsedTime dim I
degree I

