F=ZZ/32003;
--F=QQ;
R=F[x,y];
-- All code for the family (a+1,a,1)
initialize = (a) -> (
    Q = j -> {{x^j-(-1)^j*((a+1)*x+a*y)^j}, {y^j-(-1)^j*((a+1)*x+a*y)^j}};
    P = j -> (a+1)*x^j+a*y^j+(-1)^j*((a+1)*x+a*y)^j;
    )

a=random(1,20)/random(1,10);
initialize(a);
compare = (H1, H2) -> ((numerator H1) * (value denominator H2) == (numerator H2) * (value denominator H1))

-- Global variable
f = 0;
S = 0;
mat = 0;

presentR = (N) -> (
    S = F[z_1..z_N, Degrees=>{2..N+1}];
    f = map(R, S, apply(N, i->
	    P(i+2)));
    ker f )

generateR = () -> (
    count = 2;
    H := reduceHilbert hilbertSeries presentR(1);
    H' := reduceHilbert hilbertSeries presentR(2);
    while not compare(H, H') do (
	H = H';
	H' = reduceHilbert hilbertSeries presentR(count+1);
	count = count+1;
	);
    << "Hilbert series for R: " << (toString H) << endl;
    count - 1
    )

presentM = (N) -> (
    mat = matrix(Q(1));
    for i from 2 to N do mat=(mat|matrix(Q(i)));
    g = map(R^2, S^(-toList(1..N)), f, mat);
    coimage g )

generateM = () -> (
    count = 3;
    H := reduceHilbert hilbertSeries presentM(2);
    H' := reduceHilbert hilbertSeries presentM(3);
    while not compare(H, H') do (
	H = H';
	count = count+1;
	H' = reduceHilbert hilbertSeries presentM(count);
	);
    << "Hilbert series for M: " << (toString H) << endl;
    count - 1
    )
