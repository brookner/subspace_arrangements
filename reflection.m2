F=ZZ/32003;
c=random(1,20)/random(1,10);
R=F[x,y];
Q = j -> {{x^j-(-1)^j*((c+1)*x+c*y)^j}, {y^j-(-1)^j*((c+1)*x+c*y)^j}};
P = j -> (c+1)*x^j+c*y^j+(-1)^j*((c+1)*x+c*y)^j;
compare = (H1, H2) -> ((numerator H1) * (value denominator H2) == (numerator H2) * (value denominator H1))

-- Global variable
f = 0;
S = 0;
mat = 0;

presentP = (N) -> (
    S = F[z_1..z_N, Degrees=>{2..N+1}];
    f = map(R, S, apply(N, i->P(i+2)));
    ker f )

generateP = () -> (
    count = 2;
    H := reduceHilbert hilbertSeries presentP(1);
    H' := reduceHilbert hilbertSeries presentP(2);
    while not compare(H, H') do (
	H = H';
	H' = reduceHilbert hilbertSeries presentP(count+1);
	count = count+1;
	);
    << "Hilbert series for R: " << (toString H) << endl;
    count - 1
    )

--m=6;
--mat=matrix(Q(2));
--for i from 3 to m do mat=(mat|(matrix(Q(i))));
--g=map(R^2, S^(-toList(2..m)), f, mat);
--reduceHilbert hilbertSeries coimage g

presentM = (N) -> (
    mat = matrix(Q(2));
    for i from 3 to N do mat=(mat|matrix(Q(i)));
    g = map(R^2, S^(-toList(2..N)), f, mat);
    coimage g )

generateM = () -> (
    count = 3;
    H := reduceHilbert hilbertSeries presentP(2);
    H' := reduceHilbert hilbertSeries presentP(3);
    while not compare(H, H') do (
	H = H';
	count = count+1;
	H' = reduceHilbert hilbertSeries presentP(count);
	);
    << "Hilbert series for M: " << (toString H) << endl;
    count - 1
    )
