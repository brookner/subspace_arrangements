--TODO: in all isCM functions, we show pdim == codim. Make sure codim == 0 and then just replace that with 0 everywhere.

--use notation from cmpart.pdf
F=ZZ/32003;
R=F[w,x,y,z];

--global variables to analyze runtime and where the bottlenecks of the program are
globalCount = 1;
globalHilb;
state = "Initializing";

--compare reduced Hilbert series H1,H2
--TODO: make sure compareNew is strictly good and deprecate compare
compare = (H1,H2) -> (numerator H1 == numerator H2) and (value denominator H1 == value denominator H2)
compareNew = (H1, H2) -> ((numerator H1) * (value denominator H2) == (numerator H2) * (value denominator H1))


--output: presentation of surface whose coordinate ring is generated by 
--P_s with parameters p,q,r with s=1,...,N. The P_s are normalized
--such that P_1=0, reducing the number of variables by one.

surf = (p, q, r, N) -> (
    S := F[v_1..v_N, Degrees=>toList(1..N)];
    a := p/r;
    b := q/r;
    eqns := apply(N, i-> a * x^(i+1) + b * y^(i+1) + (-1)^(i+1)*(a*x + b*y)^(i+1));
    f := map(R, S, eqns);
    ker f
    )

-- output: same as above, with an extra element in the partition

surf4 = (p, q, r, s, N) -> (
    S := F[v_1..v_N, Degrees=>toList(1..N)];
    a := p/s;
    b := q/s;
    c := r/s;
    
    eqns := apply(N, i-> a * x^(i+1) + b * y^(i+1) + c * z^(i+1) + (-1)^(i+1) * (a*x + b*y + c*z)^(i+1));
    f := map(R, S, eqns);
    state = "Calculating kernel for presentation";
    ker f
    )

surf5 = (p, q, r, s, t, N) -> (
    S := F[v_1..v_N, Degrees=>toList(1..N)];
    a := p/t;
    b := q/t;
    c := r/t;
    d := s/t;
    
    eqns := apply(N, i-> a * x^(i+1) + b * y^(i+1) + c * z^(i+1) + d * w^(i+1) + (-1)^(i+1) * (a*x + b*y + c*z + d*w)^(i+1));
    f := map(R, S, eqns);
    state = "Calculating kernel for presentation";
    ker f
    )

surfGen = (N) -> (
    S := F[v_1..v_N, Degrees=>toList(1..N)];
    
    eqns := apply(N, i-> (A+1)*x^(i+1) + A*y^(i+1) + z^(i+1) + (-1)^(i+1) * ((A+1)*x + A*y + z)^(i+1));
    f := map(R, S, eqns);
    state = "Calculating kernel for presentation";
    ker f
    )

--output: smallest N such that P_(N+1) is in the subalgebra generated by P_1,...,P_N

guess = (p, q, r) -> (
    count := 2;
    H := reduceHilbert hilbertSeries surf(p, q, r, 1);
    H' := reduceHilbert hilbertSeries surf(p, q, r, 2);
    while not compare(H, H') do (
	H = H';
	count = count+1;
	H' = reduceHilbert hilbertSeries surf(p, q, r, count);
	);
    count - 1
    )

guess4 = (p, q, r, s) -> (
    count := 3;
    H := reduceHilbert hilbertSeries surf4(p, q, r, s, 1);
    H' := reduceHilbert hilbertSeries surf4(p, q, r, s, 2);
    H'' := reduceHilbert hilbertSeries surf4(p, q, r, s, 3);
    while not (compare(H, H') and compare(H, H'')) do (
	H = H';
	H' = H'';
	count = count+1;
	globalCount = count;
	H'' = reduceHilbert hilbertSeries surf4(p, q, r, s, count);
	);
    globalHilb = H;
    count - 2
    )

guess5 = (p, q, r, s, t) -> (
    count := 2;
    H := reduceHilbert hilbertSeries surf5(p, q, r, s, t, 1);
    H' := reduceHilbert hilbertSeries surf5(p, q, r, s, t, 2);
    H'' := reduceHilbert hilbertSeries surf5(p, q, r, s, t, 3);
    while not (compare(H, H') and compare(H, H'')) do (
	H = H';
	H' = H'';
	count = count+1;
	globalCount = count;
	H'' = reduceHilbert hilbertSeries surf5(p, q, r, s, t, count);
	);
    << "guess4 for (" << p << "," << q << "," << r << "," << s << "," << t << ") " << (count-2) << endl;
    globalHilb = H;
    count - 2
    )

-- Computing CM-ness of (a+1, a, 1, 1) for a as a variable

guessGen = () -> (
    count := 2;
    H := reduceHilbert hilbertSeries surfGen(1);
    H' := reduceHilbert hilbertSeries surfGen(2);
    H'' := reduceHilbert hilbertSeries surfGen(3);
    while not (compare(H, H') and compare(H, H'')) do (
	H = H';
	H' = H'';
	count = count+1;
	globalCount = count;
	H'' = reduceHilbert hilbertSeries surfGen(count);
	);
    << "guessGen is " << (count - 2) << end;
    globalHilb = H;
    count - 2
    )

--assuming guess produces meaningful results, 
-- i.e., if P_(N+1) being generated by P_1,...,P_N implies the same for all P_j with j>N, 
-- then this returns true iff surface defined by p,q,r is CM

isCM = (p, q, r) -> (
    g := guess(p, q, r);
    I := surf(p, q, r, g);
    length res I == codim I
    )

isCM4 = (p, q, r, s) -> (
    state = "Guessing";
    g := guess4(p, q, r, s);
    state = "Generating surface presentation";
    I := surf4(p, q, r, s, g);
    state = "Calculating pdim";
    length res I == codim I
    )

isCM5 = (p, q, r, s, t) -> (
    state = "Guessing";
    g := guess5(p, q, r, s, t);
    state = "Generating surface presentation";
    I := surf5(p, q, r, s, t, g);
    state = "Calculating pdim";
    length res I == codim I
    )

printMe = () -> (
    for i from 1 to 5 do (
    	time result = isCM4(6,i,i,i);
	"helloworld.txt" << "isCM(6," << i << "," << i << "," << i << "): " << result << endl; );
    "helloworld.txt" << close;)
