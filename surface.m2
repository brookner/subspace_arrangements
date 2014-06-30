--use notation from cmpart.pdf
F=ZZ/32003;
R=F[x,y,z];

--compare reduced Hilbert series H1,H2
compare = (H1,H2) -> (numerator H1 == numerator H2) and (value denominator H1 == value denominator H2)

--output: presentation of surface whose coordinate ring is generated by 
--P_s with parameters p,q,r with s=1,...,N

surf = (p, q, r, N) -> (
    S := F[v_1..v_N, Degrees=>toList(1..N)];
    a := p/r;
    b := q/r;
    eqns := apply(N, i-> a * x^(i+1) + b * y^(i+1) + (-1)^(i+1)*(a*x + b*y)^(i+1));
    f := map(R, S, eqns);
    ker f
    )

-- Computes for partition (q+r, q, r, r)
surf4 = (p, q, r, s, N) -> (
    S := F[v_1..v_N, Degrees=>toList(1..N)];
    a := p/s;
    b := q/s;
    c := r/s;
    
    eqns := apply(N, i-> a * x^(i+1) + b * y^(i+1) + c * z^(i+1) + (-1)^(i+1) * (a*x + b*y + c*z)^(i+1));
    f := map(R, S, eqns);
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
    count := 2;
    H := reduceHilbert hilbertSeries surf4(p, q, r, s, 1);
    H' := reduceHilbert hilbertSeries surf4(p, q, r, s, 2);
    while not compare(H, H') do (
	H = H';
	count = count+1;
	H' = reduceHilbert hilbertSeries surf4(p, q, r, s, count);
	);
    << "guess4 for (" << p << "," << q << "," << r << "," << s << ") " << count << "\n";
    count - 1
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
    g := guess4(p, q, r, s);
    I := surf4(p, q, r, s, g);
    length res I == codim I
    )

printMe = () -> (
    for i from 1 to 4 do (
	for j from 1 to i do (
	    for k from 1 to j do (
		"helloworld.txt" << "isCM(4," << i << "," << j << "," << k << "): " << isCM4(4,i,j,k) << endl; )));
    "helloworld.txt" << close;)
