n=8; 
F=ZZ/32003;
R=F[x_1..x_n]/ideal(sum(apply(n,i->x_(i+1))));
mu={x_1=>x_2,x_2=>x_1};
nu=append(apply(n-1,i->(x_(i+1)=>x_(i+2))), x_n=>x_1);

--creates the ideal x_1=...=x_(lambda_1), 
--x_(lambda_1+1)=...=x_(lambda_2), ..., 
--which is the subspace E_lambda from cmpart.pdf
base = lambda -> (
    I:={};
    counter:=0;
    for i in lambda do (
	counter=counter+1;
	anchor:=counter;
	for j from 0 to i-2 do (
	    counter=counter+1;
	    I=append(I,x_anchor-x_counter);
	    );
	);
    ideal mingens ideal I
    )

--input: set L of ideals
--output: closure under S_n of this set of ideals
orbit = L -> (
    Lnew = L;
    flag = true;
    while flag do (
     	temp = {};
     	for i in Lnew do (
	    iM = ideal(mingens(sub(i,mu)));
	    iN = ideal(mingens(sub(i,nu)));
	    if not(member( iM, L)) then temp = append(temp, iM);
	    if not(member( iN, L)) then temp = append(temp, iN);
     	    );
     	L = join(L,temp);
     	if #temp==0 then flag = false;
     	Lnew=temp;
     	);
    unique(L)
    )

--output: ideal of subspace arrangement X_lambda
test = lambda -> (
    intersect orbit({base(lambda)})
    )

end
restart
load "part.m2"
---------------------
betti res test({3,2})
       0  1  2 3
total: 1 10 15 6
    0: 1  .  . .
    1: .  .  . .
    2: . 10 15 6

I = test({6}); betti res I, codim(I)
I = test({5, 1}); betti res I, codim(I)
I = test({4, 2}); betti res I, codim(I)
I = test({4, 1, 1}); betti res I, codim(I)
I = test({3, 3}); betti res I, codim(I)
I = test({3, 2, 1}); betti res I, codim(I)
        0  1  2  3 4
(total: 1 25 45 22 1, 3)
     0: 1  .  .  . .
     1: .  .  .  . .
     2: .  .  .  . .
     3: .  .  .  . .
     4: .  .  .  . .
     5: . 15 14  1 .
     6: . 10 31 20 .
     7: .  .  .  1 .
     8: .  .  .  . 1
I = test({3, 1, 1, 1}); betti res I, codim(I)
        0 1 2 3
(total: 1 5 5 1, 2)
     0: 1 . . .
     1: . . . .
     2: . . . .
     3: . . . .
     4: . . . .
     5: . 5 . .
     6: . . 5 .
     7: . . . 1
I = test({2, 2, 2}); betti res I, codim(I)
I = test({2, 2, 1, 1}); betti res I, codim(I)

I = test({6, 1}); betti res I, codim(I)
I = test({5, 2}); betti res I, codim(I)
I = test({5, 1, 1}); betti res I, codim(I)
I = test({4, 3}); betti res I, codim(I)
I = test({4, 2, 1}); betti res I, codim(I)
        0  1   2   3  4 5
(total: 1 56 140 120 36 1, 4)
     0: 1  .   .   .  . .
     1: .  .   .   .  . .
     2: .  .   .   .  . .
     3: .  .   .   .  . .
     4: . 21  35  15  . .
     5: . 20  50  34  1 .
     6: . 15  55  71 35 .
     7: .  .   .   .  . .
     8: .  .   .   .  . 1
I = test({4, 1, 1, 1}); betti res I, codim(I)
I = test({3, 3, 1}); betti res I, codim(I)
I = test({3, 2, 2}); betti res I, codim(I)
        0  1  2   3  4  5
(total: 1 21 69 105 70 14, 4)
     0: 1  .  .   .  .  .
     1: .  .  .   .  .  .
     2: .  .  .   .  .  .
     3: .  .  .   .  .  .
     4: . 21 14   .  .  .
     5: .  . 20   .  .  .
     6: .  . 35 105 70 14
I = test({3, 2, 1, 1}); betti res I, codim(I) --too hard
I = test({3, 1, 1, 1, 1}); betti res I, codim(I)
        0  1  2 3
(total: 1 14 14 1, 2)
     0: 1  .  . .
     1: .  .  . .
     2: .  .  . .
     3: .  .  . .
     4: .  .  . .
     5: .  .  . .
     6: .  .  . .
     7: .  .  . .
     8: . 14 14 .
     9: .  .  . .
    10: .  .  . .
    11: .  .  . 1
I = test({2, 2, 2, 1}); betti res I, codim(I)
I = test({2, 2, 1, 1, 1}); betti res I, codim(I) --too hard

I = test({8}); betti res I, codim(I)
I = test({7, 1}); betti res I, codim(I)
I = test({6, 2}); betti res I, codim(I)
I = test({6, 1, 1}); betti res I, codim(I)
I = test({5, 3}); betti res I, codim(I)
I = test({5, 2, 1}); betti res I, codim(I)
I = test({5, 1, 1, 1}); betti res I, codim(I)
I = test({4, 4}); betti res I, codim(I)
I = test({4, 3, 1}); betti res I, codim(I)
I = test({4, 2, 2}); betti res I, codim(I)
I = test({4, 2, 1, 1}); betti res I, codim(I)
I = test({4, 1, 1, 1, 1}); betti res I, codim(I)
I = test({3, 3, 2}); betti res I, codim(I)
I = test({3, 3, 1, 1}); betti res I, codim(I)
I = test({3, 2, 2, 1}); betti res I, codim(I)
I = test({3, 2, 1, 1, 1}); betti res I, codim(I)
I = test({3, 1, 1, 1, 1, 1}); betti res I, codim(I)
I = test({2, 2, 2, 2}); betti res I, codim(I)
I = test({2, 2, 2, 1, 1}); betti res I, codim(I)
I = test({2, 2, 1, 1, 1, 1}); betti res I, codim(I)
I = test({2, 1, 1, 1, 1, 1, 1}); betti res I, codim(I)
I = test({1, 1, 1, 1, 1, 1, 1, 1}); betti res I, codim(I)
