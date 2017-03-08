D =  transpose(reshape(sym('d', [1 9]), [3,3]))
epsilon = transpose(sym('e', [1 3]))

Dsymetric = [D(1,1) D(1,2) D(1,3);
            D(1,2) D(2,2) D(2,3);
            D(1,3) D(2,3) D(3,3)]

equation =transpose(epsilon)*Dsymetric*epsilon % negative because we are trying to maximize = minimize the negative

lambda=sym('lambda')
constraint = ((epsilon(1) +epsilon(2) +epsilon(3))-1)


L = equation-lambda*constraint

dLde1 = diff(L, epsilon(1))
dLde2 = diff(L, epsilon(2))
dLde3 = diff(L, epsilon(3))

%temp={dLde1==0,dLde2==0,dLde3==0,constraint==0}
    
answer = solve([dLde1==0,dLde2==0,dLde3==0,constraint==0],[epsilon(1),epsilon(2),epsilon(3),lambda])
simplify(answer.e1)
answer.e2
answer.e3
answer.lambda

% Sub in the lambda value
% dLde1 = subs(dLde1,lambda
% answer = solve([dLde1==0,dLde2==0,dLde3==0],[epsilon(1),epsilon(2),epsilon(3),lambda])