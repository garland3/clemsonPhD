x = rand(100,100);
d = ones(1,10)
x = diag(d)
x=x+flip(x)
x(3,:)=1;


[t, numChanged] = CheckForConerElements(x, size(x,2), size(x,1), 0.3)