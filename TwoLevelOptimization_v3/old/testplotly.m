x = 80 * randn(1, 30);
y = 80 * randn(size(x));
r = randi(1500, size(x));
c = randi(10, size(x));

fig = figure;

scatter(x, y, r, c, 'filled', 'MarkerEdgeColor', 'k')

%--PLOTLY--%

response = fig2plotly(fig, 'filename', 'matlab-bubble-chart',  'strip', false);
plotly_url = response.url;


% -----------------
[X1,Y1]= meshgrid(-5:.2:5,-5:.2:5);
syms x y
f=((x^2-1)+(y^2-4)+(x^2-1)*(y^2-4))/(x^2+y^2+1)^2
zfun = @(x, y) eval(vectorize(f))
Z1=zfun(X1,Y1);
plot3(X1,Y1,Z1)
fig2plotly()

