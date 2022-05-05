function x = skewnormal(mu,sigma,alpha,dim)

x1 = random('Normal',0,1,dim);
x2 = random('Normal',0,1,dim);

x = mu + sigma/sqrt(1+alpha^2)*(alpha*abs(x1) + x2);

