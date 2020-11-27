function W = Expected_Revenue_reneging(Reward, n,lambda,mu,p,q)
%tic
N = (n+1)*(n+2)/2;
g = ones(N,1)*(-1/(lambda+mu));
for ind = 1:n+1
    g(ind*(ind-1)/2+1) = g(ind*(ind-1)/2+1) + mu*q*Reward/(lambda+mu);
end
[A_1, A0, A1] = Construct_P_reneging(lambda, mu, p, q, n);
[U,R,G] = URG_reneging(A_1, A0, A1, n);
W = Poisson_reneging(A_1, A0, A1, U, R, G, g, n);
%toc