function W = Expected_Waiting_Time(n,lambda,mu,p,q)
%tic
g = 1/(lambda+mu);
[A_1, A0, A1] = Construct_P(lambda, mu, p, q, n);
[U,R,G] = URG(A_1, A0, A1, n);
W = Poisson(A_1, A0, A1, U, R, G, g, n);
%toc