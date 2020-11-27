n = 4;
lambda = 1;
mu = 0.2;
q = 0.1;

W1 = Expected_Waiting_Time(n,lambda,mu,0.5,q);

W2 = Expected_Waiting_Time(n,lambda,mu,0.7,q);

(W2{n+1}(end)-W1{n+1}(end)) - (W2{n}(end)-W1{n}(end))

(W2{n}(end)) - (W1{n+1}(end))
