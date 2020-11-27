function [A_1, A0, A1] = Construct_P(lambda, mu, p, q, n)

% initialization
A_1 = cell(1, n+1);
A0 = cell(1, n+2);
A1 = cell(1, n+1);

% construction of A_{-1}
Temp_A_1 = zeros(n+2,n+1);
for k = 2:n+2
    Temp_A_1(k,k-1) = (mu*q)/(lambda+mu);
end
A_1{n+1} = Temp_A_1;

for k = n:-1:1
    Temp_A_1 = [eye(k+1), zeros(k+1,n+1-k)] * A_1{n+1} * [eye(k); zeros(n+1-k, k)];
    A_1{k} = Temp_A_1;
end

% construction of A_{0}
for k = 1:n-1
    Temp_A0 = zeros(k);
    Temp_A0(1,end) = (mu*(1-q))/(lambda+mu);
    for l = 2:k
        Temp_A0(l,l-1) = (mu*(1-q))/(lambda+mu);
    end
    A0{k} = Temp_A0;
end

A0{n} = eye(n)*lambda*(1-p)/(lambda+mu) + [zeros(1,n-1), (mu*(1-q))/(lambda+mu); eye(n-1)*(mu*(1-q))/(lambda+mu) zeros(n-1,1)];

for k = n:n+1
    Temp_A0 = eye(k+1)*lambda/(lambda+mu) + [zeros(1,k), (mu*(1-q))/(lambda+mu); eye(k)*(mu*(1-q))/(lambda+mu) zeros(k,1)];
    A0{k+1} = Temp_A0;
end

% construction of A_{1}
for k = 1:n-1
    Temp_A1 = [eye(k) * lambda/(lambda+mu) zeros(k,1)];
    A1{k} = Temp_A1;
end

A1{n} = [eye(n) * lambda*(p)/(lambda+mu), zeros(n,1)];
A1{n+1} = zeros(n+1,n+2);