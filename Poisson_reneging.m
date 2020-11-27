function W = Poisson_2(A_1, A0, A1, U, R, G, g, n)
n = n-1;
g_div = cell(1,n+2);
for ind = 1:n+2
    g_div{ind} = g(ind*(ind-1)/2+1:ind*(ind+1)/2);
end
y = cell(1,n+2);
y{1,1} = 0;
for r = 1:n+1
    sum_tmp = zeros(r+1,1);
    for k = r:n
        R_tmp = 1;
        for l = r+1:k+1
            R_tmp = R_tmp*R{l};
        end
        sum_tmp = sum_tmp + R_tmp*g_div{k+2};
    end
    y_tmp = inv(eye(r+1)-U{r})*(g_div{r+1}+sum_tmp)+G{r}*y{1,r};
    y{1,r+1} = y_tmp;
end

y{1} = g_div{1}+A1{1}*y{2};

% calculation of expected waiting time
W = cell(1,n+2);
W{1} = y{1}/(1-(A0{1}+A1{1}*G{1}));
for k = 1:n+1
    G_tmp = 1;
    for l = k:-1:1
        G_tmp = G_tmp*G{l};
    end
    W{k+1} = y{k+1}+G_tmp*W{1};
end
