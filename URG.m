function [U,R,G] = URG(A_1, A0, A1, n)

U = cell(1, n+1);
R = cell(1, n+1);
G = cell(1, n+1);

U{n+1} = A0{n+2};
R{n+1} = A1{n+1}*inv(eye(n+2)-U{n+1});
G{n+1} = inv(eye(n+2)-U{n+1})*A_1{n+1};

for k = n:-1:1
    U{k} = A0{k+1} + A1{k+1}*G{k+1};
    R{k} = A1{k}*inv(eye(k+1)-U{k});
    G{k} = inv(eye(k+1)-U{k})*A_1{k};
end