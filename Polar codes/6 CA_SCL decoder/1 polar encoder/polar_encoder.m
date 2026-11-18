function v = polar_encoder(c, N, Q) 
Q = Q(Q<=N);
n = log2(N);
K = length(c);
u = zeros(1, N);
free = Q(end-K+1:end); 
frozen = Q(1:end-K);
    
u(free) = c;
u(frozen) = 0;

% Polar Transform
m = 1; % number of bits combined
for d = n-1:-1:0
    for i =1:2*m:N
        a = u(i:i+m-1);     % first part
        b = u(i+m:i+2*m-1); % second part
        u(i:i+2*m-1)=[mod(a+b, 2), b]; 
    end
    m = m*2; 
end
v = u;


