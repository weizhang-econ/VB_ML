T = 5000;
mu = [1,2,3]';
sig2 = [.01, .05, .1]';
phi = [.9 .95 .96]';
h = zeros(T,n);
for t = 1:T
    if t == 1
        h(t,:) = mu' + sqrt(sig2'./(1-phi'.^2)).*randn(1,n); 
    else        
        h(t,:) = mu' + phi'.*(h(t-1,:)-mu') + sqrt(sig2').*randn(1,n);
    end
end