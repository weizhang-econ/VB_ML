% this function simulate a standard SV model with 
% random walk and varianc omegah2

function h = genSVRW(h1,omegah2,T)

h = cumsum([h1; sqrt(omegah2)*randn(T-1,1)]);

end