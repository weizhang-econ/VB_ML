function lden = ligampdf(x,a,b)
    lden = a.*log(b) - gammaln(a) - (a+1).*log(x) - b./x;
end