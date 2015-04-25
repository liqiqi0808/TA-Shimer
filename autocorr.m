function r=autocorr(y)

Z = (y - mean(y))/std(y);
N = length(Z);
r = zeros(N-1,1);

    for k=0:N-2
        r(k+1) = sum(Z(1:end-k).*Z(1+k:end))/(N-k-1);
    end

end
