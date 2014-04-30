function WdtWd = dataW(sigma)

n     =     length(sigma);
d     = 1./(sigma.*sigma);

WdtWd = spdiags(d, 0, n, n);

end


