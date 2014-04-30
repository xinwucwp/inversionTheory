function gz = basinResponse(xobs, zobs, x, b, zt, zb, rho)

k  = length(xobs);  % the # of observations = 
                    % the length of the output array gz

m  = length(zb)  ;  % the number of prizms in a basin 

gz = zeros(k,1)  ;

for i= 1:m
  gz = gz + vdyke(xobs,zobs,x(i),b,zt(i),zb(i),rho);
end

end
