function J = jacobian(model,data,delta, xobs, zobs, x, b, zt, rho)

% model here is a 1D array 

lm = length(model);
ld = length(data);

J  = zeros(ld,lm); % for Java or C - take care of the right order
                   % of for-loops (inner loop over the fast dimension)

for im=1:lm
  model(im) = model(im) + delta;
  if im>1
    model(im-1) = model(im-1) - delta; % to undo the perturbation
  end                                  % created at the previous step
   
  
  datap = basinResponse(xobs, zobs, x, b, zt, model, rho); 
                                       % calculating the perturbed data
                                       
                                       
  for id=1:ld
    J(id,im) = (datap(id)-data(id))/delta;
  end
end
