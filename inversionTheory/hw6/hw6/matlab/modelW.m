function WmtWm = modelW(nm, dx, alphaS, alphaX)
% Here we assume that our model is uniformely sampled
% Therefore, we can just pass the sampling interval (dx)
% which is a number here

WstWs = speye(nm,nm)*dx;

% WxtWx is a tridiagonal matrix
%
% |  1 -1  0  0 |   m=4 in this example;
% | -1  2 -1  0 |
% |  0 -1  2  0 |
% |  0  0 -1  1 |
Wx    = spdiags([ones(nm-1,1), -ones(nm-1,1)], [1, 0], nm-1, nm);


WmtWm = alphaS*WstWs+alphaX*(Wx'*Wx)/dx; % It is not necessary to divide
                                         % by dx at the end


end

