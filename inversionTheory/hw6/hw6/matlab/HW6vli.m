%% Loading data

clear all;
close all;

load dat.txt;

xobs  = dat(:,1);
zobs  = dat(:,2);
dobs  = dat(:,3);
sigma = dat(:,4);

rho   = -1.0;

N     = length(xobs);

%% Discretizing the basin into a set of contiguous vertical 2D prizms

b      =  100; % the width of a single vertical 2D prizm

xleft  = -6000;
xright =  6000;
x      =  (xleft-50:b:xright-50)'; % since x is a coordinate of the center
                                   % of the dyke, we subtract 50
M      = length(x);

zt     = zeros(M,1);

%% Choosing an initial model
zb     = zeros(M,1) + 800;         % Initial guess, z of the bottom = 200
                                   % z-axis points downwards

%% Calculating the gravity response of the discretized basin 
                                   %(initial guess)
dpre   = basinResponse(xobs, zobs, x, b, zt, zb, rho);
                                   % this is an initial guess data
%% Calculating the initial Jacobian
delta  = 25;                       % the amount of perturbation
J      = jacobian(zb,dpre,delta, xobs, zobs, x, b, zt, rho);  
                                   % zb=initial model, dpre=initial data
                                   
%% Calculating weighting matrices (Wd, Ws, Wx)
alphaS = 0.00002; % for the model function term
alphaX = 1.0    ; % for the first derivative of the model function
dx     = b      ;

WdtWd  =  dataW(sigma);
WmtWm  = modelW(M, dx, alphaS, alphaX);
                                   
%% Solving the non-linear inverse problem

% Defining the range of beta
%B      = logspace(-8,-4,30);
B = zeros(1,1);
B(1)=0.1;
lB     = length(B);

DATA   = zeros(length(xobs),lB);     % to store the inversion 
MODEL  = zeros(length(zb)  ,lB);     % results for different betas
MISFIT = zeros(lB, 1);               % Data misfit
MNORM  = zeros(lB, 1);               % Model norm


zsol   = zb;           % this is what we are inverting for

bool   = false;        % As long as data misfit > 2%, bool = false
                       % this is a convergence criterion

imax   = 25;           % Tentative max # of iterations
                       % Not necessary, for not to change 
                       % the size of the arrays that store the history

DMOD  = zeros(N,imax);
ZBOT  = zeros(M,imax);
ALPHA = zeros(1,imax);
 
for k=1:lB                                         % TIKHONOV BETA FOR-LOOP

  iterations = 0;   
  
  %while bool==false                             % MODEL UPDATING WHILE-LOOP
  for i=1:10
    iterations = iterations +1;
    
    
    % Calculating data misfit
    dmod       = basinResponse(xobs, zobs, x, b, zt, zsol, rho);
    misfit0    = dobs - dmod; % data misfit before model updating

    
    % Calculating Jacobian for this iteration   
    J          = jacobian(zb,dmod,delta, xobs, zobs, x, b, zt, rho);
                     % notice that the Jacobian changes at every iteration
                     % of model updating
                         
    
    % Obtaining the search direction by using CG method
    rhs        = J'*WdtWd*misfit0-B(k)*WmtWm*zsol;
    %IsSPD     = all(eig(J'*J) > 0)
    p          = CG_solver(J, B(k), WmtWm, WdtWd, rhs);
                     
    % Objective function before model update
    obj0       = misfit0'*WdtWd*misfit0 + B(k)*zsol'*WmtWm*zsol;
                 
    obj_ls     = 2*obj0;          % to make sure that while-condition
                                  % is not satisfied                 
                                  
    alpha      = 1;               % take a full step for the first try
    
    
    while obj_ls >= obj0                           % LINE SEARCH WHILE-LOOP
        zsol_new_ls = zsol + alpha*p;
        dmod_ls       = basinResponse(xobs, zobs, x, b, zt, zsol_new_ls, rho);
        %dmod_ls     = grav(zsol_new_ls);
        
        misfit_ls   = dobs - dmod_ls;
        
        obj_ls      = misfit_ls'          * WdtWd * misfit_ls + ...
                      B(k) * zsol_new_ls' * WmtWm * zsol_new_ls;
        
        if obj_ls > obj0
          alpha = alpha/2;
        end
        
    end % END OF LINE SEARCH WHILE-LOOP
    disp(alpha)
    ALPHA(iterations) = alpha;    % This is just to store a history
                                  % of step lengths for all iterations

    % Updating the model
    zsol    = zsol + alpha*p;
    dmod1   = basinResponse(xobs, zobs, x, b, zt, zsol, rho);

    
    misfit1 = dobs - dmod1;
    
    obj1    = misfit1'  * WdtWd * misfit1 + ...
              B(k)*zsol'* WmtWm * zsol;
    
    if (obj0-obj1)/obj0 <= 2.0
      disp('Convergence criterion has been met');
      bool = true;      
    end
    
  DMOD(:, iterations) = dmod1;
  ZBOT(:, iterations) = zsol;  
                                  
    
  end % END OF THE MODEL UPDATING WHILE-LOOP

DATA  (:,k) = dmod1;
MODEL (:,k) = zsol;
MISFIT  (k) = misfit1'* WdtWd * misfit1; 
MNORM   (k) = zsol'   * WmtWm * zsol;  
  
end % END OF THE TIKHONOV BETA FOR-LOOP 


%% Plotting
figure(1);                % observed and predicted data 
                          % (predicted with the initial guess) 
plot(xobs,dobs);hold on;
plot(xobs,dpre);

figure(2);                % Jacobian
pcolor(J);
