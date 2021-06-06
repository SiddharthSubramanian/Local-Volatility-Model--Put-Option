%% Implicit Scheme to solve Black Scholes PDE to price a European PUT option


% Strike Price
K1 = 120;
K2 = 80;
% Interest rate
r = 0.03;
% Dividend yield
q = 0.05;
% Volatility (Constant)
a = 0.3;

% Set the minimal and maximal stock prices
% for K1 = 120
Smin1 = 0;
Smax1 = 4*K1;

% for K2 = 80
Smin2 = 0;
Smax2 = 4*K2;


% Setup our grid in stock price direction
N = 200;
% for K1 = 120
S1 = linspace(Smin1,Smax1,N+1)';
dS1 = S1(2) - S1(1);

% for K2 = 80
S2 = linspace(Smin2,Smax2,N+1)';
dS2 = S2(2) - S2(1);

% S stores all the prices except boundary points
% for K1 = 120
S_1 = S1(2:N);

% for K2 = 80
S_2 = S2(2:N);

% Setup our grid in time direction
T = 1;
M = 2000;
tau = linspace(0,T,M+1);
dtau = tau(2) - tau(1);


sigma1 = zeros(N-1,M+1);  % for K1 = 120
sigma2 = zeros(N-1,M+1);  % for K2 = 80

%Solve BS PDE (IMPLICIT)

Vold_imp1 = max(K1 - S_1,0);  % for K1 = 120

Vold_imp2 = max(K2 - S_2,0); % for K2 = 80

% Vector to store the option prices at time k+1
Vnew_imp1 = Vold_imp1;   % for K1 = 120

Vnew_imp2 = Vold_imp2;   % for K2 = 80


% for K1 = 120
alpha1 =zeros(N-1,M+1);
beta1 =zeros(N-1,M+1);


% for K2 = 80
alpha2 =zeros(N-1,M+1);
beta2 =zeros(N-1,M+1);

for k = 1:M+1
%     For K1 = 120 
    sigma1(:,k) = 0.25*exp(-tau(k))*(100./S_1).^a;

     for m = 1:N-1
        alpha1(m,k) = 0.5 * sigma1(m,k).^2.* S_1(m).^2 * dtau / (dS1^2);
        beta1(m,k) = (r - q) * S_1(m) * dtau / (2*dS1);
    end

    l1 = -alpha1(:,k) + beta1(:,k);
    d1 = 1 + r*dtau + 2*alpha1(:,k);
    u1 = - (alpha1(:,k) + beta1(:,k));
    
    
    d1(1) = d1(1) + 2*l1(1);
    u1(1) = u1(1) - l1(1);

    d1(N-1) = d1(N-1) + 2*u1(N-1);
    l1(N-1) = l1(N-1) - u1(N-1) ;

        
    Vnew_imp1= tridiag(d1,u1,l1,Vold_imp1);
    Vold_imp1 = Vnew_imp1   ;
    
%     For K2 = 80

    sigma2(:,k) = 0.25*exp(-tau(k))*(100./S_2).^a;
    for m = 1:N-1
        alpha2(m,k) = 0.5 * sigma2(m,k).^2.* S_2(m).^2 * dtau / (dS2^2);
        beta2(m,k) = (r - q) * S_2(m) * dtau / (2*dS2);
    end
    
    l2 = -alpha2(:,k) + beta2(:,k);
    d2 = 1 + r*dtau + 2*alpha2(:,k);
    u2 = - (alpha2(:,k) + beta2(:,k));
    
    d2(1) = d2(1) + 2*l2(1);
    u2(1) = u2(1) - l2(1);

    d2(N-1) = d2(N-1) + 2*u2(N-1);
    l2(N-1) = l2(N-1) - u2(N-1) ;

    Vnew_imp2= tridiag(d2,u2,l2,Vold_imp2);
    Vold_imp2 = Vnew_imp2   ;
end

% Interpolation to find the put price when S0 = 1.0
S0 = 100.0;
put_fdm_imp1 = interp1(S_1,Vold_imp1,S0);
put_fdm_imp2 = interp1(S_2,Vold_imp2,S0);


fprintf("\n Implicit Method European Put Price \n")

fprintf("\n Price at K = 120 \n")
disp(put_fdm_imp1)

fprintf("\n Price at K = 80 \n")
disp(put_fdm_imp2)


