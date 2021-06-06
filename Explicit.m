%% Explicit Scheme to solve Black Scholes PDE to price a European PUT option

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


% for K1 = 120
alpha1 =zeros(N-1,M+1);
beta1 =zeros(N-1,M+1);


% for K2 = 80
alpha2 =zeros(N-1,M+1);
beta2 =zeros(N-1,M+1);



Vold1 = max(K1 - S_1,0);   % for K1 = 120
Vold2 = max(K2 - S_2,0);   % for K2 = 80

% Vector to store the option prices at time k+1
Vnew1 = Vold1;
Vnew2 = Vold2;


for k = 1:M+1
    %     For K1 = 120 
    sigma1(:,k) = 0.25*exp(-tau(k))*(100./S_1).^a;
    
    for m = 1:N-1
        alpha1(m,k) = 0.5 * sigma1(m,k).^2.* S_1(m).^2 * dtau / (dS1^2);
        beta1(m,k) = (r - q) * S_1(m) * dtau / (2*dS1);
    end
    
    
    l1 = alpha1(:,k) - beta1(:,k);
    d1 = 1 - r*dtau - 2*alpha1(:,k);
    u1 = alpha1(:,k) + beta1(:,k);

    
    % Boundary condition for European call option
    boundary1 = [l1(1)*(2*Vold1(1)-Vold1(2));zeros(N-3,1);u1(N-1)*(2*Vold1(N-1)-Vold1(N-2))];

    % Explicit iteration scheme
    for j=1:N-1
        if(j==1)
            Vnew1(j) = d1(j)*Vold1(j) + u1(j)*Vold1(j+1);
        elseif(j<N-1)
            Vnew1(j) = l1(j)*Vold1(j-1) + d1(j)*Vold1(j) + u1(j)*Vold1(j+1);
        else
            Vnew1(j) = l1(j)*Vold1(j-1) + d1(j)*Vold1(j);
        end
    end

    % Update the vectors from time k to time k+1
    Vold1 = Vnew1 + boundary1 ; 
    
    
    %     For K2 = 80 
    sigma2(:,k) = 0.25*exp(-tau(k))*(100./S_2).^a;
    
    for m = 1:N-1
        alpha2(m,k) = 0.5 * sigma2(m,k).^2.* S_2(m).^2 * dtau / (dS2^2);
        beta2(m,k) = (r - q) * S_2(m) * dtau / (2*dS2);
    end
    
    
    l2 = alpha2(:,k) - beta2(:,k);
    d2 = 1 - r*dtau - 2*alpha2(:,k);
    u2 = alpha2(:,k) + beta2(:,k);

    
    % Boundary condition for European call option
    boundary2 = [l2(1)*(2*Vold2(1)-Vold2(2));zeros(N-3,1);u2(N-1)*(2*Vold2(N-1)-Vold2(N-2))];

    % Explicit iteration scheme
    for j=1:N-1
        if(j==1)
            Vnew2(j) = d2(j)*Vold2(j) + u2(j)*Vold2(j+1);
        elseif(j<N-1)
            Vnew2(j) = l2(j)*Vold2(j-1) + d2(j)*Vold2(j) + u2(j)*Vold2(j+1);
        else
            Vnew2(j) = l2(j)*Vold2(j-1) + d2(j)*Vold2(j);
        end
    end

    % Update the vectors from time k to time k+1
    Vold2 = Vnew2 + boundary2 ; 

end


S0 = 100.0;
put_fdm_exp1 = interp1(S_1,Vold1,S0);

put_fdm_exp2 = interp1(S_2,Vold2,S0);

fprintf("\n Explicit Method European Put Price \n")

fprintf("\n Price at K = 120 \n")
disp(put_fdm_exp1)

fprintf("\n Price at K = 80 \n")
disp(put_fdm_exp2)


% [call_bs,put_bs] =   blsprice(S,K,r,T,sigma(:,5001),q);
% 
% % Plot of FDM compared with Built in BS function
% subplot(2,1,1)
% plot(S,Vold,S,put_bs,'LineWidth',2)
% title('European Put price, Implicit - BS formula')
% xlabel('Stock price')
% ylabel('Put price')
% legend('Implicit','BS formula','Location','SouthEast')
% 
% % Plot of Diference of FDM compared with Built in BS function
% subplot(2,1,2)
% plot(S,Vold - put_bs, 'LineWidth',2)
% title('Difference of European Put prices, Implicit - BS formula')
% xlabel('Stock price')
% ylabel('Differences in price')
