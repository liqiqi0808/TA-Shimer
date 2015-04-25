% PS4 Q1: Convex Adjustment Costs: Numerical
% Qi
% April 24 2015
clear all
close all
clc
%% 1. Calculate Steady State Capital under Deterministic Demand
u_grid=[1;2];
pi=0.1;
R=1.03;
delta=0.05;

K_determ=log(u_grid)-log(R*delta+(R-1)*delta^2)

%% 2. Solve for Value Function and Policy Function under Stochastic Demand
% use discretization method
dif=99;
tol=1e-6;
K_lb=K_determ(1);
K_ub=K_determ(2);
N_K=1000;
K_grid=linspace(K_lb,K_ub,N_K);
V=R*(delta+delta^2)*K_determ*ones(1,N_K);
V_new=V;
K=ones(size(V));
P=[1-pi,pi;pi,1-pi];

% value function interation
while dif>tol
    for i_u=1:length(u_grid)
        for i_K=1:N_K
            s=u_grid*(1-exp(-K_grid(i_K)));
            gamma=K_grid/K_grid(i_K);
            c=delta*gamma./(1+delta*(1-gamma));
            [V_new(i_u,i_K),K(i_u,i_K)]=max(s(i_u)*ones(1,N_K)-c.*K_grid+P(i_u,:)*V/R);
        end
    end
    dif=max(max(abs(V_new-V)));
    V=V_new;
    % dif % monitor speed of convergence
end

% use projection method

grow_rate=K_grid(K)./(ones(size(u_grid))*K_grid);

figure(1)
plot(K_grid,grow_rate(1,:),K_grid,grow_rate(2,:))
hold on
plot([K_lb;K_ub],[1;1],'Color','y')
hold off
title('Capital/Output Growth Rate (Gross)')
xlabel('K')
legend('Low Demand','High Demand','Location','NorthEast')

%% 3. Transition from Low Steady State to High Steady State
figure(2)
plot(K_grid,K_grid(K(1,:)),K_grid,K_grid(K(2,:)))
hold on
plot([K_lb;K_ub],[K_lb;K_ub],'Color','y')
hold off
title('Policy Function')
xlabel('K')
legend('Low Demand','High Demand','Location','NorthWest')

T=10;
K_transit=zeros(T,1);
[~,K_lss]=max(K(1,:)==1:N_K);
K_transit(1)=K_lss;
K_H=K(2,:);
gamma_transit=zeros(T,1);
for t=2:T
    K_transit(t)=K_H(K_transit(t-1));
    gamma_transit(t-1)= K_grid(K_transit(t))/ K_grid(K_transit(t-1));
    
end
c_transit=delta*gamma_transit./(1+delta*(1-gamma_transit));
figure(3)
subplot(2,1,1)
plot(1:T,K_grid(K_transit))
title('Transition from Low Steady State to High Steady State')
xlabel('time')
ylabel('K')
subplot(2,1,2)
plot(1:T-1,c_transit(1:T-1))
xlabel('time')
ylabel('c')
%% 4. Calculate Variance and Autocorrelation of Investments
% simulate the path of demand shock
burnin=500;
T=500;
u_simu=zeros(burnin+T,1);
u_simu(1)=1;
s = RandStream('mt19937ar','Seed',123);
RandStream.setGlobalStream(s);
for t=2:burnin+T
    if rand()<pi
        u_simu(t)=3-u_simu(t-1);
    else
        u_simu(t)=u_simu(t-1);
    end
end
% simulate the path of K and c
K_simu=zeros(burnin+T,1);
gamma_simu=zeros(burnin+T,1);
K_simu(1)=K_lss;
for t=2:burnin+T
    K_u=K(u_simu(t),:);
    K_simu(t)=K_u(K_simu(t-1));
    gamma_simu(t-1)= K_grid(K_simu(t))/ K_grid(K_simu(t-1));
end
c_simu=delta*gamma_simu./(1+delta*(1-gamma_simu));
figure(4)
subplot(3,1,1)
plot(burnin+1:burnin+T-1,gamma_simu(burnin+1:burnin+T-1))
title('Time-Series of Simulated gamma, c and c*K')
xlabel('time')
ylabel('gamma')
subplot(3,1,2)
plot(burnin+1:burnin+T-1,c_simu(burnin+1:burnin+T-1))
xlabel('time')
ylabel('c')
subplot(3,1,3)
plot(burnin+1:burnin+T-1,c_simu(burnin+1:burnin+T-1).*K_simu(burnin+1:burnin+T-1))
xlabel('time')
ylabel('c*K')

% calculate variance and autocorrelation
% for investment rate gamma-1
var_gamma=var(gamma_simu(burnin+1:burnin+T-1))
autocor_gamma=autocorr(gamma_simu(burnin+1:burnin+T-1));
% for investment expenditure per capital c
var_c=var(c_simu(burnin+1:burnin+T-1))
autocor_c=autocorr(c_simu(burnin+1:burnin+T-1));
% for total investment expenditure c*K
var_cK=var(c_simu(burnin+1:burnin+T-1).*K_simu(burnin+1:burnin+T-1))
autocor_cK=autocorr(c_simu(burnin+1:burnin+T-1).*K_simu(burnin+1:burnin+T-1));
figure(5)
lag=10;
subplot(3,1,1)
plot(1:lag,autocor_gamma(1:lag))
title('Autocorrelations of Simulated gamma, c and c*K')
xlabel('lag')
ylabel('gamma')
subplot(3,1,2)
plot(1:lag,autocor_c(1:lag))
xlabel('lag')
ylabel('c')
subplot(3,1,3)
plot(1:lag,autocor_cK(1:lag))
xlabel('lag')
ylabel('c*K')