function SolveProblem_1_Cheb(N)

%p13.m - solve linear BVP     u_xx = exp(4x) , u(-1) = u(1)= 0
warning('off','all')


 f_rhs = @(x) exp(4*x);
 u_exact = @(x) (exp(4*x)-sinh(4)*x-cosh(4))/16;


N_max =N;


figure('NumberTitle', 'off', 'Name', 'Chebyshev collocation method');


[D,x]=cheb(N); % Diferentiation matrix

D2 = D^2;    % No eficient O(N^3)

D2=D2(2:N,2:N);%remove elements related to the boundary

f = f_rhs(x(2:N)); 

condition_number(N)=cond(D2);
u=D2\f ; %Poisson eq.solved here 


u=[0;u;0]; % Imposing boundary condition

%% Spectral Expantion
xx=-1:.01:1;
uu=polyval(polyfit(x,u,N),xx); %Spentral Expantion. 
exact = u_exact(xx);

Error(N)=norm(uu-exact);


%% Plotting 
Plot_Solution(x,u,xx,uu,exact)

end