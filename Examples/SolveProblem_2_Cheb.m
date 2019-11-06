function SolveProblem_2_Cheb(N)

%p13.m - solve linear BVP     u_xx = exp(4x) , u(-1) = u(1)= 0
warning('off','all')




u_exact = @(x)  exp(cos(2*pi*x))-exp(1);   %u_prima = - 2*pi*sin(*2*pix) exp(cos(2*pi*x))
f_rhs = @(x) 2*pi*2*pi*exp(cos(2*pi*x)).*(-cos(2*pi*x)+sin(2*pi*x).^2);   



N_max =N;


figure('NumberTitle', 'off', 'Name', 'Chebyshev collocation method');


[D,x]=cheb(N); % Diferentiation matrix

D2 = D^2;    % No eficient O(N^3)

D2=D2(2:N,2:N);%remove elements related to the boundar

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
