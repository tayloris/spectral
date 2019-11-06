%p13.m - solve linear BVP     u_xx = exp(4x) , u(-1) = u(1)= 0

function Problem_1_Cheb(N)

warning('off','all')


 f_rhs = @(x) exp(4*x);
 u_exact = @(x) (exp(4*x)-sinh(4)*x-cosh(4))/16;


%u_exact = @(x)  exp(cos(2*pi*x))-exp(1);   %u_prima = - 2*pi*sin(*2*pix) exp(cos(2*pi*x))
%f_rhs = @(x) 2*pi*2*pi*exp(cos(2*pi*x)).*(-cos(2*pi*x)+sin(2*pi*x).^2);   



N_max =N;


figure('NumberTitle', 'off', 'Name', 'Chebyshev collocation method');

for N=1:N_max

[D,x]=cheb(N); % Diferentiation matrix

D2 = D^2;    % No eficient O(N^3)

D2=D2(2:N,2:N);%remove elements related to the boundary

f = f_rhs(x(2:N)); 

condition_number(N)=cond(D2);
u=D2\f ; %Poisson eq.solved here 


u=[0;u;0]; % Imposing boundary condition


xx=-1:.01:1;
uu=polyval(polyfit(x,u,N),xx); %Spentral expansion. 
exact = u_exact(xx);

Error(N)=norm(uu-exact);


subplot(2,2,1); plot(x,u,'o', xx , uu,'b',xx,exact) %,'linewidth',.8)
xlabel('x');
legend('u^N(x_i)','u^N(x)','u(x)')
title('u(x)\approx u^N(x)');
ylabel('u');
grid on 



end
% 
N=N_max;
subplot(2,2,2); loglog(1:N,Error,1:N,(1:N).^(-1),1:N,(1:N).^(-2))
legend('Cheb O(N^{-N})','O(N^{-1})','O(N^{-2})');
xlabel('N');
title('Convergence');
grid on
% 
% 
subplot(2,2,3);  plot(xx, abs(u_exact(xx)-uu))
xlabel('x');
ylabel('abs(u-u^N)')
title('Error');
grid on

subplot(2,2,4);  loglog(1:N,condition_number,1:N,1:N,1:N,(1:N).^2)
legend('Cheb','O(N^{1})','O(N^{2})');
xlabel('N');
title('Condition Number');
grid on

end