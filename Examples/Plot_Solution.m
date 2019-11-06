%% Plotting
function Plot_Solution(x,u,xx,uu,exact)


subplot(2,1,1); plot(x,u,'o', xx , uu,'b',xx,exact) %,'linewidth',.8)
xlabel('x');
legend('u^N(x_i)','u^N(x)','u(x)')
title('u(x) \approx u^N(x)');
ylabel('u');
grid on 


% 
subplot(2,1,2);  plot(xx, abs(exact-uu))
xlabel('x');
ylabel('abs(u-u^N)')
title('Error');
grid on
end