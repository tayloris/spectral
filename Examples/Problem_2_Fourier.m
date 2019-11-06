function Problem_2_Fourier(N)


% u_exact = @(x)  exp(cos(2*pi*x))-exp(1);     %u_prima = - 2*pi*sin(*2*pix) exp(cos(2*pi*x))
% f_rhs = @(x) 2*pi*2*pi*exp(cos(2*pi*x)).*(-cos(2*pi*x)+sin(2*pi*x).^2);   
warning('off','all')


%t = pi*x + pi;   [-1, 1] ==>  [0, 2pi]

u_exact = @(t)  exp(cos(2*t))-exp(1);     %u_prima = - 2*sin(*2*t) exp(cos(2*t))
f_rhs = @(t) 2*2*exp(cos(2*t)).*(-cos(2*t)+sin(2*t).^2);   




N_max =N;

figure('NumberTitle', 'off', 'Name', 'Fourier collocation method');
for N=2:2:N_max
        x = linspace(0,2*pi-2*pi/N,N); 


        f =f_rhs(x);


        ik = 1i*[0:N/2 -N/2+1:-1];   %i * wave number vector (matlab ordering) 

           ik2 = ik.*ik;                %multiplication factor for second derivative   
           ii  = find(ik~= 0);           %indices where ik is nonzero   
          ik2inverse = ik2;             %initialize zeros in same locations as in ik2  
          ik2inverse(ii) = 1./ik2(ii);  %multiplier factor to solve u'' = f


        f_hat = fft(f);  % Fourier Transform of f O(N log(N))
        
          u_hat = ik2inverse.*f_hat;

          u = real(ifft(u_hat));   %imaginary parts should be round off level

          xx = 0: 0.01:2*pi;

          u_expand_fine = expand_fourier(xx,u_hat,ik);

          u_expande_coarse = expand_fourier(x,u_hat,ik);

            Error(N)=norm(u_expand_fine-u_exact(xx));
            
            
          subplot(2,2,1); plot(x,u_expande_coarse,'bo',xx,u_expand_fine,'b',xx,u_exact(xx),'r')
          title('u(x)\approx u^N(x)');
          xlabel('x');
          ylabel('u');
            grid on
          legend('u^N(x_i)','u^N(x)','u(x)')
end

subplot(2,2,2); loglog(1:N,real(Error),'b*',1:N,(1:N).^(-1),1:N,(1:N).^(-2))
legend('Fourier O(N^{-N})','O(N^{-1})','O(N^{-2})');
xlabel('N');
title('Convergence');
grid on
  

subplot(2,2,3);  plot(xx, abs(u_exact(xx)-u_expand_fine))
xlabel('x');
ylabel('abs(u-u^N)')
title('Error');
grid on


subplot(2,2,4);  loglog(1:N, (1:N).*log(1:N),1:N,1:N,1:N,(1:N).^2)
legend('O(N log(N)) ','O(N^{1})','O(N^{2})');
xlabel('N');
title('Computaitonal cost');
grid on

end
 % norm(u2-v, inf)
  
  function u_expand = expand_fourier(x,u_hat , ik)
  N =length(u_hat);
  %u_expand(1)=0;
  for i= 1:length(x)
      u_expand(i) = 0;
        for m = 1:length(u_hat)     
            u_expand(i) =  real( u_expand(i) + (1./N)*u_hat(m)*exp(ik(m)*x(i)));                           
        end      
  end
  %u_expand(length(x))=0;
  u_expand = u_expand - u_expand(1);
  end