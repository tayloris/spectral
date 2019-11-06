
function SolveProblem_2_Fourier(N)


% u_exact = @(x)  exp(cos(2*pi*x))-exp(1);     %u_prima = - 2*pi*sin(*2*pix) exp(cos(2*pi*x))
% f_rhs = @(x) 2*pi*2*pi*exp(cos(2*pi*x)).*(-cos(2*pi*x)+sin(2*pi*x).^2);   


%t = pi*x + pi;   [-1, 1] ==>  [0, 2pi]

u_exact = @(t)  exp(cos(2*t))-exp(1);     %u_prima = - 2*sin(*2*t) exp(cos(2*t))
f_rhs = @(t) 2*2*exp(cos(2*t)).*(-cos(2*t)+sin(2*t).^2);   



figure('NumberTitle', 'off', 'Name', 'Fourier collocation method');


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

          uu = expand_fourier(xx,u_hat,ik);

          u = expand_fourier(x,u_hat,ik);

            
            exact = u_exact(xx);
            
            
            Error(N)=norm(uu-exact);
          %% Plotting solution
          
          
        Plot_Solution(x,u,xx,uu,exact)



end
  function u_expand = expand_fourier(x,u_hat , ik)
  N =length(u_hat);
  %u_expand(1)=0;
  for i= 1:length(x)
      u_expand(i) = 0;
        for m = 1:length(u_hat)     
            u_expand(i) = real( u_expand(i) + (1./N)*u_hat(m)*exp(ik(m)*x(i)));                           
        end      
  end
  %u_expand(length(x))=0;
  u_expand = u_expand - u_expand(1);
  end