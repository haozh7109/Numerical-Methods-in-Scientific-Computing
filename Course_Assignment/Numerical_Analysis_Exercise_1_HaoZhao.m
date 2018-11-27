%% ========================================================================
%Exercise 1 of Numerical Analysis Course, Coded by Hao Zhao, Feb 15,2016
%% ========================================================================


close all;

%% ========================================================================
%--------------------Problems and Computer Exercises : 4.1.1 (a) ----------
%% ========================================================================

for n = 2:2:16
 
    % (1) option 1: define the equidistant points
    x = zeros(1,n);
    f = zeros(1,n);
    for i = 1:n
        x(i) = -1 + 2*(i-1)/(n-1);
        f(i) = 1/(3+x(i));
    end
    
    % define the Vandermonde Matrix for the linear system V'c = F
    V = zeros(n,n);
    c = zeros(n,1);
    p = zeros(n,1); 
    for i = 1:n
        for j = 1:n
            V(i,j) = x(i)^(j-1);
        end
    end
    % solve the linear sytem
    c = V\f';
    % derive the numerical solution
    P = V * c;
    
    % derive the calculation difference
    diff = abs(f-P');
    
    % derive the the log difference
    diff_log = diff;
    diff_log(diff_log == 0) = 1e-20; 
    
    figure;
    subplot(2,2,1);semilogy(x,diff_log,'--o');
    axis([-1 1 1e-20 1e-15])
    xlabel('x')
    ylabel('log of   |f(x)- p(x)|')
%     title(sprintf('Interpolation Error of f(x) = 1/(3+x):  (equidistant points n = %d)',n));
    title(sprintf('Interpolation Error: |f(x)-p(x)| :  (equidistant points n = %d)',n));
    
    % (2) option 2: define the Chebyshev points
    x = zeros(1,n);
    f = zeros(1,n);
    for i = 1:n
        x(i) = cos((2*i-1)*pi/(2*n));
        f(i) = 1/(3+x(i));
    end
    
    % define the Vandermonde Matrix for the linear system V'c = F
    V = zeros(n,n);
    c = zeros(n,1);
    p = zeros(n,1); 
    for i = 1:n
        for j = 1:n
            V(i,j) = x(i)^(j-1);
        end
    end
    % solve the linear sytem
    c = V\f';
    % derive the numerical solution
    P = V * c;
    
    % derive the calculation difference
    diff = abs(f-P');
    
    % derive the the log difference
    diff_log = diff;
    diff_log(diff_log == 0) = 1e-20;
    
    subplot(2,2,2);semilogy(x,diff_log,'--o');
    axis([-1 1 1e-20 1e-15])
    xlabel('x')
    ylabel('log of   |f(x)- p(x)|')
%     title(sprintf('Interpolation Error of f(x) = 1/(3+x):  (Chebyshev points n = %d)',n));
    title(sprintf('Interpolation Error:  |f(x)-p(x)|: (Chebyshev points n = %d)',n));

   %% ----------------------------------------------------------------------------------------
    Pertb = 0.2;
    
    % (1) option 1: define the equidistant points
    x = zeros(1,n);
    f = zeros(1,n);
    for i = 1:n
        x(i) = -1 + 2*(i-1)/(n-1);
        f(i) = 1/(3+x(i));
    end
    
    % define the Vandermonde Matrix for the linear system V'c = F
    V = zeros(n,n);
    c = zeros(n,1);
    q = zeros(n,1); 
    for i = 1:n
        for j = 1:n
            V(i,j) = (x(i)- Pertb)^(j-1);
        end
    end
    % solve the linear sytem
    c = V\f';
    % derive the numerical solution
    q = V * c;
    
    % derive the calculation difference
    diff = abs(q'-P');
    
    % derive the the log difference
    diff_log = diff;
%     diff_log(diff_log == 0) = 1e-20; 
    

    subplot(2,2,3);semilogy(x,diff_log,'--o');
%     axis([-1 1 1e-20 1e-15])
    xlabel('x')
    ylabel('log of   |q(x)- p(x)|')
    title(sprintf('Interpolation Error of |q(x)-p(x)|:  (equidistant points n = %d with Random Perturbation))',n));
    
    % (2) option 2: define the Chebyshev points
    x = zeros(1,n);
    f = zeros(1,n);
    for i = 1:n
        x(i) = cos((2*i-1)*pi/(2*n));
        f(i) = 1/(3+x(i));
    end
    
    % define the Vandermonde Matrix for the linear system V'c = F
    V = zeros(n,n);
    c = zeros(n,1);
    q = zeros(n,1); 
    for i = 1:n
        for j = 1:n
            V(i,j) = (x(i)- Pertb)^(j-1);
        end
    end
    % solve the linear sytem
    c = V\f';
    % derive the numerical solution
    q = V * c;
    
    % derive the calculation difference
    diff = abs(q'-P');
    
    % derive the the log difference
    diff_log = diff;
    diff_log(diff_log == 0) = 1e-20;
    
    subplot(2,2,4);semilogy(x,diff_log,'--o');
    %     axis([-1 1 1e-20 1e-15])
    xlabel('x')
    ylabel('log of   |q(x)- p(x)|')
%     title(sprintf('Interpolation Error of f(x) = 1/(3+x):  (Chebyshev points n = %d with Random Perturbation)',n));
    title(sprintf('Interpolation Error:  |q(x)-p(x)|: (Chebyshev points n = %d  with Random Perturbation))',n));
    
end

%% ========================================================================
%--------------------Problems and Computer Exercises : 4.1.1 (b) ----------
%% ========================================================================

for n = 16:-8:8
    
    % (1) option 1: define the equidistant points
    x = zeros(1,n);
    f = zeros(1,n);
    for i = 1:n
        x(i) = -1 + 2*(i-1)/(n-1);
%         f(i) = 1/(x(i)^2 +1);
        f(i) = sin(x(i)*pi);
    end
    
    % define the Vandermonde Matrix for the linear system V'c = F
    V = zeros(n,n);
    p = zeros(n,1); 
    for i = 1:n
        for j = 1:n
            V(i,j) = x(i)^(j-1);
        end
    end
    % solve the linear sytem
    c = V\f';
    % derive the numerical solution
    P = V * c;
    
    % define the actual solution of f(x), the denser sample ratio used.        
    X  = -1: 0.05 :1;
    F  = sin(X.*pi);
    
    % plot the f(x) series.
    figure;
    subplot(1,2,1);plot(X,F,'b',x,P,'r--o');
    axis([-1 1 -1 1])
    xlabel('x')
    ylabel('f(x)')
    title(sprintf('Interpolation of random vector f :  (equidistant points n = %d)',n));
    
    
    % (2) option 2: define the Chebyshev points
    x = zeros(1,n);
    f = zeros(1,n);
    for i = 1:n
        x(i) = cos((2*i-1)*pi/(2*n));
        f(i) = sin(x(i)*pi);
    end
    % define the Vandermonde Matrix for the linear system V'c = F
    V = zeros(n,n);
    c = zeros(n,1);
    p = zeros(n,1); 
    for i = 1:n
        for j = 1:n
            V(i,j) = x(i)^(j-1);
        end
    end
    % solve the linear sytem
    c = V\f';
    % derive the numerical solution
    P = V * c;
    
    % define the actual solution of f(x), the denser sample ratio used.        
    X  = -1: 0.05 :1;
    F  = sin(X.*pi);
    
    % plot the f(x) series.
    subplot(1,2,2);plot(X,F,'b',x,P,'r--o');
    axis([-1 1 -1 1])
    xlabel('x')
    ylabel('f(x)')
    title(sprintf('Interpolation of random vector f :  (Chebyshev points n = %d)',n));

end
