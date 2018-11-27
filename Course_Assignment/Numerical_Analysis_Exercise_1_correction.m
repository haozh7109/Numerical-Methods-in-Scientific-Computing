%% ========================================================================
%Exercise 1 of Numerical Analysis Course, Coded by Hao Zhao, Feb 15,2016
%% ========================================================================

close all;

%========================================================================
%% --------------------Problems and Computer Exercises : 4.1.1 (a-1) ------
%========================================================================

for n = 2:2:16
 
    %% (1) option 1: define the equidistant points
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
    
    % define the actual solution of f(x), using denser sample (1000) used.
    m = 1000;
    x      = -1: 2/m :1;
    f_act  = 1 ./(3 + x);
    
    % define the numerical solution of f(x), using denser sample (1000) used.
    V_num  = zeros(size(x,2),size(c,1));
    for i = 1:size(x,2)
        for j = 1:size(c,1)
            V_num(i,j) = x(i)^(j-1);
        end
    end
    f_num  = V_num * c;
    
    % plot the f(x) series.
    figure;subplot(2,2,1);plot(x,f_act,'b',x,f_num,'r');
    axis([-1 1 -1 1])
    xlabel('x')
    ylabel('f(x)')
    title(sprintf('Interpolation of f(x) by p(x) :  (equidistant points n = %d)',n));
    
    % derive the calculation difference
    diff = abs(f_act-f_num');

    % plot the error.
    subplot(2,2,2);semilogy(x,diff);
    xlabel('x')
    ylabel('log of   |f(x)- p(x)|')
    title(sprintf('Interpolation Error: |f(x)-p(x)| :  (equidistant points n = %d)',n));
    
    %% (2) option 2: define the Chebyshev points
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
    
    % define the actual solution of f(x), using denser sample (1000) used.
    m = 1000;
    x      = -1: 2/m :1;
    f_act  = 1 ./(3 + x);
    
    % define the numerical solution of f(x), using denser sample (1000) used.
    V_num  = zeros(size(x,2),size(c,1));
    for i = 1:size(x,2)
        for j = 1:size(c,1)
            V_num(i,j) = x(i)^(j-1);
        end
    end
    f_num  = V_num * c;
    
    % plot the f(x) series.
    subplot(2,2,3);plot(x,f_act,'b',x,f_num,'r');
    axis([-1 1 -1 1])
    xlabel('x')
    ylabel('f(x)')
    title(sprintf('Interpolation of f(x) by p(x) :  (Chebyshev points n = %d)',n));
    
    % derive the calculation difference
    diff = abs(f_act-f_num');

    % plot the error.
    subplot(2,2,4);semilogy(x,diff);
    xlabel('x')
    ylabel('log of   |f(x)- p(x)|')
    title(sprintf('Interpolation Error: |f(x)-p(x)| :  (Chebyshev points n = %d)',n));
end

%========================================================================
%% --------------------Problems and Computer Exercises : 4.1.1 (a-2) ------
%========================================================================

for n = 16
    
    %% (1) option 1: define the equidistant points
    x = zeros(1,n);
    f = zeros(1,n);
    for i = 1:n
        x(i) = -1 + 2*(i-1)/(n-1);
        f(i) = 1/(3+x(i));
    end
    
    % add random perturbation for f(x)
    f_pert = f;
    f_pert(round(size(f,2)/2)) = f_pert(round(size(f,2)/2)) * 1.00001;
    
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
    c_pert = V\f_pert';
    
    
    % define the actual solution of f(x), using denser sample (1000) used.
    m = 1000;
    x      = -1: 2/m :1;
    f_act  = 1 ./(3 + x);
    
    % define the numerical solution of f(x), using denser sample (1000) used.
    V_num  = zeros(size(x,2),size(c,1));
    for i = 1:size(x,2)
        for j = 1:size(c,1)
            V_num(i,j) = x(i)^(j-1);
        end
    end

    f_num  = V_num * c;
    f_num_pert  = V_num * c_pert;
    
    % derive the calculation difference
    diff = abs(f_num_pert'-f_num');

    % plot the error.
    figure;subplot(2,1,1);semilogy(x,diff);
    xlabel('x')
    ylabel('log of   |q(x)- p(x)|')
    title(sprintf('Interpolation Error after single point small perturbation: |q(x)-p(x)| :  (equidistant points n = %d)',n));
     
    %% (2) option 2: define the Chebyshev points
    x = zeros(1,n);
    f = zeros(1,n);
    for i = 1:n
        x(i) = cos((2*i-1)*pi/(2*n));
        f(i) = 1/(3+x(i));
    end
    
    % add random perturbation for f(x)
    f_pert = f;
    f_pert(round(size(f,2)/2)) = f_pert(round(size(f,2)/2)) * 1.00001;
    
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

    % define the actual solution of f(x), using denser sample (1000) used.
    m = 1000;
    x      = -1: 2/m :1;
    f_act  = 1 ./(3 + x);
    
    % define the numerical solution of f(x), using denser sample (1000) used.
    V_num  = zeros(size(x,2),size(c,1));
    for i = 1:size(x,2)
        for j = 1:size(c,1)
            V_num(i,j) = x(i)^(j-1);
        end
    end
    
    f_num  = V_num * c;
    f_num_pert  = V_num * c_pert;
    
    % derive the calculation difference
    diff = abs(f_num_pert'-f_num');

    % plot the error.
    subplot(2,1,2);semilogy(x,diff);
    xlabel('x')
    ylabel('log of   |q(x)- p(x)|')
    title(sprintf('Interpolation Error after single point small perturbation: |q(x)-p(x)| :  (Chebyshev points n = %d)',n));
end


%========================================================================
%% --------------------Problems and Computer Exercises : 4.1.1 (b) ------
%========================================================================

for n = 8:8:16
    
    %% (1) option 1: define the equidistant points
    x = zeros(1,n);
    for i = 1:n
        x(i) = -1 + 2*(i-1)/(n-1);
        x_disc = x;
    end
    f = random('Normal',0,1,1,n);
    
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
    
    
    % define the actual solution of f(x), using denser sample (1000) used.
    m = 1000;
    x      = -1: 2/m :1;
    f_act  = 1 ./(3 + x);
    
    % define the numerical solution of f(x), using denser sample (1000) used.
    V_num  = zeros(size(x,2),size(c,1));
    for i = 1:size(x,2)
        for j = 1:size(c,1)
            V_num(i,j) = x(i)^(j-1);
        end
    end

    f_num  = V_num * c;
    
    % plot the f(x) series.
    figure;subplot(2,1,1);plot(x_disc,f,':bs',x,f_num,'r');
    axis([-1 1 -10 10])
    xlabel('x')
    ylabel('f(x)')
    title(sprintf('Interpolation of f(x) by p(x) :  (equidistant points n = %d)',n));

    

     
    %% (2) option 2: define the Chebyshev points
    x = zeros(1,n);
    for i = 1:n
        x(i) = cos((2*i-1)*pi/(2*n));
        x_disc = x;
    end
%     f = random('Normal',0,1,1,n);
    
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

    % define the actual solution of f(x), using denser sample (1000) used.
    m = 1000;
    x      = -1: 2/m :1;
    f_act  = 1 ./(3 + x);
    
    % define the numerical solution of f(x), using denser sample (1000) used.
    V_num  = zeros(size(x,2),size(c,1));
    for i = 1:size(x,2)
        for j = 1:size(c,1)
            V_num(i,j) = x(i)^(j-1);
        end
    end
    
    f_num  = V_num * c;
    
    % plot the f(x) series.
    subplot(2,1,2);plot(x_disc,f,':bs',x,f_num,'r');
    axis([-1 1 -10 10])
    xlabel('x')
    ylabel('f(x)')
    title(sprintf('Interpolation of f(x) by p(x) :  (Chebyshev points n = %d)',n));
end