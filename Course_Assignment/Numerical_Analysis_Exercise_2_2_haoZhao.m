%% ========================================================================
%Exercise 2 of Numerical Analysis Course, Coded by Hao Zhao, March 13,2016
%% ========================================================================

close all;clc;clear;

%========================================================================
%% --------------------Problems and Computer Exercises : 4.3.12 (a) ------
%========================================================================


    %% (1)  define the input
    x  = (1:1:5)*pi/180;
    f  = cot(x);
  
    % define the Vandermonde Matrix for the linear system V'c = F
    n = size(x,2);
    V1 = zeros(n,n-2);
    V2 = zeros(n,2);
    c = zeros(n,1);
    p = zeros(n,1); 
    for i = 1:n
        for j = 1: n-2
            V1(i,j) = x(i)^(j-1);
        end
    end
    for i = 1:n
        for j = 1: 2
            V2(i,j) = -f(i)*x(i)^(j);
        end
    end
    
    % combine the two matrices into coefficient matrix    
    V = cat(2,V1,V2);

    % derive the F vector, assuming the q0 = 1.
    q0 = 1;
    F  = q0 .* f;
    
    % solve the linear sytem. The rational interp follwing below
    % expression: r = (c(1) + c(2)*x + c(3)*x^2)/(q0 + c(4)*x + c(5)*x^2);
    c = V\F';

    % define the actual solution of f(x)
    x      = (-10:0.1:10)*pi/180;
    f_act  = cot(x);
    
    % define the numerical solution of f(x), using derived C vector.
    f_num  = (c(1) + c(2).*x + c(3).*x.^2)./(q0  + c(4).*x + c(5).*x.^2);
    
    % plot the f(x) series.
    subplot(3,1,1); plot(x,f_act,'r','LineWidth',2);grid on;
    axis([-0.15 0.15 -500 500])
    xlabel('x (radian)')
    ylabel('f(x)')
    title('Plot of function: f = cot(x) :  (x = -10 to 10 degree)');
    
    % plot the interploated series.

    subplot(3,1,2);plot((1:1:5)*pi/180,cot((1:1:5)*pi/180),':ks','LineWidth',2);hold on
    subplot(3,1,2);plot(x,f_num,'--b','LineWidth',2);grid on;
    axis([-0.15 0.15 -500 500])
    xlabel('x (radian)')
    ylabel('f(x)')
    title('Plot of function by Rational Interploation:  ( discrete point x = 1 to 5 degree used)');

    
    % plot the error.
    subplot(3,1,3);semilogy(x,abs(f_num-f_act),'--k','LineWidth',2);grid on;
    axis([-0.15 0.15 0 1e5])
    xlabel('x (radian)')
    ylabel('f(x)')
    title('Plot of Error of Rational interpolation function and cot(x) function');
    
    
    
    % display the function and testing values.

    Funct = sprintf('Derived Rational Interpolation function r = (%e + %e *x + %e *x^2)/(%d + %e *x + %e *x^2)',c(1),c(2),c(3),q0,c(4),c(5));
    disp(Funct)
    
    x_test      = 2.5 *pi/180;
    f_test_act  = cot(x_test);
    f_test_num  = (c(1) + c(2)*x_test + c(3)*x_test^2)./(q0  + c(4)*x_test + c(5)*x_test^2);
    diff_test   = abs(f_test_act - f_test_num);
    Test = sprintf(' For x = 2.5 deg, cot(x) is %e, F(x) from rational interpolation is: %e , error is: %e',f_test_act,f_test_num,diff_test);
    disp(Test)  
    

%========================================================================
%% --------------------Problems and Computer Exercises : 4.3.12 (b) ------
%========================================================================
    
     %% (2) option 1: define the input
    x  = (1:1:5)*pi/180;
    f  = cot(x);
  
    % define the Vandermonde Matrix for the linear system V'c = F
    n = size(x,2);
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
    x      = (-10:0.1:10)*pi/180;
    f_act  = cot(x);
    
    % define the numerical solution of f(x), using denser sample (1000) used.
    V_num  = zeros(size(x,2),size(c,1));
    for i = 1:size(x,2)
        for j = 1:size(c,1)
            V_num(i,j) = x(i)^(j-1);
        end
    end
    f_num  = (V_num * c)';
    
    % plot the f(x) series.
    figure;
    subplot(3,1,1); plot(x,f_act,'r','LineWidth',2);grid on;
    axis([-0.15 0.15 -500 500])
    xlabel('x (radian)')
    ylabel('f(x)')
    title('Plot of function: f = cot(x) :  (x = -10 to 10 degree)');
    
    % plot the interpolated series.

    subplot(3,1,2);plot((1:1:5)*pi/180,cot((1:1:5)*pi/180),':ks','LineWidth',2);hold on
    subplot(3,1,2);plot(x,f_num,'--b','LineWidth',2);grid on;
    axis([-0.15 0.15 -500 500])
    xlabel('x (radian)')
    ylabel('f(x)')
    title('Plot of function by polynomials interpolation:  ( discrete point x = 1 to 5 degree used)');

    
    % plot the error.
    subplot(3,1,3);semilogy(x,abs(f_num-f_act),'--k','LineWidth',2);grid on;
    axis([-0.15 0.15 0 1e5])
    xlabel('x (radian)')
    ylabel('f(x)')
    title('Plot of Error of polynomials interpolation function and cot(x) function');
    
    Funct = sprintf('Derived Polynomial Interpolation function P(x) = (%e + %e *x + %e *x^2 + %e *x^3 + %e *x^4)',c(1),c(2),c(3),c(4),c(5));
    disp(Funct)
    
    % define the numerical solution of f(x), using derived C vector.
    x_test      = 2.5 *pi/180;
    f_test_act  = cot(x_test);
   
    f_test_num  = c(1) + c(2)*x_test + c(3)*x_test^2  + c(4)*x_test^3 + c(5)*x_test^4;
    diff_test   = abs(f_test_act - f_test_num);
    Test = sprintf(' For x = 2.5 deg, cot(x) is %e, F(x) from polynomials interpolation is: %e , error is: %e',f_test_act,f_test_num,diff_test);
    disp(Test)  
    
    
    
    
    
    
    
    
    
    
    
    
    
  