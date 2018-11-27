%% ========================================================================
%Exercise 2 of Numerical Analysis Course, Coded by Hao Zhao, March 13,2016
%% ========================================================================

close all;clc;clear;

%========================================================================
%% --------------------Problems and Computer Exercises : 4.3.9 (b) ------
%========================================================================


    %% (1)  define the input; f(x,y) = a_0 + a_1 * x + a_2 * y + a_3 * xy
    x  = [0 , 1 ];
    y  = [0 , 1 ];
    
    % the know cornner point: f(0,0)= 1,f(1,0)= 2,f(0,1)= 3,f(1,1)= 5.
    n  = size(x,2);
    f  = zeros(size(x,2),size(x,2));
    for i = 1 : n
        for j = 1 : n
            if( i==1 && j==1)
                f(i,j)=1;
            elseif (i==1 && j==n)
                f(i,j)=2;
            elseif (i==n && j==1)
                f(i,j)=3;
            else
                f(i,j)=5;
            end
        end
    end
    
    % define the Vandermonde Matrix for the linear system V'c = F
    

    V = zeros(4,4);
    c = zeros(n,1);
    p = zeros(n,1);
    
    for i = 1 : n
        for j = 1 : n
            V((i-1)*2+j,:) = [ 1 , x(i) , y(j) , x(i)*y(j) ];
        end
    end

    % derive the F vector, assuming the q0 = 1.
    F  = [f(1),f(2),f(3),f(4)]';
    
    % solve the linear sytem. The bilineal interp follwing below
    % expression: f(x,y) = a_0 + a_1 * x + a_2 * y + a_3 * xy;
    c = V\F;


    
    % define the numerical solution of f(x), using derived C vector.
    
    x  = 0:0.01:1;
    y  = 1:-0.01:0;
    f_num_mat = zeros(size(x,2),size(y,2));
    
    for j = 1 : size(y,2)
        for i = 1 : size(x,2)
             f_num_mat(j,i)  = c(1) + c(2) * x(i) + c(3) * y(j) + c(4) * x(i) * y(j);
        end
    end
    
    imagesc(f_num_mat);
            
    % display the function and testing values.
    Funct = sprintf('Derived bilinearal Interpolation function f(x,y) = %d + %d * x + %d * y + %d * xy',c(1),c(2),c(3),c(4));
    disp(Funct);
    xlabel('x')
    ylabel('y')
    title('Plot of Bilinearal interploation: (4 conner points used for interplation f(x,y), x = 0:0.01:1,y = 0:0.01:1)');
    
    set(gca, 'XTick', 1:10:101);
    set(gca, 'XTickLabel', 0:0.1:1);
    set(gca, 'yTick', 1:10:101);
    set(gca, 'yTickLabel', 1:-0.1:0);

    
    % define the numerical solution of f(x), using derived C vector.
    x   = 0.5;
    y   = 0.25;
    f_num  = c(1) + c(2) * x + c(3) * y + c(4) * x*y;
    
    Test = sprintf(' For x = 0.5,y=0.25 The bilinear interpolated value is: %f',f_num);
    disp(Test)  
    
    
% %     % plot the f(x) series.
% %     subplot(3,1,1); plot(x,f_act,'r','LineWidth',2);grid on;
% %     axis([-0.15 0.15 -500 500])
% %     xlabel('x (radian)')
% %     ylabel('f(x)')
% %     title('Plot of function: f = cot(x) :  (x = -10 to 10 degree)');
% %     
% %     % plot the interploated series.
% % 
% %     subplot(3,1,2);plot((1:1:5)*pi/180,cot((1:1:5)*pi/180),':ks','LineWidth',2);hold on
% %     subplot(3,1,2);plot(x,f_num,'--b','LineWidth',2);grid on;
% %     axis([-0.15 0.15 -500 500])
% %     xlabel('x (radian)')
% %     ylabel('f(x)')
% %     title('Plot of function by Rational Interploation:  ( discrete point x = 1 to 5 degree used)');
% % 
% %     
% %     % plot the error.
% %     subplot(3,1,3);semilogy(x,abs(f_num-f_act),'--k','LineWidth',2);grid on;
% %     axis([-0.15 0.15 0 1e5])
% %     xlabel('x (radian)')
% %     ylabel('f(x)')
% %     title('Plot of Error of Rational interpolation function and cot(x) function');
% %     
% %     
% %     
% %     % display the function and testing values.
% % 
% %     Funct = sprintf('Derived Rational Interpolation function r = (%e + %e *x + %e *x^2)/(%d + %e *x + %e *x^2)',c(1),c(2),c(3),q0,c(4),c(5));
% %     disp(Funct)
% %     
% %     x_test      = 2.5 *pi/180;
% %     f_test_act  = cot(x_test);
% %     f_test_num  = (c(1) + c(2)*x_test + c(3)*x_test^2)./(q0  + c(4)*x_test + c(5)*x_test^2);
% %     diff_test   = abs(f_test_act - f_test_num);
% %     Test = sprintf(' For x = 2.5 deg, cot(x) is %e, F(x) from rational interpolation is: %e , error is: %e',f_test_act,f_test_num,diff_test);
% %     disp(Test)  
% %   