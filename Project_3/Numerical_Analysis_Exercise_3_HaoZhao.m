%% ========================================================================
%Exercise 3 of Numerical Analysis Course, Coded by Hao Zhao, May 9,2016
%Exericise of Gauss–Legendre quadrature 
%% ========================================================================

close all;clc;clear;

%========================================================================
%% --------------------Problems and Computer Exercises : 5.3.8 (a) ------
%========================================================================


    %% (1)  define the nodes and associated weights based on Gauss-Legendre quadrature,derived from Table 5.3.1
    x = [-0.906179845938664 -0.538469310105683 0 0.538469310105683 0.906179845938664];
    w = [ 0.236926885056189  0.478628670499366 0.568888888888889 0.478628670499366 0.236926885056189];
    
    I  = 0;
    fx = 0;
    
    for i=1:size(x,2)
        fx = (x(i)+1)^4*(sin((x(i)+1)*pi/2))^2; 
        I  = I + w(i)* fx; 
    end
    
    I = I/16;
    
    Result = sprintf('Integration reuslt by Gauss–Legendre quadrature rule P = %f ',I);
    disp(Result)