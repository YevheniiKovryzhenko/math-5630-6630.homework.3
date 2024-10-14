% Author: Yevhenii Kovryzhenko / yzk0058@auburn.edu
% Date: 2024-09-01
% Assignment Name: hw03

classdef hw03
    methods (Static)

        function y = p1(data, eval)
            % Lagrange interpolation method
            %
            % :param: data: a matrix of size n x 2, where n is the number of data points
            %         data(:,1) is the x values
            %         data(:,2) is the y values
            % :param: eval: a vector of x values to evaluate the interpolating polynomial
            % :return: a vector, the evaluations of resulting interpolating polynomial at x values in eval

            n_eval = length(eval);
            
            [n_poly, ~] = size(data);
            x_interp = data(1:n_poly, 1);
            y_interp = data(1:n_poly, 2);
            y = zeros(n_eval, 1);
            
            for i_eval = 1:n_eval
                x = eval(i_eval);

                q = 1;
                for i = 1:n_poly
                    q = q*(x - x_interp(i));
                end
                
                tmp = 0;
                for j = 1:n_poly
                    kj = 1;
                    for s = 1:n_poly
                        if s == j
                            continue
                        end
                        kj = kj/(x_interp(j) - x_interp(s));
                    end

                    tmp = tmp + kj*y_interp(j) / (x - x_interp(j));
                end
                y(i_eval) = tmp*q;
            end
        end

        function p2()
            % Use equally spaced nodes and Chebyshev nodes to compute the Lagrange polynomial interpolation of 
            % f(x) = 1/(1 + 25x^2) and g(x) = sin(pi * x) on [-1, 1]. 
            % This code uses your implementation of p1 to compute the
            % interpolation, and record the maximum interpolation error at 1000 equally spaced points in [-1, 1].
            % ----------------------------------------------------------------------------------

            % First, run this function and tabulate your result in the table below. 
            % Then, make a comment/explanation on the trend of the error for **equally spaced nodes** as n increases for each function.

            % Write your comments here.
            % It is clear that as N grows, the interpolation error for
            % equally spaced nodes grows and is much greater than
            % interpolation that uses Chebyshev nodes. This is due to the
            % Runge's phenomenon.
            %
            %
            %
            %
            %
            

            % |n  |                        Function       | Error with equally spaced nodes | Error with Chebyshev nodes  |
            % |---|---------------------------------------|---------------------------------|-----------------------------|
            % |  5|                @(x)1./(1+25*x.^2)     | 4.3267e-01                      |  5.5589e-01                 |
            % | 10|                @(x)1./(1+25*x.^2)     | 1.9156e+00                      |  1.0915e-01                 |
            % | 15|                @(x)1./(1+25*x.^2)     | 2.1069e+00                      |  8.3094e-02                 |
            % | 20|                @(x)1./(1+25*x.^2)     | 5.9768e+01                      |  1.5333e-02                 |
            % | 25|                @(x)1./(1+25*x.^2)     | 7.5764e+01                      |  1.1411e-02                 |
            % | 30|                @(x)1./(1+25*x.^2)     | 2.3847e+03                      |  2.0613e-03                 |
            % | 35|                @(x)1./(1+25*x.^2)     | 3.1708e+03                      |  1.5642e-03                 |
            % | 40|                @(x)1./(1+25*x.^2)     | 1.0438e+05                      |  2.8935e-04                 |
            % | 45|                @(x)1./(1+25*x.^2)     | 1.4243e+05                      |  2.1440e-04                 |
            % | 50|                @(x)1./(1+25*x.^2)     | 4.8178e+06                      |  3.9629e-05                 |
            % | 55|                @(x)1./(1+25*x.^2)     | 6.6475e+06                      |  2.9383e-05                 |
            % |---|---------------------------------------|---------------------------------|-----------------------------|
            % |  5|                     @(x)sin(pi*x)     | 2.6754e-02                      |  1.3193e-02                 |
            % | 10|                     @(x)sin(pi*x)     | 5.1645e-05                      |  6.0348e-06                 |
            % | 15|                     @(x)sin(pi*x)     | 9.2162e-10                      |  2.0995e-11                 |
            % | 20|                     @(x)sin(pi*x)     | 3.0112e-13                      |  1.3323e-15                 |
            % | 25|                     @(x)sin(pi*x)     | 1.2485e-11                      |  1.1102e-15                 |
            % | 30|                     @(x)sin(pi*x)     | 1.6042e-10                      |  1.3323e-15                 |
            % | 35|                     @(x)sin(pi*x)     | 5.8851e-09                      |  1.6653e-15                 |
            % | 40|                     @(x)sin(pi*x)     | 1.2282e-07                      |  1.7764e-15                 |
            % | 45|                     @(x)sin(pi*x)     | 7.5566e-06                      |  2.1094e-15                 |
            % | 50|                     @(x)sin(pi*x)     | 1.7865e-04                      |  2.1094e-15                 |
            % | 55|                     @(x)sin(pi*x)     | 5.0905e-03                      |  2.2204e-15                 |
            % |---|---------------------------------------|---------------------------------|-----------------------------|

            eval_pts = linspace(-1, 1, 1000)';
            funcs = { @(x) 1 ./ (1 + 25 * x.^2), @(x) sin(pi * x) };
            fprintf('|n  |                        Function       | Error with equally spaced nodes | Error with Chebyshev nodes  |\n');
            fprintf('|---|---------------------------------------|---------------------------------|-----------------------------|\n')
            for i = 1:2
                func = funcs{i};
                for n = 5:5:55
                    % Equally spaced nodes
                    x = linspace(-1, 1, n+1);
                    y = func(x);
                    y_eval = hw03.p1([x', y'], eval_pts);
                    eq_error_f = max(abs(func(eval_pts) - y_eval));
                    % Chebyshev nodes
                    x = cos((2 * (1:(n+1)) - 1) * pi / (2 * n + 2));
                    y = func(x);
                    y_eval = hw03.p1([x', y'], eval_pts);
                    cheby_error_f = max(abs(func(eval_pts) - y_eval));
                    fprintf('| %2d|    %30s     | %6.4e                      |  %6.4e                 |\n', n, functions(func).function, eq_error_f, cheby_error_f);
                end
                    fprintf('|---|---------------------------------------|---------------------------------|-----------------------------|\n')
            end
        end

        function p3()
            % **Only for 6630**.
            % Use the extreme Chebyshev nodes to compute the Lagrange polynomial interpolation of 
            % f(x) = 1/(1 + 25x^2) and g(x) = sin(pi * x) on [-1, 1]. 
            % This code uses your implementation of p1 to compute the
            % interpolation, and record the maximum interpolation error at 1000 equally spaced points in [-1, 1].
            % ----------------------------------------------------------------------------------

            % Run this function and tabulate your result in the table below. 
            % Then, make a comment on the performance of the extreme Chebyshev nodes compared to Chebyshev nodes.

            % Write your comment here.
            %
            % Given the two example functions, it does look like extreme
            % Chebyshev nodes are have interpolation errors on the same 
            % order as the regular Chebyshev nodes, but slightly worse. 
            %
            %
            %

            % |n  |                        Function       | Error with extreme Chebyshev nodes  |
            % |---|---------------------------------------|-------------------------------------|
            % |  5|                @(x)1./(1+25*x.^2)     | 6.3862e-01                          |
            % | 10|                @(x)1./(1+25*x.^2)     | 1.3219e-01                          |
            % | 15|                @(x)1./(1+25*x.^2)     | 9.9308e-02                          |
            % | 20|                @(x)1./(1+25*x.^2)     | 1.7738e-02                          |
            % | 25|                @(x)1./(1+25*x.^2)     | 1.3649e-02                          |
            % | 30|                @(x)1./(1+25*x.^2)     | 2.4252e-03                          |
            % | 35|                @(x)1./(1+25*x.^2)     | 1.8710e-03                          |
            % | 40|                @(x)1./(1+25*x.^2)     | 3.3987e-04                          |
            % | 45|                @(x)1./(1+25*x.^2)     | 2.5645e-04                          |
            % | 50|                @(x)1./(1+25*x.^2)     | 4.6187e-05                          |
            % | 55|                @(x)1./(1+25*x.^2)     | 3.5147e-05                          |
            % |---|---------------------------------------|-------------------------------------|
            % |  5|                     @(x)sin(pi*x)     | 1.3357e-02                          |
            % | 10|                     @(x)sin(pi*x)     | 1.1730e-05                          |
            % | 15|                     @(x)sin(pi*x)     | 2.1000e-11                          |
            % | 20|                     @(x)sin(pi*x)     | 1.6653e-15                          |
            % | 25|                     @(x)sin(pi*x)     | 1.5543e-15                          |
            % | 30|                     @(x)sin(pi*x)     | 1.3323e-15                          |
            % | 35|                     @(x)sin(pi*x)     | 1.3323e-15                          |
            % | 40|                     @(x)sin(pi*x)     | 2.2204e-15                          |
            % | 45|                     @(x)sin(pi*x)     | 1.8874e-15                          |
            % | 50|                     @(x)sin(pi*x)     | 1.9984e-15                          |
            % | 55|                     @(x)sin(pi*x)     | 1.8874e-15                          |
            % |---|---------------------------------------|-------------------------------------|

            eval_pts = linspace(-1, 1, 1000)';
            funcs = { @(x) 1 ./ (1 + 25 * x.^2), @(x) sin(pi * x) };
            fprintf('|n  |                        Function       | Error with extreme Chebyshev nodes  |\n');
            fprintf('|---|---------------------------------------|-------------------------------------|\n')
            for i = 1:2
                func = funcs{i};
                for n = 5:5:55
                    % extreme Chebyshev nodes
                    x = cos(((1:n+1) - 1) * pi / (n));
                    y = func(x);
                    y_eval = hw03.p1([x', y'], eval_pts);
                    ex_cheby_error_f = max(abs(func(eval_pts) - y_eval));
                    fprintf('| %2d|    %30s     | %6.4e                          |\n', n, functions(func).function, ex_cheby_error_f);
                end
                fprintf('|---|---------------------------------------|-------------------------------------|\n')
            end
        end
    end
end