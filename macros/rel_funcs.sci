//
// Limit state functions from
// J. F. Unger, D. Roos, "Investigation and benchmark of algorithms for reliability analysis"
//

//
// Limit state function 1
// Inputs are normally distributed
//

function y = get_mean_limit_state_1(n)
y = zeros(n,1);
endfunction

function y = get_var_limit_state_1(n)
y = ones(n,1);
endfunction

function y = limit_state_1(x,beta)
y = beta*sqrt(length(x)) + sum(x);
endfunction

function y = failure_limit_state_1(x,beta)
y = exp(-beta);
endfunction

//
// Limit state function 2 (convex)
// Inputs are exponentially distributed
//

function y = get_lambda_limit_state_2(n)
y = ones(n,1);
endfunction

function y = limit_state_2(x,C)
y = C - sum(x);
endfunction

function y = failure_limit_state_2(x,C)
y = 0.01;
endfunction

//
// Limit state function 3 (non-convex)
// Inputs are exponentially distributed
//

function y = get_lambda_limit_state_3(n)
y = ones(n,1);
endfunction

function y = limit_state_3(x,C)
y = - C + sum(x);
endfunction

function y = failure_limit_state_3(x,C)
y = 0.01;
endfunction

//
// Limit state function 4
// Inputs are normally distributed
// This function has several failure points
//

function y = get_mean_limit_state_4()
y = [0;0];
endfunction

function y = get_var_limit_state_4()
y = [1;1];
endfunction

function y = limit_state_4(x,C)
y = 5 - abs(x(1) + x(2));
endfunction

function y = failure_limit_state_4(x,C)
y = 0.01;
endfunction

//
// Limit state function 5
// Inputs are normally distributed
// This function has several failure points
//

function y = get_mean_limit_state_5()
y = [0;0];
endfunction

function y = get_var_limit_state_5()
y = [2;2];
endfunction

function y = limit_state_5(x)
y = 5 - abs(x(1) + x(2));
endfunction

function y = failure_limit_state_5(x)
y = 0.01;
endfunction

//
// Limit state function 6
// Inputs 1 to 5 are log-normal distributed (mean = 60, var = 6)
// Input 6 is Gumbel distributed (mean = 20, var = 6)
// Input 7 is Gumbel distributed (mean = 25, var = 7.5)
// This function has several failure points
//

function y = get_mean_limit_state_6()
y = [60;60;60;60;60;20;25];
endfunction

function y = get_var_limit_state_6()
y = [6;6;6;6;6;6;7.5];
endfunction

function y = limit_state_5(x)
aux_a = x(1) + 2*x(3) + 2*x(4) + x(5) - 5*x(6) - 5*x(7);
aux_b = x(1) + 2*x(2) + x(4) + x(5) - 5*x(6);
aux_c = x(2) + 2*x(3) + x(4) - 5*x(7);
y = min([aux_a aux_b aux_c]);
endfunction

function y = failure_limit_state_5(x)
y = 0.631;
endfunction

//
// Limit state function 10
// Inputs are normally distributed (mean = 0, var = 1)
//

function y = get_mean_limit_state_10()
y = [0;0];
endfunction

function y = get_var_limit_state_10()
y = [1;1];
endfunction

function y = limit_state_10(x)
y = 2 - (x(1) + x(2)) + 0.05*(sin(100*x(1)) + cos(100*x(2)));
endfunction

function y = failure_limit_state_10(x)
y = 0.01;
endfunction

//
// Limit state function 11
// Inputs are normally distributed (mean = 0, var = 1)
//

function y = get_mean_limit_state_11()
y = [0;0];
endfunction

function y = get_var_limit_state_11()
y = [1;1];
endfunction

function y = limit_state_11(x)
y = 1 - abs(x(1) + x(2) + 0.5) + 0.9*(x(1) + x(2))
endfunction

function y = failure_limit_state_11(x)
y = 0.01;
endfunction

