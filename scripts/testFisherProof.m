theta0 = [1 .9 .2 ; 
          .9 1 -.5;
          .2 -.5 1];
theta1 = nearcorr(theta0);
x = eye(3) + (1-eye(3))*.1;
eta1 = atanh(theta1)
y = atanh(x) 
%yhat1 = eta1 + (1/(1-tanh(eta1)^2))*(x - tanh(eta1))
%eps1 = y - yhat1

t = .9
eta2 = t*eta1 + (1-t)*y
theta2 = tanh(eta2)

disp("===========")
disp(abs(x-theta2))
disp("<=")
disp(abs(x-theta0))