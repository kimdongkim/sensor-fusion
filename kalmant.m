function [xh yh] = kalmant(ax, ay, xm, ym)
%
%
persistent A H Q R B
persistent x u P
persistent firstRun

if isempty(firstRun)
  dt = 0.1;
  
  A = [ 1  dt  0   0
        0  1   0   0
        0  0   1  dt
        0  0   0   1 ];

  B = [0.5*dt^2, 0;
     dt, 0;
     0, 0.5*dt^2;
     0, dt];
  
  H = [ 1  0  0  0
        0  0  1  0 ];
 
  Q = 0.05*eye(4);
  R = [ 1  0
        0  1 ];

  x = [0, 0, 0, 0]';
  P = 1*eye(4);
  
  firstRun = 1;
end

u = [ax ay]';
xp = A*x+B*u;
Pp = A*P*A' + Q;

K = Pp*H'*inv(H*Pp*H' + R);

z = [xm ym]';
x = xp + K*(z - H*xp);
P = Pp - K*H*Pp;


xh = x(1);
yh = x(3);