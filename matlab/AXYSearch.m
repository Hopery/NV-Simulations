function AXY = AXYSearch(x0)
    options = optimset('Diagnostics','off', 'Display','off');
    [x] = fsolve(@myfunKDD5even,x0, options);
    x(end+1) = x(2) + x(4) - (x(1) + x(3)) + pi/2;
    AXY = x./(2*pi);
end

function F = myfunKDD5even(x)
global f1e f2e f3e f4e


theta1 = x(1);
theta2 = x(2);
theta3 = x(3);
theta4 = x(4);

theta5 = theta2 + theta4 -(theta1 + theta3) + pi/2;

F = [f1e    - 4/(1*pi)*(  sin(1*theta1) + sin(1*theta3) + sin(1*theta5) - sin(1*theta2) - sin(1*theta4) );  
     f2e    - 4/(2*pi)*(  sin(2*theta1) + sin(2*theta3) + sin(2*theta5) - sin(2*theta2) - sin(2*theta4) );      
     f3e    - 4/(3*pi)*(  sin(3*theta1) + sin(3*theta3) + sin(3*theta5) - sin(3*theta2) - sin(3*theta4) );
     f4e    - 4/(4*pi)*(  sin(4*theta1) + sin(4*theta3) + sin(4*theta5) - sin(4*theta2) - sin(4*theta4) );                  
        ];


end
