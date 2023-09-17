

function zdot=MaxRadiusVarExtEoms(t,z)

global mu T

phi=atan2(-z(5),-z(6));

%states
r=z(1); u=z(2); v=z(3);
zdot(1,1)=u;
zdot(2,1)=v^2/r-mu/r^2+T*sin(phi);
zdot(3,1)=-u*v/r+T*cos(phi);

%costates  
pr=z(4);  pu=z(5);  pv=z(6);
zdot(4,1)=-pu*(-v^2/r^2+2*mu/r^3)-pv*(v*u/r^2);
zdot(5,1)=-pr+pv*v/r;
zdot(6,1)=-pu*2*v/r+pv*u/r;

return
