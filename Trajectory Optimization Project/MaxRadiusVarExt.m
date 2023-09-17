% Uses variation of extremals w/ free final state (final states in the obj fcn)

clear all

global mu T control

mu=3.986e5;     %gravity constant
ro=6800;        %initial radius
T=.0002;          %thrust acceleration (m/sec^2)
tfinal=4*3600;  %total time in sec
tstep=20;      %simulation step size
H11=4e4; H22=4e4;   %weighting factors for perfomance index

% time=[0:tstep:tfinal]';
% numsteps=size(time,1);

%Set initial guess for the costates and delta for finite difference
% Guess from Steepest descent
po=1000*[-0.0047   -0.9518   -4.9564]; nump=size(po,2); delp=.00001;

%Bootstrap method (Start with T=.01T)
%po=zeros(1,3); nump=size(po,2); delp=.00001;

xo=[ro 0 sqrt(mu/ro)];  numx=size(xo,2);
options=odeset('RelTol',1e-08,'AbsTol',1e-08*ones(6,1));

error=1e9;
itercount=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start of iteration
while error>.1

itercount=itercount+1;
% if itercount==4 || itercount==6 || itercount==8 
%     H11=H11*2;  H22=H22*2;
% end
    
[t,z]=ode45('MaxRadiusVarExtEoms', [0 tfinal], [xo po],options);

%final desired value of costates
rf=z(end,1); uf=z(end,2);  vf=z(end,3);
pf=[-1+.5*H22*(vf*sqrt(mu/rf^3)-mu/rf^2)   H11*uf    H22*(vf-sqrt(mu/rf))];
error=norm(z(end,4:6)-pf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%estimate new initial costates

%calculate influence matricies
for i=1:nump
    popert=po;
    popert(1,i)=po(1,i)+delp;
    [tpert,zpert]=ode45('MaxRadiusVarExtEoms', [0 tfinal], [xo popert],options);
    Pp(1:nump,i)=( zpert(end,4:6)'-z(end,4:6)' ) / delp;
    Px(1:numx,i)=( zpert(end,1:numx)'-z(end,1:numx)' ) / delp;
end
d2hdx2=[H22*(-.75*vf*sqrt(mu/rf^5)+mu/rf^3)  0  .5*H22*sqrt(mu/rf^3);
    0  H11  0;
    .5*H22*sqrt(mu/rf^3)  0 H22];
po=po+ (  inv(d2hdx2*Px - Pp) * (z(end,4:6)-pf)' )';


end  %End of iteration

'# of iterations'
itercount

%Calculate control
for i=1:size(z,1)
    phi(i,1)=atan2(-z(i,5),-z(i,6));
end

figure(1)
subplot(2,1,1)
plot(t,phi*180/pi)
xlabel('t (s)')
ylabel('phi (deg)')
grid on

subplot(2,1,2)
plot(t,z(:,1))
xlabel('t (s)')
ylabel('r (m)')
grid on


'final results'
disp(['r(tf)=' num2str(z(end,1)) ])
disp(['u(tf)=' num2str(z(end,2)) ])
disp(['v(tf)' num2str(z(end,3)) ])
disp(['v circular=' num2str(sqrt(mu/z(end,1))) ])


