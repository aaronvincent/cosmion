Omega_SHO = 0.006452
xi = [2705710525.4906921,       -3873534938.3634562 ,      -2681433813.0402393 ]
vi = [-11372871.430080282,        73591.840957018765,       -16765518.336228890     ]

phase_i = atan2(-vi,xi*Omega_SHO);
A_i = xi./cos(phase_i);

t = linspace(0,3e3,1000)

for k = 1:length(t)
x = A_i.*cos(Omega_SHO*t(k) + phase_i);
r(k) = sqrt(sum(x.^2,2));
v = -A_i.*Omega_SHO.*sin(Omega_SHO*t(k) + phase_i);
xhat = x./sqrt(x*x');
vr(k) = v*xhat';
vmag(k) = sqrt(v*v');

end

plot(t,r)