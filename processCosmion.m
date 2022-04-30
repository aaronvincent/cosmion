%units CGS everywhere

load positions.dat

%% 
m = 10* 1.78e-24;
Rsun = 69.57d9;
kB = 1.38e-16;
r = sqrt(sum(positions(:,1:3).^2,2));
v = sqrt(sum(positions(:,4:6).^2,2));
vout = sqrt(sum(positions(:,7:9).^2,2));

T  = .5*m*v.^2/kB*(2/3);

t = positions(:,10);
dE = .5*m*(vout.^2-v.^2);
nbins = 200; %start here, see what happens
bins = linspace(0,Rsun,nbins)
for i = 2:nbins    
    L(i) = sum(dE((r>bins(i-1))& (r<=bins(i))));
end
L = L/sum(t);
plot(bins,L)