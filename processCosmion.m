%units CGS everywhere

load positions.dat

%% 
m = 1* 1.78e-23;
Rsun = 69.57d9;
r = sqrt(sum(positions(:,1:3).^2,2));
v = sqrt(sum(positions(:,4:6).^2,2));
vout = sqrt(sum(positions(:,7:9).^2,2));
t = positions(:,10);
dE = .5*m*(vout.^2-v.^2);
nbins = 50; %start here, see what happens
bins = linspace(0,Rsun,nbins)
for i = 2:nbins    
    L(i) = sum(dE((r>bins(i-1))& (r<=bins(i))));
end
L = L/sum(t);
plot(bins,L)