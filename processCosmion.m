%units CGS everywhere

load positions.dat
% positions = positions_SHO

%% 
GeV = 1.78e-24;
m = 2* GeV;
Rsun = 69.57d9;
kB = 1.38e-16;
r = sqrt(sum(positions(:,1:3).^2,2));
v = sqrt(sum(positions(:,4:6).^2,2));
vout = sqrt(sum(positions(:,7:9).^2,2));

T  = .5*m*v.^2/kB*(2/3);

t = positions(:,10);
dE = .5*m*(vout.^2-v.^2);
nbins = 70; %start here, see what happens
bins = linspace(0,Rsun,nbins)
L =bins*0;
Tav = bins*0
stdL = bins*0;
for i = 2:nbins        
    binned = (r>bins(i-1))& (r<=bins(i));
    ninbin = sum(binned);
    L(i) = sum(dE(binned));
    stdL(i) = sqrt(sum(dE(binned).^2)-L(i).^2/ninbin);
    Tav(i) = sum(T(binned))/ninbin;
end


%% plots 
% histogram(r)
L = L/sum(t);
stdL = stdL/sum(t);
% plot(bins,L,'linestyle','none')
errorbar(bins,L,stdL,'linestyle','none','markersize',10,'marker','.')
% plot(bins,Tav)