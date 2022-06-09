%units CGS everywhere

load positions_SHO.dat
% positions = positions_SHO2;

%% 
GeV = 1.78e-24;
mx = 4; %* GeV;
m = mx*GeV;
Rsun = 69.57d9;
sigma = 1e-37;
% kB = 1.38e-16;
% r = sqrt(sum(positions(:,1:3).^2,2));
% v = sqrt(sum(positions(:,4:6).^2,2));
% vout = sqrt(sum(positions(:,7:9).^2,2));
% 
% T  = .5*m*v.^2/kB*(2/3);
% 
% t = positions(:,10);
% dE = .5*m*(vout.^2-v.^2);
% nbins = 70; %start here, see what happens
% bins = linspace(0,Rsun,nbins)
% L =bins*0;
% Tav = bins*0
% stdL = bins*0;
% for i = 2:nbins        
%     binned = (r>bins(i-1))& (r<=bins(i));
%     ninbin = sum(binned);
%     L(i) = sum(dE(binned));
%     stdL(i) = sqrt(sum(dE(binned).^2)-L(i).^2/ninbin);
%     Tav(i) = sum(T(binned))/ninbin;
% end


%% plots 
% histogram(r)
% L = L/sum(t);
% stdL = stdL/sum(t);
% plot(bins,L,'linestyle','none')
% errorbar(bins,L,stdL,'linestyle','none','markersize',10,'marker','.')
% plot(bins,Tav)



xvec = positions_SHO(:,1:3);
vvec = positions_SHO(:,4:6);
potan = positions_SHO(:,end);

r = sqrt(sum(xvec.^2,2));

for j = 1:length(xvec)
xdotv(j) = xvec(j,:)*vvec(j,:)';
end


for j = 1:length(xvec)
xhat(j,:) = xvec(j,:)./r(j);
end
% size(xhat)
for j = 1:length(xvec)
    vr(j) = (vvec(j,:)*xhat(j,:)');
vperp(j,:) = vvec(j,:) - vr(j)*xhat(j,:);
end
vperpAmp = sqrt(sum(vperp.^2,2));

E = .5*sum(vvec.^2,2) + potan;

histogram(vr,[-1:.01:1]*1e8,'normalization','pdf')
hold on



load positions.dat

xvec = positions(:,1:3);
vvec = positions(:,4:6);
potnum = positions(:,end);

rnum = sqrt(sum(xvec.^2,2));

for j = 1:length(xvec)
xdotv(j) = xvec(j,:)*vvec(j,:)';
end


for j = 1:length(xvec)
xhat(j,:) = xvec(j,:)./rnum(j);
end
% size(xhat)
for j = 1:length(xvec)
    vrNum(j) = (vvec(j,:)*xhat(j,:)');
vperpNum(j,:) = vvec(j,:) - vrNum(j)*xhat(j,:);
end
vperpAmpNum = sqrt(sum(vperpNum.^2,2));
Enum = .5*sum(vvec.^2,2) + potnum;



histogram(vrNum,[-1:.01:1]*1e8,'normalization','pdf')
xlabel('$v_r$','fontsize',16,'interpreter','latex')
legend('SHO', 'NUMERICAl')
figure
histogram(vperpAmp,[-1:.01:1]*1e8,'normalization','pdf')    
hold on
histogram(vperpAmpNum,[-1:.01:1]*1e8,'normalization','pdf')

xlabel('$v_\theta$','fontsize',16,'interpreter','latex')
 %% 
figure
histogram(r,[0:.001:.2]*Rsun,'normalization','pdf')
hold on
rnum = sqrt(sum(positions(:,1:3).^2,2));
histogram(rnum,[0:.001:.2]*Rsun,'normalization','pdf')
xlabel('$r$','fontsize',16,'interpreter','latex')


[R, Etrans,Q, K, nx, sigsOut,nxIso,nxLTE, Ltrans,LPS,LLTE,Rchi] = luminosity_constrho_slim(sigma,mx ,0, 0,220e5,1,1);
plot(R*Rsun,R.^2*Rsun.^2.*nxIso./trapz(R*Rsun,R.^2*Rsun.^2.*nxIso),'linewidth',2)
set(gca,'xlim',[0,.2]*Rsun)


figure
histogram(E,'normalization','pdf')
hold on
histogram(Enum,'normalization','pdf')
xlabel('$E$','fontsize',16,'interpreter','latex')
