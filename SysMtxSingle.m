function DA = SysMtxSingle(k,b0,b1,dt,tb0,Positions)
%for single channel only now!
% k in (g ms)/cm.
k = k/1000; % k in (g s)/cm.
k = k*(42.57E6*1E-4); %k in Hz/cm
Nt = size(k,1);
Ns = numel(b1);
PulseDuration = dt*(Nt-1);
Dr = diag(b1);% in Hz
%GAMMA = 42.577E6; % in Hz/T
GAMMA = 2*pi;
%GAMMA = 1;

m0 = 1;
AA = zeros(Ns, Nt);
Positions(:,3) = [];

for ii = 1:Ns
    for jj = 1:Nt
        AA(ii,jj) = (-1i*GAMMA*dt*m0) * exp(-1i*GAMMA*b0(ii)*(tb0(jj)-PulseDuration)) ...
            * exp(1i*GAMMA*k(jj,:)*Positions(ii,:)');
    end
end
DA = Dr*AA;