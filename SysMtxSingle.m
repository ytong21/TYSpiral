function DA = SysMtxSingle(k,b0,b1,mask,dt,tb0,Positions)
%for single channel only now!
Nt = size(k,1);
Ns = numel(b1);
PulseDuration = dt*(Nt-1);

Dr = diag(b1);
GAMMA = 42.577E6; % in Hz/T

m0 = 1;
AA = zeros(Ns, Nt);
Positions(3,:,:) = [];
Positions = Positions(:,:).';
Positions = Positions(mask,:);
Positions = conj(Positions);
for ii = 1:Ns
    for jj = 1:Nt
        AA(ii,jj) = (1i*GAMMA*dt*m0) * exp(1i*GAMMA*b0(ii)*(tb0(jj)-PulseDuration)) ...
            * exp(1i*k(jj,:)*Positions(ii,:)');
    end
end
DA = Dr*AA;