function DA = SysMtx(k,mask,Nc,b1,dt,tb0,b0,Positions)
%for single channel only now!
Nt = size(k,1);
Ns = numel(b1);
PulseDuration = dt*(Nt-1);

b0 = b0.*repmat(mask,[1 1 Nc]);
b0 = permute(b0,[3 1 2]);b0 = b0(:,:).';
b0 = b0(mask,:);
b0 = conj(b0);

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