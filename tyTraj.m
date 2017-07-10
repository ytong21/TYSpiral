function [kTraj,MaxRadius] = tyTraj(gHzpercm,dt)
% g in Hz/cm
kTraj = zeros(size(gHzpercm));
kTraj(1,:) = gHzpercm(1,:);
for ii = 2:size(gHzpercm,1)
    kTraj(ii,1) = kTraj(ii-1,1) + gHzpercm(ii,1);
    kTraj(ii,2) = kTraj(ii-1,2) + gHzpercm(ii,2);
end
kTraj = kTraj*dt; % to convert to 1/cm;  
kRadius = hypot(kTraj(:,1),kTraj(:,2));
MaxRadius = max(kRadius(:));