function MakeDotH(RFOut,gOut,dt)
rf_max = max(abs(RFOut(:)));
rf_normalized = RFOut/rf_max;
g01_max = max(abs(gOut(:,1)));
g01_norm = gOut(:,1)/g01_max;
g02_max = max(abs(gOut(:,2)));
g02_norm = gOut(:,2)/g02_max;
rf_phase = angle(RFOut);
for ii = 1:length(rf_phase)
    if rf_phase(ii) < 0
        rf_phase(ii) = rf_phase(ii) + pi*2;
    end
end

% Declaring amp, phase and max amp. 
str_header = 'SingleSpiral.h';

str_grad_x = sprintf('%.5ff, ', g01_norm);
str_grad_x = strcat('const float grissom_grad01_norm [] = { ', str_grad_x, ' };');
str_grad_x_neg = sprintf('%.5ff, ', -g01_norm);
str_grad_x_neg = strcat('const float grissom_grad01_norm_neg [] = { ', str_grad_x_neg, ' };');
str_grad_y = sprintf('%.5ff, ', g02_norm);
str_grad_y = strcat('const float grissom_grad02_norm [] = { ', str_grad_y, ' };');
str_grad_y_neg = sprintf('%.5ff, ', -g02_norm);
str_grad_y_neg = strcat('const float grissom_grad02_norm_neg [] = { ', str_grad_y_neg, ' };');
grad_pts =  size(gOut,1);
str_grad_pts = sprintf('long grissom_grad_pts = %d;',grad_pts);
str_grad_duration = sprintf('long grissom_grad_duration = %d;',round(dt*1e6*grad_pts));


str_rf_mag = sprintf('%.5ff, ', abs(rf_normalized));
str_rf_mag = strcat('const float grissom_rf_norm [] = { ', str_rf_mag, ' };');
str_rf_phase = sprintf('%.5ff, ', rf_phase);
str_rf_phase = strcat('const float grissom_rf_phase [] = { ', str_rf_phase, ' };');
%str_rf_max = sprintf('const float grissom_rf_amp = %f;', rf_max);
str_rf_pts = sprintf('const long grissom_rf_pts = %d;',length(RFOut));
str_rf_AmpInt = sprintf('float grissom_rf_ampint = %sf;', num2str(sum(abs(rf_normalized))));
str_rf_duration = sprintf('const long grissom_rf_duration = %d;',round(dt*1e6*length(RFOut)));


str_grad01_amp = sprintf('const float grissom_grad01_amp = %sf;', num2str(g01_max));
str_grad02_amp = sprintf('const float grissom_grad02_amp = %sf;', num2str(g02_max));

% Write arrays into C++ header files.
FID = fopen(strcat('/Users/ytong/Documents/XP/',str_header), 'w+');
fprintf(FID, '%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s', str_grad_x, str_grad_y,...
    str_grad_x_neg, str_grad_y_neg, str_grad_pts, str_grad_duration, str_rf_mag, str_rf_phase,...
    str_rf_AmpInt, str_rf_pts, ...
    str_grad01_amp, str_grad02_amp, str_rf_duration);
fclose(FID);
