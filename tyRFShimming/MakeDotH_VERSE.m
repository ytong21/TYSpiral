function MakeDotH_VERSE(RFOut,gOut,dt,phase_vector_1mm)
rf_max = max(abs(RFOut(:)));
rf_normalized = RFOut/rf_max;
if numel(gOut) == size(gOut,1) ||numel(gOut) == size(gOut,2)
    g01_max = max(abs(gOut));
else
    disp('1D gradient only')
    return
end
if isstruct(dt)
    dt_rf = dt.rf;
    dt_grad = dt.grad;
elseif isnumeric(dt)
    dt_rf = dt;
    dt_grad = dt;
end
g01_norm = gOut/g01_max;
% rf_phase = angle(RFOut);
% for ii = 1:length(rf_phase)
%     if rf_phase(ii) < 0
%         rf_phase(ii) = rf_phase(ii) + pi*2;
%     end
% end
rf_phase = phase_vector_1mm;
% Declaring amp, phase and max amp. 
str_header = 'VERSE.h';

str_grad_x = sprintf('%.5ff, ', g01_norm);
str_grad_x = strcat('const float verse_grad01_norm [] = { ', str_grad_x, ' };');
grad_pts =  numel(gOut);
str_grad_pts = sprintf('long verse_grad_pts = %d;',grad_pts);
str_grad_duration = sprintf('long verse_grad_duration = %d;',round(dt_grad*1e6*grad_pts));


str_rf_mag = sprintf('%.5ff, ', abs(rf_normalized));
str_rf_mag = strcat('const float verse_rf_norm [] = { ', str_rf_mag, ' };');




str_rf_phase = sprintf('%.5ff, ', rf_phase);
str_rf_phase = strcat('const float verse_rf_phase_one_mm [] = { ', str_rf_phase, ' };');
%str_rf_max = sprintf('const float verse_rf_amp = %f;', rf_max);
str_rf_pts = sprintf('const long verse_rf_pts = %d;',length(RFOut));
str_rf_AmpInt = sprintf('float verse_rf_ampint = %sf;', num2str(sum(abs(rf_normalized))));
str_rf_duration = sprintf('const long verse_rf_duration = %d;',round(dt_rf*1e6*length(RFOut)));


str_grad01_amp = sprintf('const float verse_grad01_amp = %.5ff;', g01_max);

% Write arrays into C++ header files.
FID = fopen(strcat('/Users/ytong/Documents/XP/',str_header), 'w+');
fprintf(FID, '%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s', str_grad_x,...
    str_grad_pts, str_grad_duration, str_rf_mag, str_rf_phase,...
    str_rf_AmpInt, str_rf_pts, ...
    str_grad01_amp, str_rf_duration);
fclose(FID);
