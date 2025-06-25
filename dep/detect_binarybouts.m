function [bouts, lens] = detect_binarybouts(binary_vector)
% Safer version: handles edge cases and avoids size mismatch errors

diff_v = diff(binary_vector);
v_s = find(diff_v > 0) + 1;  % bout starts (rising edges)
v_e = find(diff_v < 0);      % bout ends (falling edges)

% Handle if bout starts at frame 1
if binary_vector(1) == 1
    v_s = [1; v_s(:)];  % start fix
end

% Handle if bout ends at final frame
if binary_vector(end) == 1
   v_e = [v_e(:); length(binary_vector)];
end

% Now, make sure the number of starts matches number of ends
if length(v_s) ~= length(v_e)
    % Try to fix it, otherwise return empty
    min_len = min(length(v_s), length(v_e));
    v_s = v_s(1:min_len);
    v_e = v_e(1:min_len);
end

bouts = [v_s(:), v_e(:)];
lens = v_e(:) - v_s(:);
end
