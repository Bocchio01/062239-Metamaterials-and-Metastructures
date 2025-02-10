function [force, final_time] = tone_burst(t_vector, fe, dfe)

number_cycles = floor(2 * fe / dfe);
final_time = number_cycles / fe;

signal = @(t) 1 / 4.6360 * sin(2*pi*fe * t) .* ( ...
    1 - ...
    1.93 * cos(2 * pi*fe/number_cycles * t) + ...
    1.29 * cos(4 * pi*fe/number_cycles * t) - ...
    0.388 * cos(6 * pi*fe/number_cycles * t) + ...
    0.028 * cos(8 * pi*fe/number_cycles * t));

force = signal(t_vector) .* (t_vector < final_time) + 0 .* (t_vector >= final_time);

% if (final_time > max(t_vector))
%     warning('Tone burst not fully developed. Increase simulation time vector');
% end

end

