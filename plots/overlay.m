function overlay(modulation_label, modulation_frequency, alpha_channel)

if (nargin == 2)
    alpha_channel = 0.2;
end

label = sprintf('%s %+.0f', modulation_label, abs(modulation_frequency));
base_path = 'plots/screenshots/';

switch (label)

    case {'ON-ON-ON +0' 'C- +0'}
        imagesc([-5 0], [0 20], flipud(imread([base_path 'ON.png'])));

    case {'OFF-OFF-OFF +0' 'OFF +0'}
        imagesc([-5 0], [0 20], flipud(imread([base_path 'OFF.png'])));

    case {'ON-OFF-OFF +0', 'Sinusoidal (discrete) +0'}
        imagesc([-5 0], [0 20], flipud(imread([base_path 'freezed.png'])));

    case 'Sinusoidal (discrete) +1000'
        imagesc([-5 -3.5], [6 13], flipud(imread([base_path 'Sinusoidal (discrete) +1000.png'])));
        imagesc([3.25 4.75], [6 13], flipud(imread([base_path 'Sinusoidal (discrete) -1000.png'])));

    case 'Sinusoidal (discrete) +2000'
        imagesc([-5 -3.5], [6 13], flipud(imread([base_path 'Sinusoidal (discrete) +2000.png'])));
        imagesc([3 4.5], [6 13], flipud(imread([base_path 'Sinusoidal (discrete) -2000.png'])));

    case 'Sinusoidal (discrete) +3000'
        imagesc([-5 -3.5], [6 13], flipud(imread([base_path 'Sinusoidal (discrete) +3000.png'])));
        imagesc([3 4.5], [6 13], flipud(imread([base_path 'Sinusoidal (discrete) -3000.png'])));

    otherwise
        return
end

image_objs = findobj(gca, 'Type', 'image');
if ~isempty(image_objs)
    for image_obj = image_objs(:)'
        alphaData = ones(size(image_obj, 1), size(image_obj, 2)) * alpha_channel;
        set(image_obj, 'AlphaData', alphaData);
    end
end

set(gca, 'YDir', 'normal')

end