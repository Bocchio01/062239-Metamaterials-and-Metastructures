function overlay(file_name, frequency)

bgImage = imread(['plots/' file_name '.png']);

imagesc([-5 0], [0 20], flipud(bgImage));
alphaData = ones(size(bgImage, 1), size(bgImage, 2)) * 0.2;

set(findobj(gca, 'Type', 'Image'), 'AlphaData', alphaData);
set(gca, 'YDir', 'normal')

end