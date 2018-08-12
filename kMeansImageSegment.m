% K means algorithm for Image segmentation

num_cent = 4; % assigning the number of centroids
num_dp = 100; % assigning the number of datapoints
colorChart = ['b', 'r', 'g', 'y', 'c', 'm'];
points_type=['o' 'x' 's' '+' 'p' 'h','*'];

filename = 'wolf';



function [indices, mask] = get_image(filename);

IMG = imread(strcat(filename, '.jpg' ));
img = rgb2gray(IMG);
indices = double(img);

msk = imread(strcat(strcat(filename,'_mask' ),'.jpg'));
    if filename == 'wolf'
        msk(find(100>msk)) = 0;
    end

end