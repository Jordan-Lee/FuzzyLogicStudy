% K means algorithm for Image segmentation

num_cent = 4; % assigning the number of centroids
num_dp = 100; % assigning the number of datapoints
colorChart = ['b', 'r', 'g', 'y', 'c', 'm'];
points_type=['o' 'x' 's' '+' 'p' 'h','*'];

filename = 'wolf';
[data, mask_loc] = get_image(filename); %mask_loc is to check the error rate
[col, row] = size(data);
if filename=='wolf'
    num_cent = 4;
end
num_dp = col*row;



% K means Clustering
function [resultMat] = k_means_1d(data, col, row, num_dp. num_cent)
    dp = data(:);
    indices_min = 0;
    indices_max = 255;
    cent = rand([0 255], num_cent, 1);
    cluster_index = zeros(num_dp, 2);
    mean = zeros(num_cent, 2);
    
    distBuff = zeros(num_dp,num_cent);
    
    
    for i = 1:num_dp
        for j = 1:num_dp
            distBuff(i,j) = abs(cent(j) - data(i));
        end
    end
    % identifying the cluster of each datapoints 
    for i=1:num_dp
        dist_min = distBuff(i, 1);
        for j=1:num_cent
            if dist_min > distBuff(i, j)
                cluster_index(i, 1) = j;
            end
        end
        
        % assigning the color for datapoints to distinguish them
        cluster_index(i,2) = colorChart(cluster_index(i,1));
        
    end
    % relocating the centroids
    for i=1:num_dp
        mean(cluster_index(i, 1), 1) = mean(cluster_index(i,1), 1) + dp(i); % adding datapoints to their 
        mean(cluster_index(i, 1), 2) = mean(cluste_index(i, 1), 2) + 1; % count
    end
    for i=1:num_cent
        cent(i) = mean(i, 1) / mean(i, 2);
    end
    


end

function [indices, mask] = get_image(filename)

    IMG = imread(strcat(filename, '.jpg' ));
    img = rgb2gray(IMG);
    indices = double(img);

    MSK = imread(strcat(strcat(filename,'_mask' ),'.jpg'));
    msk = rgb2gray(MSK);
    msk = double(msk);
    mask_loc = zeros(size(msk));
    % figure;
    % imshow(msk)
    if filename == 'wolf'
        
        mask_loc(find(100>msk))=0;
        mask_loc(find(msk>200))=1;
        mask_loc(find(((200>=msk).*msk)>=100)) = 0.5;
        
    end
    
    figure
    imshow(mask_loc);
    figure;
    imshow(img);
    
    mask = mask_loc;
    
end