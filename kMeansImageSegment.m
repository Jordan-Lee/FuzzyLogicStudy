% K means algorithm for Image segmentation

num_cent = 4; % assigning the number of centroids
num_dp = 100; % assigning the number of datapoints
%colorChart = ['b', 'r', 'g', 'y', 'c', 'm'];
points_type=['o' 'x' 's' '+' 'p' 'h','*'];

filename = 'wolf';
[data, mask_loc] = get_image(filename); %mask_loc is to check the error rate
[col, row] = size(data);
if filename=='wolf'
    num_cent = 3;
end
num_dp = col*row;

% Algorithm Testing

%img_kmean = k_means_1d(data, col, row, num_dp, num_cent);
img_3d_kmean = k_means_3d(data, col, row, num_dp, num_cent);

%figure;
%title('k_means_algorithm')
%imshow(img_kmean)

figure;
title('3d k means')
imshow(img_3d_kmean)


% accuracy Calculation


% K means Clustering
function [resultMat] = k_means_1d(data, col, row, num_dp, num_cent)
    dp = data(:);
    indices_min = 0;
    indices_max = 255;
    cent = indices_min + (indices_max-indices_min)*rand(num_cent, 2); % 2nd row is to store new value
    cluster_index = zeros(num_dp, 2); % 2nd row is to store color info
    mean = zeros(num_cent, 2);
    colorChart = ['b', 'r', 'g', 'y', 'c', 'm'];
    proceed = 1;
    iter_count = 0;
    
    distBuff = zeros(num_dp,num_cent);
    
    while proceed
        iter_count = iter_count+1;
    
        for i = 1:num_dp
            for j = 1:num_cent
                distBuff(i,j) = abs(cent(j, 1) - dp(i));
            end
        end
        
        % identifying the cluster of each datapoints 
        for i=1:num_dp
            dist_min = distBuff(i, 1);
            cluster_index(i, 1) = 1;
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
            mean(cluster_index(i, 1), 2) = mean(cluster_index(i, 1), 2) + 1; % count
        end
        for i=1:num_cent
            if mean(i,2) == 0
                cent(i,2) = 0;
            else    
               cent(i, 2) = mean(i, 1) / mean(i, 2); 
            end
        end
    
    % check continuity
    
    if round(cent(:,1)) == round(cent(:,2))
        proceed = 0; % terminate iteration
    end
    cent(:,1) = cent(:,2); % copying new centroid to 1st row
        
    end % while end
    
    %returning result in form of image
    temp = reshape(cluster_index(:,1), [col, row]);
    resultMat = temp./3;
    
end

function [resultMat] = k_means_3d(data, col, row, num_dp, num_cent)
    %performs K means clustering based on 3 dimensional data
    indices_min = 0;
    indices_max = 255;
    
    cent = zeros(num_cent, 3); % row1:col, row2:row, row3:indices
    cent(:,1) = col*rand(num_cent, 1); % assigning random value for column
    cent(:,2) = row*rand(num_cent, 1); % assigning random value for row
    cent(:,3) = indices_min + (indices_max-indices_min)*rand(num_cent,1); %assigning random value for indices
    cent_new = zeros(num_cent, 3); % matrix to hold newly calculated centroid coordinate
    cluster_index = zeros(num_dp, 2); % 2nd row is to store color info
    mean = zeros(num_cent, 4); % col, row, indices, count
    colorChart = ['b', 'r', 'g', 'y', 'c', 'm'];
    proceed = 1;
    iter_count = 0;
    
    distBuff = zeros(num_dp,num_cent);
    
    while proceed
        iter_count = iter_count + 1; % counting iteration
        
        %calculate the distance
        for i=1:num_cent
            for j=1:col
                for k=1:row
                    index = (j-1)*row + k;
                    distBuff(index, i) = sqrt( (j-cent(i, 1))^2 + (k-cent(i, 2))^2 + (data(j,k)-cent(i, 3))^2 );
                end
            end
        end
        % assign the cluster based on the distance
        for i = 1:num_dp
            min = distBuff(i,1);
            cluster_index(i, 1) = 1;
            for k=2:num_cent
                if(min > distBuff(i, k))
                   min = distBuff(i, k);
                   cluster_index(i, 1) = k;
                end
            end
            % assigning the color for datapoints to distinguish them
            cluster_index(i,2) = colorChart(cluster_index(i,1));
        end
        %reassign the centroid with the mean of datapoints
        for i=1:col
            for j=1:row
                index = (i-1)*row + j;
                index = cluster_index(index,1);
                mean(index,1) = mean(index,1) + i;
                mean(index,2) = mean(index,2) + j;
                mean(index,3) = mean(index,3) + data(i, j);
                mean(index,4) = mean(index,4) + 1;
            end
        end
        for i=1:num_cent
            if(mean(i,4) ~= 0)
                for j=1:3
                    cent_new(i, j) = mean(i,j) / mean(i,4);
                end
            end
        end
        
        %compare the coordinates and determine to finish or not
        if( round(cent) == round(cent_new) )
           proceed = 0;
        end
        cent = cent_new;
        
        %plot...
        figure;
        scatter3(cent(:,1), cent(:,2), cent(:,3),'filled', 'b','r','g');
        hold on
        for i=1:50
            for j=1:50
                scatter3(i,j,data(i, j), char(cluster_index((i-1)*row + j, 2)) )
            end
        end
        hold off 
        
    end % end of iteration
    
    % return result
    resultMat = reshape(cluster_index(:,1), [col, row] );
    
    
end

function [resultMat] = fuzzy_c_means_1st(data, col, row, num_dp)

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
    
    %figure
    %imshow(mask_loc);
    %figure;
    %imshow(img);
    
    mask = mask_loc;
    
end