% K means algorithm implementation


num_cent = 5;
num_dp = 100;
ws = 100;

colorChart = ['b', 'r', 'g', 'y', 'c', 'm']; 

dataPoints = randi([0 100], num_dp, 2);
centroids = randi([0 100], num_cent, 2);
clusterIndex = zeros(num_dp, 2);

% inital plot 
figure;

set(gca, 'xlim',[0, ws], 'ylim', [0, ws])
title('randomly scattered datapoints & initial centroids');
scatter(dataPoints(:,1), dataPoints(:,2), 'black') 

hold on
for i = 1:num_cent
    scatter(centroids(i,1), centroids(i,2), colorChart(i) , 'filled');
end
hold off


% clustering starts

distBuff = zeros(num_dp, num_cent);
% distance calculation
for i = 1:num_dp
    for j = 1:num_cent
    
        distBuff(i, j) = sqrt( (dataPoints(i, 1) - centroids(j,1))^2 + (dataPoints(i, 2) - centroids(j,2))^2 );
        
    end
end

% cluster back

for i=1:num_dp
    
    minimum = distBuff(i, 1);
    clusterIndex(i) = 1;
    
   for j = 2:num_cent
       
       if ( distBuff(i, j) < minimum)
           minimum = distBuff(i, j);
           clusterIndex(i, 1) = j;
       end
       
   end
   
       clusterIndex(i, 2) = colorChart(clusterIndex(i, 1)); 
end



figure;

set(gca, 'xlim',[0, ws], 'ylim', [0, ws])
title('initial custering');
scatter(dataPoints(1,1), dataPoints(1,2), char(clusterIndex(1,2)))
hold on
for i = 2:num_dp
    scatter(dataPoints(i,1), dataPoints(i,2), char(clusterIndex(i,2)))
end


for i = 1:num_cent
    scatter(centroids(i,1), centroids(i,2), colorChart(i) , 'filled');
end
hold off



% moving centroid

mean = zeros(num_cent, 3);
for i=1:num_dp
        mean(clusterIndex(i, 1), 1) = mean(clusterIndex(i, 1), 1) + dataPoints(i, 1);
        mean(clusterIndex(i, 1), 2) = mean(clusterIndex(i, 1), 2) + dataPoints(i, 2);
        mean(clusterIndex(i, 1), 3) = mean(clusterIndex(i, 1), 3) + 1;
end
for i=1:num_cent
   
    centroids(i, 1) = mean(i, 1) / mean(i, 3);
    centroids(i, 2) = mean(i, 2) / mean(i, 3);
    
end

figure;

set(gca, 'xlim',[0, ws], 'ylim', [0, ws])
title('modified centroids');
scatter(dataPoints(1,1), dataPoints(1,2), char(clusterIndex(1,2)))
hold on
for i = 2:num_dp
    scatter(dataPoints(i,1), dataPoints(i,2), char(clusterIndex(i,2)))
end


for i = 1:num_cent
    scatter(centroids(i,1), centroids(i,2), colorChart(i) , 'filled')
end
hold off


prevCentroids = zeros(num_cent, 2);

%iternation

while prevCentroids ~= centroids
    
    prevCentroids = centroids;
   
       distBuff = zeros(num_dp, num_cent);
    % distance calculation
    for i = 1:num_dp
        for j = 1:num_cent

            distBuff(i, j) = sqrt( (dataPoints(i, 1) - centroids(j,1))^2 + (dataPoints(i, 2) - centroids(j,2))^2 );

        end
    end

    % cluster back

    for i=1:num_dp

        minimum = distBuff(i, 1);
        clusterIndex(i) = 1;

       for j = 2:num_cent

           if ( distBuff(i, j) < minimum)
               minimum = distBuff(i, j);
               clusterIndex(i, 1) = j;
           end

       end

           clusterIndex(i, 2) = colorChart(clusterIndex(i, 1)); 
    end

    mean = zeros(num_cent, 3);
    for i=1:num_dp
            mean(clusterIndex(i, 1), 1) = mean(clusterIndex(i, 1), 1) + dataPoints(i, 1);
            mean(clusterIndex(i, 1), 2) = mean(clusterIndex(i, 1), 2) + dataPoints(i, 2);
            mean(clusterIndex(i, 1), 3) = mean(clusterIndex(i, 1), 3) + 1;
    end
    for i=1:num_cent

        centroids(i, 1) = mean(i, 1) / mean(i, 3);
        centroids(i, 2) = mean(i, 2) / mean(i, 3);

    end
    
end

figure;

set(gca, 'xlim',[0, ws], 'ylim', [0, ws])
title('10th iteration');
scatter(dataPoints(1,1), dataPoints(1,2), char(clusterIndex(1,2)))
hold on
for i = 2:num_dp
    scatter(dataPoints(i,1), dataPoints(i,2), char(clusterIndex(i,2)))
end


for i = 1:num_cent
    scatter(centroids(i,1), centroids(i,2), colorChart(i) , 'filled')
end
hold off