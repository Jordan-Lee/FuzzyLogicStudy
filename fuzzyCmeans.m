% Fuzzy C-means Algorithm implementation
% 
% 

clear all;

num_cent = 3;
num_dp = 200;

fuzzCoeff = 2;

dataPoints = randi([0 100],num_dp, 2);
centroids = zeros(num_cent, 2);
deg_of_mem = zeros(num_dp, num_cent);
maxMemb = [0,0];
proceed = 1;
countIter = 1;
THR = 0.000001;

colorChart = ['r', 'g', 'b', 'y', 'c', 'm'];
ws = 100;

% initializing deg_of_mem

for i=1:num_dp
   
    deg_of_mem(i,:) = rand(1, num_cent);
    norm = sum(deg_of_mem(i,:));
    deg_of_mem(i,:) = deg_of_mem(i,:) / norm;
    
end
maxMemb(1) = max(deg_of_mem(:));


while proceed
    countIter = countIter+1;
    %compute each centroids

    for cent=1:num_cent

        num = zeros(2);
        denom = 0;

        % numerator : sum of all deg_of_mem * datapoints
        for i = 1:num_dp

            num(1) = num(1) + (deg_of_mem(i, cent)^fuzzCoeff) * dataPoints(i, 1);
            num(2) = num(2) + (deg_of_mem(i, cent)^fuzzCoeff) * dataPoints(i, 2);

            denom = denom + deg_of_mem(i, cent)^fuzzCoeff;

        end
        centroids(cent,1) = num(1) / denom;
        centroids(cent,2) = num(2) / denom;
    end

    % compute degree of membership 

    for i = 1:num_dp
       for j = 1:num_cent

           dist_j = sqrt( (dataPoints(i,1) - centroids(j,1) )^2  + ( dataPoints(i,2) - centroids(j,2) )^2 );

           tempSum = 0;
           tempDist = 0;
           for k = 1:num_cent
              tempDist = sqrt( (dataPoints(i,1) - centroids(k,1) )^2  + ( dataPoints(i,2) - centroids(k,2) )^2 );
              tempSum = tempSum + (dist_j/tempDist)^(2/(fuzzCoeff-1)); 
           end

           deg_of_mem(i,j) = 1 / tempSum;

       end
    end
    
    maxMemb(2) = max(deg_of_mem(:));
    
    % convergence check
    if (abs(maxMemb(1) - maxMemb(2)) < THR)
        proceed = 0;
    end
    maxMemb(1) = maxMemb(2);
    
%     figure;
% 
%     set(gca, 'xlim',[0, ws], 'ylim', [0, ws])
%     title('randomly scattered datapoints & initial centroids');
%     scatter(dataPoints(i,1), dataPoints(i,2), 'MarkerEdgeColor', [deg_of_mem(1,1),deg_of_mem(1,2),deg_of_mem(1,3)])
% 
%     hold on
%     
%     for i=1:num_dp
%         scatter(dataPoints(i,1), dataPoints(i,2), 'MarkerEdgeColor', [deg_of_mem(i,1),deg_of_mem(i,2),deg_of_mem(i,3)]) 
%     end
%   
%     for i = 1:num_cent
%         scatter(centroids(i,1), centroids(i,2), colorChart(i) , 'filled');
%     end
%     hold off
    
    
    
end  
    
    
    
% 3 dimensional plot
    
figure;
    
set(gca, 'xlim',[0, ws], 'ylim', [0, ws])
title('randomly scattered datapoints & initial centroids');
scatter(dataPoints(i,1), dataPoints(i,2), 'MarkerEdgeColor', [deg_of_mem(1,1),deg_of_mem(1,2),deg_of_mem(1,3)])
    
hold on
    
    
    
for i=1:num_dp
    scatter(dataPoints(i,1), dataPoints(i,2), 'MarkerEdgeColor', [deg_of_mem(i,1),deg_of_mem(i,2),deg_of_mem(i,3)]) 
end

for i = 1:num_cent
    scatter(centroids(i,1), centroids(i,2), colorChart(i) , 'filled');
end
hold off
