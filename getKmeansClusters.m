% GETKMEANSCLUSTERS Classifies identified spikes into different clusters.
%
% Syntax: spclusters = getKmeansClusters(spikes, stimes);
%
% Inputs:
%     spikes - M x N matrix of N spikes with M samples per spike
%     stimes - Array of N cells, each containing M time values for the
%              samples in spikes.
% Outputs:
%     spclusters - Clusters
%
% Created by Artemio - 24/June/2019
function varargout = getKmeansClusters(spikes,stimes,varargin)
   if nargin > 2, debug = varargin{1}; else, debug = 'none'; end
   if nargin > 3
      cluster_range = varargin{2};
      min_k = cluster_range(1);
      max_k = cluster_range(2);
      extracting = 'spikes';
   else
      max_k = 30;
      min_k = 5;
      extracting = 'families';
   end
   NS = size(spikes,2);
%    pos_peak = zeros(1,NS);
%    neg_peak = zeros(1,NS);
%    pos_duration = zeros(1,NS);
%    neg_duration = zeros(1,NS);
%    diff_peaks = zeros(1,NS); % Distance between positive and negative peaks
   % Display message of progress
   str = sprintf('\tFinding cluster classification of %d spikes...\n',NS);
   printMessage('off','Keywords',str);
   for i = 1:NS
      [pos_peak(i), posidx] = max(spikes(:,i));
      [neg_peak(i), negidx] = min(spikes(:,i));
      diff_peaks(i) = abs(negidx - posidx);
      pos_duration(i) =  sum(diff(find(spikes(:,i)>0))==1);
      neg_duration(i) = sum(diff(find(spikes(:,i)<0))==1);
   end
   % Principal component analysis
   [~,score,l] = pca(spikes'); % Principal component analysis
   score = score(:,1:5); % Get only the projection on the highest latent component
  
   sse = zeros(max_k,1);
   idx = zeros(max_k,NS);
   % Find the clusters with a range of initial number of clusters (k)
   for k = min_k:max_k
      sse(k) = 0;
%       [idx(k,:), c] = kmeans([pos_peak;neg_peak;diff_peaks;pos_duration;neg_duration;score]',k); % K-means with k clusters
      [idx(k,:), c, e, ~] = kmeans(score(:,1:5),k, 'EmptyAction','drop','MaxIter',200); % K-means with k clusters
      sse(k) = sum(e.^2);
   end
   % Find the elbow
   elbow = find(abs(diff(diff(sse(min_k:max_k)))) < mean(abs(diff(diff(sse(min_k:max_k)))))/10, 1, 'first') + 4;
   if ~strcmp('none',debug)
      figure('Name','K-means elbow test'); plot(min_k:max_k,sse(min_k:max_k),'-o');hold('on');
      plot(elbow, sse(elbow),'xr', 'LineWidth', 2, 'MarkerSize', 10);
      xlabel('K');
      ylabel('SSE');
   end
   
   % Get the clusters with the best k value.
   idx = idx(elbow,:);
   APspikes = cell(1,elbow);
   APtimes = cell(1,elbow);
   APtemplates = zeros(size(spikes,1),elbow);
   if strcmp('full',debug), ff = figure('Name','Templates');end
   for i = 1:elbow
      APspikes(i) = cell({spikes(:,idx==i)});
      APtimes(i) = cell({stimes(:,idx==i)});
      APtemplates(:,i) = mean(spikes(:,idx==i),2)';
      if strcmp('full',debug)
         % If DEBUG, plot the first 10 spikes found for each of the
         % clusters and the cluster mean, i.e. the template.
         figure(ff);subplot(ceil(sqrt(elbow)),ceil(sqrt(elbow)),i);
         plot(APspikes{i}(:,1:min(end,10)),'--r'); hold('on');
         plot(APtemplates(:,i),'LineWidth',2);
      end
   end
   
   if strcmp('full',debug)
      % If DEBUG si full, plot the first 5 clusters in 3D displaying the
      % first 3 principal components.
      figure; plot3(score(idx(1,:)==1,1), score(idx(1,:)==1,2), score(idx(1,:)==1,3),'o','LineWidth',2);hold('on');grid;
      plot3(score(idx(1,:)==2,1), score(idx(1,:)==2,2), score(idx(1,:)==2,3),'or','LineWidth',2);
      plot3(score(idx(1,:)==3,1), score(idx(1,:)==3,2), score(idx(1,:)==3,3),'og','LineWidth',2);
      plot3(score(idx(1,:)==4,1), score(idx(1,:)==4,2), score(idx(1,:)==4,3),'om','LineWidth',2);
      plot3(score(idx(1,:)==5,1), score(idx(1,:)==5,2), score(idx(1,:)==5,3),'oy','LineWidth',2);
      xlabel('X');ylabel('Y');zlabel('Z');
   end
   
   switch extracting
      case 'APS'
         varargout = {APspikes, APtemplates};
      case 'spikes'
         varargout = {APspikes, APtimes};
      otherwise
         error('Invalid option. This function either extracts spikes or AP families.');
   end
   
   % Print 'Done' message
   str = sprintf('\tDone!\n');
   printMessage('off','Keywords',str);
end