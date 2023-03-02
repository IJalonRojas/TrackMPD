function [triangleIndices,barycentric] = pointLocationQuadTree(tri,queryPoints,n)
% pointLocationQuadTree  Triangle containing specified points
%   This function mimics the behavior of triangulation/pointLocation for 2d
%   triangulations using a quadtree search to accelerate the search for
%   many queryPoints in a large triangulation. If the triangulation is a
%   delaunayTriangulation, it directly uses
%   delaunayTriangulation/pointLocation, which is already very fast.
%
% See also triangulation/pointLocation,
% delaunayTriangulation/pointLocation.
%
% Copyright (C) Sebastian Ullmann 2018
  switch nargin
    
    case 3
      
      % Do nothing.
      
    case 2
      
      % Set default value. The algorithm terminates for any n, but a
      % reasonable value is needed for efficiency. Depending on the
      % machine, n = 1e6 may be a reasonable value.
      %
      % The algorithm does a standard search with pointLocationStandard()
      % as soon as the number of query points in the cluster times the
      % number of triangles in the cluster is smaller then n.
      %
      % If n == inf, then the algorithm directly does a standard search
      % with pointLocationStandard().
      %
      % If n <= 1, then the algorithm clusters isolated (multiple) query
      % points and then does a standard search with pointLocationStandard()
      % for each distinct query point.
      
      n = 1e6;
      
    otherwise
      
      error('Wrong number of input arguments');
      
  end
  
  if isa(tri,'delaunayTriangulation')
    [triangleIndices,barycentric] = tri.pointLocation(queryPoints);
    return
  end
  
  if size(tri.ConnectivityList,1) * size(queryPoints,1) < n
    [triangleIndices,barycentric] = ...
      pointLocationStandard(tri,queryPoints);
    return
  end
  % We define the clusters according to the location of the query points.
  % This ensures that in the case of n==1, we end up with isolated
  % (multiple) query points, which is easy to detect. What is a good
  % termination criterion if we cluster according to the location of the
  % simplices?
  
  xmedian = median(queryPoints(:,1),1);
  ymedian = median(queryPoints(:,2),1);
  
  isEast  = @(x) x(:,1) >= xmedian;
  isNorth = @(x) x(:,2) >= ymedian;
  
  isQueryPointE = isEast( queryPoints);
  isQueryPointN = isNorth(queryPoints);
  
  isQueryPointInQuad = ...
    [( isQueryPointN &  isQueryPointE), ...
     (~isQueryPointN &  isQueryPointE), ...
     (~isQueryPointN & ~isQueryPointE), ...
     ( isQueryPointN & ~isQueryPointE)];
                  
  if sum(any(isQueryPointInQuad))<2
    
    % All query points are in one cluster. This means all query points are 
    % the same.
    
    assert(all(all(queryPoints(1,:)==queryPoints)))
    
    [triangleIndices,barycentric] = pointLocationStandard(tri,queryPoints(1,:));
    
    triangleIndices = repmat(triangleIndices,size(queryPoints,1),1);
    barycentric = repmat(barycentric,size(queryPoints,1),1);
    
    return
    
  end
  
  points = tri.Points;
  triangles = tri.ConnectivityList;
  
  A = triangles(:,1);
  B = triangles(:,2);
  C = triangles(:,3);
  
  isTriPointE = isEast( points);
  isTriPointW = ~isTriPointE;
  isTriPointN = isNorth(points);
  isTriPointS = ~isTriPointN;
  
  isTriangleE =  isTriPointE(A) | isTriPointE(B) | isTriPointE(C);
  isTriangleW =  isTriPointW(A) | isTriPointW(B) | isTriPointW(C);
  isTriangleN =  isTriPointN(A) | isTriPointN(B) | isTriPointN(C);
  isTriangleS =  isTriPointS(A) | isTriPointS(B) | isTriPointS(C);
    
  isTriangleInQuad = ...
    [(isTriangleN & isTriangleE), ...
     (isTriangleS & isTriangleE), ...
     (isTriangleS & isTriangleW), ...
     (isTriangleN & isTriangleW)];
   
  triangleIndices = nan(size(queryPoints,1),1);
  barycentric = nan(size(queryPoints,1),size(triangles,2));  
  
  for i = find(any(isQueryPointInQuad) & any(isTriangleInQuad))
    
    % We do not need to perform any computation if there is not any
    % query point in the i-th quad, because the result will be empty.
    %
    % We do not need to perform any computation if there is not any
    % triangle in the i-th quad, because the result will be NaNs.
    
    indTrianglesInQuad = find(isTriangleInQuad(:,i));
    
    % We create a local triangulation by picking selected triangles from the
    % global triangulation. Therefore we need a local point list and a local
    % triangle list pointing to the indices of the local point list.
    indPointsInQuad = unique(triangles(indTrianglesInQuad,:));
    localPoints = points(indPointsInQuad,:);
    localPointMap = nan(size(points,1),1);
    localPointMap(indPointsInQuad,:) = 1:length(indPointsInQuad);
    localTriangles = localPointMap(triangles(indTrianglesInQuad,:));
    localTriangles = reshape(localTriangles,length(indTrianglesInQuad),[]);
    localTriangulation = struct('ConnectivityList',localTriangles, ...
                                'Points',localPoints);
    % We perform a tree search in the local triangulation.
    
    localQueryPoints = queryPoints(isQueryPointInQuad(:,i),:);
    [localTriangleIndices,outputBarycentric] = pointLocationQuadTree(localTriangulation,localQueryPoints,n);
    % We provide the output triangle indices in terms of indices into the
    % point list of the global triangulation.
    outputTriangleIndices = nan(size(localTriangleIndices));
    isNotNan = ~isnan(localTriangleIndices);
    outputTriangleIndices(isNotNan) = indTrianglesInQuad(localTriangleIndices(isNotNan));
    
    triangleIndices(isQueryPointInQuad(:,i),:) = outputTriangleIndices;
    barycentric(isQueryPointInQuad(:,i),:) = outputBarycentric;
    
  end
end
function [triangleIndices,barycentric] = pointLocationStandard(tri,queryPoints)
  triangles = tri.ConnectivityList;
  points = tri.Points;
  A = triangles(:,1);
  B = triangles(:,2);
  C = triangles(:,3);
  xA = points(A,1);
  yA = points(A,2);
  xB = points(B,1);
  yB = points(B,2);
  xC = points(C,1);
  yC = points(C,2);
  inverseDeterminant = 1 ./ ((xB-xA).*(yC-yA) - (yB-yA).*(xC-xA));
  triangleIndices = nan(size(queryPoints,1),1);
  barycentric = nan(size(queryPoints,1),3);
  for i=1:size(queryPoints,1)
    xQ = queryPoints(i,1);
    yQ = queryPoints(i,2);
    barycentricA = ((xB-xQ).*(yC-yQ) - (yB-yQ).*(xC-xQ)) .* inverseDeterminant;
    barycentricB = ((xC-xQ).*(yA-yQ) - (yC-yQ).*(xA-xQ)) .* inverseDeterminant;
    barycentricC = ((xA-xQ).*(yB-yQ) - (yA-yQ).*(xB-xQ)) .* inverseDeterminant;
    
    isInside = barycentricA>=0 & barycentricB>=0 & barycentricC>=0;
    
    ind = find(isInside,1);
    
    if isempty(ind)  % Found no triangle containing the query point.
      
       % Find the closest triangle in terms of barycentric coordinates.
      [negativeDistance,ind] = max(min([barycentricA barycentricB barycentricC],[],2));
      % The tolerance should be large enough to have enough "overlap"
      % between the triangles to catch points which are "between" triangles
      % due to roundoff in the barycentric coordinates. The tolerance
      % should be small enough to prevent points outside the triangulation
      % to be falsely assigned to a triangle -- we would like to return NaN
      % for outliers. Choose a larger tolerance in case of heavily
      % distorted triangles or if falsely assigned outliers are not a
      % problem.
      
      distance = -negativeDistance;
      triangleLength = sqrt(1./abs(inverseDeterminant(ind)));
      tolerance = 1000*eps(triangleLength);
      
      if distance > tolerance
        
        % Leave NaNs in the output arrays.
        continue
        
      end
      
    end
    
    triangleIndices(i) = ind;
    barycentric(i,:) = [barycentricA(ind) barycentricB(ind) barycentricC(ind)];
    
  end
end