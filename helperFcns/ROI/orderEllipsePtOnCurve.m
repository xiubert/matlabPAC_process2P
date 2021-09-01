%order coordinates of ellipse along ellipse curve
%IMPORTANT NOTE: will likely only work on smooth ellipses without folds,
%otherwise the problem becomes a nearest neighbor problem like the traveling salesman

%based off of: https://www.mathworks.com/matlabcentral/answers/431568-how-to-effectively-reorder-data
function curveOrderedPtsXY = orderEllipsePtOnCurve(input2byPtNo)

a = atan2(input2byPtNo(1,:) - mean(input2byPtNo(1,:)),...
    input2byPtNo(2,:) - mean(input2byPtNo(2,:)));
[~,order] = sort(a);
reorderXY = [input2byPtNo(1,order); input2byPtNo(2,order)];

if ~all(reorderXY(:,end)==reorderXY(:,1))
reorderXY(:,numel(order)+1) = reorderXY(:,1);
end
curveOrderedPtsXY = reorderXY;
% figure;plot(reorderXY(1,:),reorderXY(2,:))
% hold on
% plot(input2byPtNo(1,:),input2byPtNo(2,:),'r.')

