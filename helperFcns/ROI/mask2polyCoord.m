function ROIxyCoord = mask2polyCoord(mask)
%takes binary image mask of shape imgWidth X imgHeight of type logical for
%a given ROI and outputs the coordinates of the ROI outline in ROIxyCoord
%Row1: x; Row2: y

[r,c] = ind2sub(size(mask),find(mask));
[uC,~,ibC] = unique(c);

[xC,yC] = deal(zeros(1,length(uC)));
iPt = 0;
for k = 1:length(uC)
    iPt = iPt+1;
    xC(iPt) = uC(k);
    yC(iPt) = min(r(ibC==k));
    
    iPt = iPt+1;
    xC(iPt) = uC(k);
    yC(iPt) = max(r(ibC==k));
end

[uR,~,ibR] = unique(r);
[xR,yR] = deal(zeros(1,length(uR)));
iPt = 0;
for k = 1:length(uR)
    iPt = iPt+1;
    
    xR(iPt) = min(c(ibR==k));
    yR(iPt) = uR(k);
    
    iPt = iPt+1;
    xR(iPt) = max(c(ibR==k));
    yR(iPt) = uR(k);
end
cXY = [xC;yC]';
rXY = [xR;yR]';
ROIxyCoord = unique([rXY;cXY],'rows')';
end