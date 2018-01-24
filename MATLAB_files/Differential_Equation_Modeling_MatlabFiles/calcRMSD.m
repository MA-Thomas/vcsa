function r = calcRMSD(simSAXSexp,realSAXSexp)
% input: simSAXSexp is curves(:,idx,:) - each col vector is Iq curve.
% realSAXSexp is expCurve{1}(:,:) - eac col vector is Iq curve.
size(simSAXSexp,3)
size(realSAXSexp,2)
assert(isequal( size(simSAXSexp,3),size(realSAXSexp,3) ));
SE = zeros(1,size(simSAXSexp,3));

for u = 1:size(simSAXSexp,3)
    diff = simSAXSexp(:,1,u) - realSAXSexp(:,1,u);
    %SE(u) = diff'*diff/length(diff); %MT: This division is part of the "mean" in RMSD. IT MUST COME AFTER both sums which comprise the mean, not btw them.
    SE(u) = diff'*diff;
end

r = (sum(SE)/numel(realSAXSexp))^0.5;
end