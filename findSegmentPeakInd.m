function [alignmentInd] = findSegmentPeakInd(spike,segmentSigns,AlignmentSegment,takeAverage)
    % findSegmentPeakInd Description: For a custom shaped template where we
    % allow the user to select how spikes are aligned, we need a method to
    % recalculate the alignment index
    % ______________
    %[alignmentInd] = findSegmentPeakInd(spike,segmentSigns,AlignmentSegment,takeAverage)
    % Inputs:
    %   - spike: Ntime x 1 or Ntime x Nspike vector/matrix of spieks
    %   - segmentSigns: vector specifying template shape based on the expected positive negative eg. [1 -1] or [1 -1 1]
    %   - AlignmentSegment: index of segment to align on. (between 1 and
    %   length(segmentSigns). The extreme value (positive negative depneds
    %   on the sign of the segment) will be used to do alignment
    %   - takeAverage: true: If spike is a matrix we can output one value
    %   for the average of the set
    %                  false: we output a vector 1 x Nspike, for each
    %                  column
    % Outputs:
    %   - alignmentInd: either scalar index or [1 x Nspike] vector of
    %   indexes
    % ______________
    % Author: Kathrine Clarke 
    % Date: 30-Nov-2023

    % Output help if there are no inputs
    if nargin == 0
       help findSegmentPeakInd
       return;
    end

    if size(spike,2)>1 && takeAverage == true
        spike = nanmean(spike,2);
    end
    alignmentInd = zeros(1,size(spike,2));
    for i = 1:size(spike,2)
        alignmentInd(i) = do_alignment(spike(:,i),segmentSigns,AlignmentSegment);
    end
end
function alignmentInd = do_alignment(spike,segmentSigns,AlignmentSegment)
    %find the segments based on zero crossings
    zcs = [0; find(abs(diff(sign(spike)))>0);length(spike)];
    % remove zero crossings before we get into the first segment (this can
    % be going from zero to negative and then finaly to positive (if first
    % segment is positive
    zcs(1:find(spike(zcs(1:end-1)+1)*segmentSigns(1)>0,1)-1) = [];
    
    
    seg_inds = zcs(AlignmentSegment)+1:zcs(AlignmentSegment+1);
    [~,mind] = max(spike(seg_inds)*segmentSigns(AlignmentSegment));
    alignmentInd = mind + zcs(AlignmentSegment);
end