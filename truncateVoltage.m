% [voltage, time] = truncateVoltage(voltage, time)
% 
% Prompt user to select either time range for voltage timeseries, or number
% of samples to include, & then truncate voltage & time before returning
function [voltage, time] = truncateVoltage(voltage, time)
   def{1,1} = 2;
   def{2,1} = 0;
   def{3,1} = length(voltage);
   
   prompt{1,1} = 'Truncate voltage by time or number of samples';
   prompt{1,2} = 'type';

   prompt{2,1} = 'Min';
   prompt{2,2} = 'min';

   prompt{3,1} = 'Max';
   prompt{3,2} = 'max';


    % fields in format structure that you can fuck with
    s1 = struct('type','list', 'style','popupmenu', 'size',-1, 'format','integer');
    s2 = struct('type','edit', 'format','float', 'style','edit', 'size', -1);
    formats = [s1; s2; s2];
    formats(1).items  = {'Time   ','Samples'};
    formats(1).limits = [1 2];
    formats(2).limits = [0 length(time)];
    formats(3).items  = [0 length(time)];
    
    % default answers
    [truncate,cancel] = inputsdlg(prompt, 'Truncate voltage', formats, def);
    if cancel
        return;
    else
       maxT = time(end);
       maxN = length(time);
       switch truncate.type % index value
          case 1 % time
             dt   = time(2) - time(1);
             sInd = ternaryOp(0<truncate.min && truncate.min<maxT, round(truncate.min / dt), 1);
             eInd = ternaryOp(0<truncate.max && truncate.max<maxT, round(truncate.max / dt), maxN);
          case 2 % samples
             sInd = ternaryOp(0<truncate.min && truncate.min<maxN, round(truncate.min), 1);
             eInd = ternaryOp(0<truncate.max && truncate.max<maxN, round(truncate.max), maxN);
       end
       voltage = voltage(sInd:eInd);
       time    = time(sInd:eInd);
    end
end