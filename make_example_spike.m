function [] = make_example_spike(Nblocks,First_sign,align_block,voltage_thresholds,duration_phases)
    % make_example_spike Description: This function makes a little plot,
    % that demonstrates how spike shape works. It should be called in the
    % set parameters section for IdentifyAP templates
    % ______________
    %[] = make_example_spike(Nblocks,First_sign,align_block)
    % Inputs:
    %   - Nblocks: number of blocks of positive negative bits of signal
    %   - First_sign: -1 or 1 <- the sign we're looking for in the first
    %   block
    %   - align_block: integer between 1 and Nblocks - align spikes based
    %   on the extrema in this block
    % Outputs:
    %   - :
    % ______________
    % Author: Kathrine Clarke 
    % Date: 30-Nov-2023

    % Output help if there are no inputs
    if nargin == 0
       help make_example_spike
       return;
    end
    %spike 1
    Blength = ceil(min(voltage_thresholds))*100;
    bl_array = ceil(Blength.*[1 duration_phases 1]);

    Nt = sum(bl_array);
    fill([-Nt*2 Nt*2 Nt*2 -Nt*2  -Nt*2],[ -1 -1 1 1 -1],[0.8 0.8 0.8],'EdgeAlpha',0)
    hold on
    
    for i = 1:5
        bl_array = ceil(Blength.*[1 duration_phases 1]);
        bl_array = bl_array + randi(ceil(Blength/2),size(bl_array));
        dummy_spike(Nblocks,First_sign,align_block,bl_array,voltage_thresholds,duration_phases,0,[0.5 0.5 0.5]);
    end
    bl_array = ceil(Blength.*[1 duration_phases 1]);
    bl_array = bl_array + randi(ceil(Blength/3),size(bl_array)); 
    [block_signs,xl] = plot_1_spike(Nblocks,First_sign,align_block,bl_array,voltage_thresholds,duration_phases,0,'r');

    xticks(0)
    xticklabels('Aligned')
    yticks(unique(sort([block_signs(2:end-1).*voltage_thresholds';0;-1;1])))
    ylabel('voltage (x std noise)')
    xlim(xl)
    title({sprintf('%d Phase spike,',Nblocks), sprintf('starting with %s phase,',ternaryOp(block_signs(2)==1,'positive','negative')),sprintf('aligned based on phase %d',align_block)})
    hold off
end
function  dummy_spike(Nblocks,First_sign,align_block,bl_array,voltage_thresholds,duration_phases,plot_offset,col)
    redo = true;
    while redo
        Nt = sum(bl_array);
        
        %make the heights of the block signs be voltage_thresholds
        block_signs = First_sign*cumprod(-1*ones(Nblocks+2,1));
        pvals = [block_signs(1)*0.01 block_signs(1)*0.01 toVec([zeros(1,Nblocks); block_signs(2:end-1)'.*(rand(1,Nblocks)+voltage_thresholds)])' 0 block_signs(end)*0.01 block_signs(end)*0.01];
        rand_pos = arrayfun(@(b) 1+randi(bl_array(b)),2:Nblocks+1);

        pinds = [0 ceil(bl_array(1)/2) toVec([zeros(1,Nblocks);rand_pos]+cumsum(bl_array(1:end-2)))' sum(bl_array(1:end-1)) Nt-ceil(bl_array(end)/2) Nt];
        pinds = pinds + linspace(0,1,length(pinds));
        sp = spline(pinds,pvals,1:Nt);
        block_signs([1 end]) = nan;
        voltage_thresholds2 = [0 voltage_thresholds 0];
        blocksignal = arrayfun(@(b) {[0;voltage_thresholds2(b).*block_signs(b).*ones(bl_array(b)-1,1)]} ,1:(Nblocks+2));
        blocksignal = vertcat(blocksignal{:});
        if ~any((sign(blocksignal)'.*sign(sp))==-1)
            redo = false;
        end
         sp(isnan(blocksignal)|blocksignal==0) = 0.5*sp(isnan(blocksignal)|blocksignal==0)./max(sp(isnan(blocksignal)|blocksignal==0));
   
        if nanmax(abs(sp))>max(voltage_thresholds)*2
            redo = true;
        end
    end
    %sp = sp./nanmax(abs(sp));
    sp_a = zeros(size(sp));
    block_inds = @(i) sum(bl_array(1:i-1))+1:sum(bl_array(1:i));
    sp_a(block_inds(align_block+1)) = abs(sp(block_inds(align_block+1)) );
    [m,a_ind] = max(sp_a);
    hold on
    plot((1:length(sp))-a_ind,sp+plot_offset,':','Color',[0.5 0.5 0.5])
end
function [block_signs,xl] = plot_1_spike(Nblocks,First_sign,align_block,bl_array,voltage_thresholds,duration_phases,plot_offset,col)
    redo = true;
    while redo
        Nt = sum(bl_array);
        
        %make the heights of the block signs be voltage_thresholds
        block_signs = First_sign*cumprod(-1*ones(Nblocks+2,1));
        pvals = [block_signs(1)*0.01 block_signs(1)*0.01 toVec([zeros(1,Nblocks); block_signs(2:end-1)'.*(rand(1,Nblocks)+voltage_thresholds)])' 0 block_signs(end)*0.01 block_signs(end)*0.01];
        rand_pos = arrayfun(@(b) 1+randi(bl_array(b)),2:Nblocks+1);

        pinds = [0 ceil(bl_array(1)/2) toVec([zeros(1,Nblocks);rand_pos]+cumsum(bl_array(1:end-2)))' sum(bl_array(1:end-1)) Nt-ceil(bl_array(end)/2) Nt];
        pinds = pinds + linspace(0,1,length(pinds));
        sp = spline(pinds,pvals,1:Nt);
        block_signs([1 end]) = nan;
        voltage_thresholds2 = [0 voltage_thresholds 0];
        blocksignal = arrayfun(@(b) {[0;voltage_thresholds2(b).*block_signs(b).*ones(bl_array(b)-1,1)]} ,1:(Nblocks+2));
        blocksignal = vertcat(blocksignal{:});
        sp(isnan(blocksignal)|blocksignal==0) = 0.5*sp(isnan(blocksignal)|blocksignal==0)./max(sp(isnan(blocksignal)|blocksignal==0));
   
        if ~any((sign(blocksignal)'.*sign(sp))==-1)
            redo = false;
        end
        if nanmax(abs(sp))>max(voltage_thresholds)*2
            redo = true;
        end
    end
        %sp = sp./nanmax(abs(sp));
    sp_a = zeros(size(sp));
    block_inds = @(i) sum(bl_array(1:i-1))+1:sum(bl_array(1:i));
    block_edges = @(i) [sum(bl_array(1:i-1))+1 sum(bl_array(1:i))];
    sp_a(block_inds(align_block+1)) = abs(sp(block_inds(align_block+1)) );
    [m,a_ind] = max(sp_a);
    
    plot([-a_ind length(sp)-a_ind],[0 0]+plot_offset,'k')
    hold on
    plot((1:length(sp))-a_ind,sp+plot_offset,'Color',col)
    plot((1:length(blocksignal))-a_ind,blocksignal+plot_offset,'k','lineWidth',3)
    plot([0 0],[0 m.*block_signs(align_block+1)]+plot_offset,'g','lineWidth',3)
    arrayfun(@(b) plot(-a_ind+block_edges(b),(voltage_thresholds2(b)+0.3).*block_signs(b).*[1 1],'k') ,2:(Nblocks+1));
    arrayfun(@(b) text(-a_ind+sum(bl_array(1:b-1)),(voltage_thresholds2(b)+0.5).*block_signs(b).*[1],sprintf('>%0.2f ms',duration_phases(b-1))) ,2:(Nblocks+1));
   
    xl = [-a_ind, length(sp)-a_ind];
end