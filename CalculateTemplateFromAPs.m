function [template,templateAPs,alignment_inds_template,tempTimes] = CalculateTemplateFromAPs(templateAPs,alignment_inds_template,template_proportion,cut_proportion,doTime,tempTimes)
    % CalculateTemplateFromAPs Description: This function standardises the
    % process of recalculating a template from variable length aligned
    % action potentials. indices that are present in less than the
    % cut_proportion of AP samples are removed completely (to prevent one
    % really long sequence making storage massive), while indices that are
    % not present in most samples (template proportion) are currently not
    % included in the template, but are stored until later. That way we
    % aren't continually cutting the sequence shorter and shorter
    % ______________
    %[template,templateAPs] = CalculateTemplateFromAPs(templateAPs,template_proportion,cut_proportion,doTime,tempTimes)
    % Inputs:
    %   - templateAPs: matrix of APs of size APlength x N_APs
    %   - template_proportion: proportion of APs that have data in order to
    %   be included in template (should be around 0.6)
    %   - alignment_inds_template: if we pad or shrink before the peak, we
    %   need to adjust the alignment index
    %   - cut_proportion: proportion of APs that have data for us to
    %   completely cut and discard those indices, should be around 0.05 (we only want to throw out when its a bit extreme) 
    %   - doTime: if we want to also keep track of the time matrix
    %   - tempTimes: times NEED TO CHECK IF I HAVE HANDLED THIS CORRECTLY
    % Outputs:
    %   - template: Template, padded with nans where we don't include
    %   indices in template
    %   - templateAPs: All the template APs, same length as template
    %   -tempTimes: the new temptimes (same size)
    % ______________
    % Author: Kathrine Clarke 
    % Date: 18-Aug-2023

     

       template  = nanmean( templateAPs, 2 );


       nAPs  = size( templateAPs, 2 );
       Nreqd = ceil( nAPs * template_proportion );
       Nreqd_cut = ceil(nAPs * cut_proportion);
       if nAPs > Nreqd
           %now that we have a sufficient number of APs we should start
           %including only indices where most have points
           Nensemble = sum( ~isnan( templateAPs ), 2 );
           ditch = Nensemble < Nreqd;
           template(ditch) = nan;

           permanent_ditch = Nensemble < Nreqd_cut;
           alignment_inds_template = alignment_inds_template-find(permanent_ditch==0,1)+1;
           template(permanent_ditch) = [];
           templateAPs(permanent_ditch,:) = [];
           if doTime
               tempTimes(permanent_ditch,:) = [];
           end
       end
 end