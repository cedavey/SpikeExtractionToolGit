function [waveform_struct] = get_spike_shapes(tseries, params)
% [waveform_struct] = get_spike_shapes(tseries, params)
% this function uses zero crossings, and estimation of noise std to find
% significant negative and positive deviations of voltage. It then extracts
% the most common combinations of positive and negative segments from the
% data, with a waveform structure to plot


dt         = double( tseries.dt );
   voltage    = double( tseries.data(:) );
   time       = double( tseries.time(:) );

posthresh  = params.positive_threshold.value;   % min pos voltage (std dev)
negthresh  = params.negative_threshold.value;   % min neg voltage (std dev)
if posthresh > 10
  str = sprintf( '\tWarning: positive threshold for spike extraction is very large (%2.2f)\n', posthresh );
  cprint( 'Keywords*', str );
end
if negthresh > 10
  str = sprintf( '\tWarning: negative threshold for spike extraction is very large (%2.2f)\n', negthresh );
  cprint( 'Keywords*', str );
end
minpostime = params.min_positive_duration.value;% min duration (ms)
minnegtime = params.min_negative_duration.value;% min duration (ms)

minpostime = round((minpostime*1e-3)/dt);       % min duration given in ms
minnegtime = round((minnegtime*1e-3)/dt);       % min duration given in ms

% Now let's find all significant postive devations
posv = voltage>0;
zerocross = find(abs(diff(posv))>0);
%first zero crossing is going between negative and positive, and last zero
%crossing is going between negative and positive - that way we get an equal
%number of negative and postive segments to work with
if posv(1)==1
    zerocross = [1;zerocross];
end
if posv(end) == 0
    zerocross = [zerocross;length(posv)];
end
dur = diff(zerocross)';

%reshape so that rows are negative and postive - can still use single
%indices for whole vector
dur = reshape(dur,2,[]);
starts = reshape(zerocross(1:end-1),2,[]);
pass_length = dur>[minpostime;minnegtime];

is_noise = ones(size(voltage));
for i = 1:length(pass_length)
    is_noise(starts(pass_length(i))+(1:dur(pass_length(i)))) = 0;
end
noise_vec = voltage(is_noise==1);
%to coursely omit glitches, lets remove most extreme 1% of noise
noise_vec(noise_vec < prctile(noise_vec,1) | noise_vec>prctile(noise_vec,99)) = [];
noise_std = std(noise_vec);

extreme_value = nan(size(dur));
extreme_value(pass_length==1) = arrayfun(@(i) max(abs(voltage(starts(i)-1+(1:dur(i))))), find(pass_length==1));
pass_height = extreme_value>[posthresh*noise_std;negthresh*noise_std];


%now we want to look at every instance where height test was passed. We
%want to find out whether length test was passed on either side, and
%if so what the height and duration of each of the values were
pass_inds = find(pass_height);
pass_inds(pass_inds==1|pass_inds==numel(dur)) = []; %omit first and last as we don't have full info
before_pass = pass_length(pass_inds-1);
after_pass = pass_length(pass_inds+1);
waveform_blocks =  sign(rem(pass_inds,2)-0.5).*[-1 1 -1]; %polarity of each block
waveform_blocks(before_pass==0,1) = 0;
waveform_blocks(after_pass==0,3) = 0;
waveform_lengths = dur(pass_inds+[-1 0 1]);
waveform_heights = extreme_value(pass_inds+[-1 0 1]);

waveform_template = struct('name','place_holder','Nsegments',0,'mags',[],'taus',[],'Nmatches',0);

waveform_struct = repmat(waveform_template,3,2,2);
first_seg = 2 - (waveform_blocks(:,1)~=0);
align_on = 1+(waveform_blocks(:,1)~=0);
for N_segs = 1:3
    for p = 1:2
        for a = 1:2
            %alignment index
           waveform_struct(N_segs,p).Nsegments = N_segs;
           %do positive first then negative
           switch N_segs
                case 1
                    waveform_struct(N_segs,p,a).name = ('monophasic');
                case 2
                    waveform_struct(N_segs,p,a).name = ('biphasic');
                case 3
                    waveform_struct(N_segs,p,a).name = ('triplet');
           end
            
           if N_segs==2 
               if p==a
                    waveform_struct(N_segs,p,a).name = [waveform_struct(N_segs,p,a).name,'_+-'];
               else
                    waveform_struct(N_segs,p,a).name = [waveform_struct(N_segs,p,a).name,'_-+'];
               end
               if p==1
                    waveform_struct(N_segs,p,a).name = [waveform_struct(N_segs,p,a).name,'_+dominant'];
                    inds = sum(waveform_blocks~=0,2)==N_segs & waveform_blocks(:,2)==1 & align_on==a;
               else
                   waveform_struct(N_segs,p,a).name = [waveform_struct(N_segs,p,a).name,'_-dominant'];
                   inds = sum(waveform_blocks~=0,2)==N_segs & waveform_blocks(:,2)==-1 & align_on==a;
               end
           else
               if p ==1
                   waveform_struct(N_segs,p,a).name = ['positive_',waveform_struct(N_segs,p,a).name];
                   inds = sum(waveform_blocks~=0,2)==N_segs & waveform_blocks(:,2)==1 & align_on==a;
               else
                   waveform_struct(N_segs,p,a).name = ['negative_',waveform_struct(N_segs,p,a).name];
                   inds = sum(waveform_blocks~=0,2)==N_segs & waveform_blocks(:,2)==-1 & align_on==a;
               end
           end
           
           
           if N_segs<3
               waveform_struct(N_segs,p,a).mags = [waveform_heights(inds==1 & first_seg==2,2+(0:N_segs-1));waveform_heights(inds==1 & first_seg==1,1+(0:N_segs-1))];
               waveform_struct(N_segs,p,a).taus = [waveform_lengths(inds==1 & first_seg==2,2+(0:N_segs-1));waveform_lengths(inds==1 & first_seg==1,1+(0:N_segs-1))];
           else
               waveform_struct(N_segs,p,a).mags = waveform_heights(inds==1,:);
               waveform_struct(N_segs,p,a).taus = waveform_lengths(inds==1,:);
           end
            waveform_struct(N_segs,p,a).Nmatches = sum(inds);

        end
    end
end



durP = dur(1:2:end);
startsP = zerocross(1:2:end-1);
pass_length
%remove those segments which are too short
startsP(durP<minpostime) = [];
durP(durP<minpostime) = [];

durN = dur(2:2:end);
startsN = zerocross(2:2:end-1);
%remove those segments which are too short
startsN(durN<minnegtime) = [];
durN(durN<minnegtime) = [];

%Now estimate noise standard deviation from all segments that are not
%candidates due to their length
[peakvalP,peakindP] = deal(zeros(size(durP)));
is_noise = ones(size(voltage));
length_status = zeros(size(voltage));
for pos_ind = 1:length(startsP) 
    is_noise(startsP(pos_ind)+(1:durP(pos_ind))) = 0;
    
    length_status(startsP(pos_ind)+(1:durP(pos_ind))) = 1;
    
    [peakvalP(pos_ind),peakindP(pos_ind)] = max(voltage(startsP(pos_ind)+(1:durP(pos_ind))));
end
[peakvalN,peakindN] = deal(zeros(size(durN)));
for neg_ind = 1:length(startsN) 
    is_noise(startsN(neg_ind)+(1:durN(neg_ind))) = 0;
    
    length_status(startsN(neg_ind)+(1:durN(neg_ind))) = -1;
    
    [peakvalN(neg_ind),peakindN(neg_ind)]  = min(voltage(startsN(neg_ind)+(1:durN(neg_ind))));
end
noise_vec = voltage(is_noise==1);
%to coursely omit glitches, lets remove most extreme 1% of noise
noise_vec(noise_vec < prctile(noise_vec,1) | noise_vec>prctile(noise_vec,99)) = [];
noise_std = std(noise_vec);


%now we can threshold remove segments that don't exceed threshold based on
%noise std

%durP


passP = peakvalP>posthresh*noise_std;
startsP = startsP(passP==1);
durP = durP(passP==1);
peakvalP = peakvalP(passP==1);
peakindP = peakindP(passP==1);

passN = peakvalN>-negthresh*noise_std;
startsN = startsN(passN==1);
durN = durN(passN==1);
peakvalP = peakvalN(passN==1);
peakindN = peakindN(passN==1);

length_and_height_status = zeros(size(voltage));
for pos_ind = 1:length(startsP) 
    length_and_height_status(startsP(pos_ind)+(1:durP(pos_ind))) = 1;
end

for neg_ind = 1:length(startsN) 
    length_and_height_status(startsN(neg_ind)+(1:durN(neg_ind))) = -1;
end
%Now that we have every individually significant positive and negative
%event, we want to look at the surrounding structure to see if there is
%some prevailing shapes

%let's align based on the positive and negative peaks

%start with positive peaks
pos_bin_image = zeros(minpostime*6,length(startsP));
pos_bin_sig_image = zeros(minpostime*6,length(startsP));
im_inds = 1:minpostime*6;
sig_ratsP = zeros(2,length(startsP));
for pos_ind = 1:length(startsP)
    ts_inds = startsP(pos_ind)+peakindP(pos_ind) + (1-minpostime*3:minpostime*3);
    keep_inds = ts_inds>0 & ts_inds<=length(length_and_height_status);
    pos_bin_image(im_inds(keep_inds==1),pos_ind) = length_status(ts_inds(keep_inds==1));
    pos_bin_sig_image(im_inds(keep_inds==1),pos_ind) = length_and_height_status(ts_inds(keep_inds==1));
    starting_part = ts_inds<startsP(pos_ind);
    ending_part = ts_inds>=startsP(pos_ind)+durP(pos_ind);
   
    
end

neg_bin_image = zeros(minpostime*6,length(startsN));
neg_bin_sig_image = zeros(minpostime*6,length(startsN));
for neg_ind = 1:length(startsN)
    ts_inds = startsN(neg_ind)+peakindN(neg_ind) + (1-minpostime*3:minpostime*3);
    keep_inds = ts_inds>0 & ts_inds<=length(length_and_height_status);
    neg_bin_image(im_inds(keep_inds==1),neg_ind) = length_status(ts_inds(keep_inds==1));
    neg_bin_sig_image(im_inds(keep_inds==1),neg_ind) = length_and_height_status(ts_inds(keep_inds==1));
end


end



