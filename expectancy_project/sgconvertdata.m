%% Convert 5 channels data: phase - rt, C3 - Cz - C4

FiveChannels3_phase_rt = zeros(13,4,3,3,100);
FiveChannels3_phase_rt(:,:,:,1,:) = singleEEGHILB(:,:,2:4,:);
FiveChannels3_phase_rt(:,:,:,2,:) = singleEEGHILB(:,:,2:4,:);
FiveChannels3_phase_rt(:,:,1,3,:) = avgRT(:,:,:);
FiveChannels3_phase_rt(:,:,2,3,:) = avgRT(:,:,:);
FiveChannels3_phase_rt(:,:,3,3,:) = avgRT(:,:,:);

cd('f:\raw_data\human_SG\')
save FiveChannels3_phase_rt FiveChannels3_phase_rt

%% Convert 5 channels data: phase - amp - lat, C3 - Cz - C4

FiveChannels3_phase_amp_lat = zeros(13,4,3,3,100);
FiveChannels3_phase_amp_lat(:,:,:,1,:) = singleEEGHILB(:,:,2:4,:);
FiveChannels3_phase_amp_lat(:,:,:,2,:) = singP3amphilbsort(:,:,2:4,:);
FiveChannels3_phase_amp_lat(:,:,:,3,:) = singP3lathilbsort(:,:,2:4,:);

cd('f:\raw_data\human_SG\')
save FiveChannels3_phase_amp_lat FiveChannels3_phase_amp_lat

%% Convert 5 channels data: phase - rt, Fz - Cz - Pz

FiveChannel3_phase_rt = zeros(13,4,3,3,100);
FiveChannel3_phase_rt(:,:,:,1,:) = singleEEGHILB(:,:,1:2:5,:);
FiveChannel3_phase_rt(:,:,:,2,:) = singleEEGHILB(:,:,1:2:5,:);
FiveChannel3_phase_rt(:,:,1,3,:) = avgRT(:,:,:);
FiveChannel3_phase_rt(:,:,2,3,:) = avgRT(:,:,:);
FiveChannel3_phase_rt(:,:,3,3,:) = avgRT(:,:,:);

cd('f:\raw_data\human_SG\')
save FiveChannel3_phase_rt FiveChannel3_phase_rt

%% Convert 5 channels data: phase - amp - lat, Fz - Cz - Pz

FiveChannel3_phase_amp_lat = zeros(13,4,3,3,100);
FiveChannel3_phase_amp_lat(:,:,:,1,:) = singleEEGHILB(:,:,1:2:5,:);
FiveChannel3_phase_amp_lat(:,:,:,2,:) = singP3amphilbsort(:,:,1:2:5,:);
FiveChannel3_phase_amp_lat(:,:,:,3,:) = singP3lathilbsort(:,:,1:2:5,:);

cd('f:\raw_data\human_SG\')
save FiveChannel3_phase_amp_lat FiveChannel3_phase_amp_lat