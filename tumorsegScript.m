% Script for loading Seg3D Tumor Segmentations 
% NOTE: Change line 12 for your segmentation file name convention
% NEW file convention: C#M# s000000 tumor mask.mat

addpath(genpath('/v/raid10/users/sjohnson/Matlab Code/MedImageReslicer/'))

%%%% Load the segmentation
batchdirs = {'/v/raid10/animal_data/IACUC21-04003/Cohort1_Batch1/20211027130238MR HIFU^Allison/Analysis/MRTI Recon/';
             '/v/raid10/animal_data/IACUC21-04003/Cohort1_Batch2/20211102073249MR HIFU^Allison/Analysis/MRTI Recon/';
             '/v/raid10/animal_data/IACUC21-04003/Cohort02_Batch01/20211214130817MR HIFU^Allison/Analysis/MRTI Recon/';
             '/v/raid10/animal_data/IACUC21-04003/Cohort02_Batch02/20211231085629MR HIFU^Allison/Analysis/MRTI Recon/';
             '/v/raid10/animal_data/IACUC21-04003/Cohort03_Batch01/20220127073620MR HIFU^Allison/Analysis/MRTI Recon/';
             '//v/raid10/animal_data/IACUC21-04003/Cohort03_Batch02/20220207181303MR HIFU^Allison/Analysis/MRTI Recon/'};

for b = 1:length(batchdirs)
    workdir = batchdirs{b};
    cd(workdir);
    load('mouse_file_info.mat');

    seg_path = [workdir 'tumor masks/'];


    for m = 1:length(fileinfo)
        anID = fileinfo(m).anID; 
        segtype = {'tumor_mask_BV';  'tumor_mask'};

        % load temperature data from anID for re-slicing
        tempdata = dir(['temps_' fileinfo(m).anID '*.mat']);
        load(tempdata(3).name,'temps','geomInfo','ims'); 
        targetImge = temps; 
        targetInfo = geomInfo;


        % load the geomInfo from the anID T2w dicom for re-slicing
        T2w = dir(['T2w_pre*' anID '*.mat']); 
        load(T2w(1).name,'geomInfo', 'imge_MRTI'); 
        sourceInfo = geomInfo; 
            
        for s = 1:length(segtype)
        
            % This string is the name of the exported SEG3D Segmentation:
            load([seg_path anID '_' segtype{s} '.mat'])
            tumor_seg = scirunnrrd.data; 

            %%%% Permute and flip segmentation to match MATLAB Dicom
            tumor_seg = permute(tumor_seg,[3,1,2]);
            tumor_seg = flip(tumor_seg,1); 
            tumor_seg = flip(tumor_seg,3); 

            %%%% Reslice tumor-segmentation to match MRTI data
            sourceImge = tumor_seg; 
            % reslice tumor_seg (into tumor_mask)
            [mask] = coregister_images('nearest',targetImge,targetInfo,sourceImge,sourceInfo); 
            overlayVolume(imge_MRTI, mask, 'title',['mouse ' num2str(m)]);
            uiwait(gcf); 

            eval([segtype{s} '= mask;']);
            
            % Append tumor segmentation to "T2w_pre....mat" we created earlier.
            save(T2w(1).name, segtype{s}, '-append'); 
        end 
    end 
end 