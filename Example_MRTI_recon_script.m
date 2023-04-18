% Add path where "MedImageReslicer" repo is saved in your folder
reporoot = '/v/raid10/users/sjohnson/Matlab Code/MedImageReslicer/';
addpath(genpath(reporoot)); 

% Set-up paths for saving reconstructed images 
% NOTE: You can save to your reporoot
workdir = ['/v/raid10/animal_data/IACUC21-04003/Cohort02_Batch01/20211214130817MR HIFU^Allison/Analysis/MRTI Recon/']; %some data may be saved here 
savedir_raw = [workdir]; %Where to save raw reconstructed data
savedir_dicom = [workdir]; %Where to save dicom reconstructed data

if ~exist(savedir_raw,'dir')
    mkdir(savedir_raw)
end 

if ~exist(savedir_dicom,'dir')
    mkdir(savedir_dicom)
end 

% Path where sorted Dicom folders are stored - DO NOT CHANGE
%dicomdir = '/v/raid10/animal_data/IACUC19-12014/R010-011/20201216125852MR HIFU^Henrik//DicomData/'; 
% Path were raw data is stored - DO NONT CHANGE
rawdatadir = '/v/raid10/animal_data/IACUC21-04003/Cohort02_Batch01/20211214130817MR HIFU^Allison/RawData/'; 
dicomdatadir = '/v/raid10/animal_data/IACUC21-04003/Cohort02_Batch01/20211214130817MR HIFU^Allison/DicomData/';


%% Define files & parameters for each mouse 
% T1 weighted images

%%% CHANGE FOR EACH STUDY %%%
fileinfo(1).anID = 'C4M15';
fileinfo(1).eartag = 'C4_2L';
fileinfo(1).MRTIdatafiles = dir([rawdatadir '*HIFU11c_C4_2L_*.dat']); 
fileinfo(1).T2wDicom = 's000004 t2_spc_1mm_iso_cor_TE300_2ave';
fileinfo(1).T2wDicom_save = 'T2w_pre_C4M15_s000004.mat';
fileinfo(1).BLtemp = 34.00;
fileinfo(1).ablationType = 'D';

fileinfo(2).anID = 'C6M23';
fileinfo(2).eartag = 'C6_2L';
fileinfo(2).MRTIdatafiles = dir([rawdatadir '*HIFU11c_C6_2L_*.dat']); 
fileinfo(2).T2wDicom = 's000041 t2_spc_1mm_C6_2L'; 
fileinfo(2).T2wDicom_save = 'T2w_pre_C6M23_s000041.mat';
fileinfo(2).BLtemp = 34.00;
fileinfo(2).ablationType = 'D';

fileinfo(3).anID = 'C4M13';
fileinfo(3).eartag = 'C4_1L';
fileinfo(3).MRTIdatafiles = dir([rawdatadir '*HIFU11c_C4_1L_*.dat']); 
fileinfo(3).T2wDicom = 's000078 t2_spc_1mm_C4_1L';
fileinfo(3).T2wDicom_save = 'T2w_pre_C4M13_s000078.mat';
fileinfo(3).BLtemp = 34.00;
fileinfo(3).ablationType = 'Sp';

fileinfo(4).anID = 'C5M19';
fileinfo(4).eartag = 'C5_2L';
fileinfo(4).MRTIdatafiles = dir([rawdatadir '*HIFU11c_C5_2L_*.dat']); 
fileinfo(4).T2wDicom = 's000111 t2_spc_1mm_C4_1L';
fileinfo(4).T2wDicom_save = 'T2w_pre_C5M19_s000103.mat';
fileinfo(4).BLtemp = 34.00;
fileinfo(4).ablationType = 'Sp';

save([workdir 'mouse_file_info.mat'],'fileinfo'); 
               
%% Reconstruct all Temperatures from Raw and Calculate Coordinates
%cd(savedir_raw)


load([workdir 'mouse_file_info.mat']); 

% For ZFI, set doZFI == 1. Specify desired resolution after ZFI (desR):
desR = [.5 .5 1]; % [mm] 
doZFI = 1 ; 

% Reconstruct all .dat files in loop
for m=1:4
    rawdata = fileinfo(m).MRTIdatafiles; 
for i = 1:length(rawdata)
    if i == 9
        i = i+1;
    else
        
    end
    %%%% Reconstruct data %%%%%
    filename = rawdata(i).name; 
    
    
    %%% CHANGE FOR EACH STUDY %%%
    eartg = fileinfo(m).eartag;
    
    anID = fileinfo(m).anID;
    ind = findstr(filename, eartg);
    savename = ['temps_' filename(ind:end-4)];
    
    
    %ind = findstr(filename,'_15s');
    %if isempty(ind); ind = findstr(filename,'_6p4s'); end 
    %savename = ['temps_M' num2str(m) '_' fileinfo(m).eartag filename(ind:end-4)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%% S. Johnson: fast_read_ve11_3 has 'Dim.Recon' variable, which
    %%% contains more information about how the image is supposed to be
    %%% Reconed from the scanner. For example, we acquired at 1x1 in-plane,
    %%% but reconned at 0.5x0.5 in-plane 
    [ks2, fileInfo] =  fast_read_ve11_3([rawdatadir filename],'NoFillToFullFourier','ShiftDCToMatrixCenter', 'NoGUI','ReadFileInfo','LastMeasDatOnly','CollapseSeg'); 
    ks2 = ks2{end}; % still has RO oversampling 
    header = fileInfo.Protocol(end); 
    ks2.dim = lower(ks2.dim); %Changes ks2.dim from capitalized to lower case strings
    [ks2.data, ks2.dim, ~] = cellpermute(ks2.data, ks2.dim,{'col' 'lin' 'par' 'rep' 'cha' 'acq' 'ave' 'eco'});
    
    %%% Determine if Slice Partial Fourier is on, perform circshift of k-space if necessary %%%
    switch header.sKSpace.ucSlicePartialFourier
        case 4
            sPFfraction = 3/4;
        case 8
            sPFfraction = 6/8; 
            % Circshift a slice around 
            sliceshift = (size(ks2.data,3) - size(ks2.data,3)/sPFfraction)/2; 
            ks2.data = circshift(ks2.data,sliceshift,3); 
        case 16
            sPFfraction = 1; 
            % Circshift a slice around 
            sliceshift = (size(ks2.data,3) - size(ks2.data,3)/sPFfraction)/2; 
            ks2.data = circshift(ks2.data,sliceshift,3); 
    end 
    
    %%% Determine data dimensionality %%%
    if header.sKSpace.lPartitions == 1
        DIM = 2; 
        ks2.data = permute(ks2.data,[1:2 ndims(ks2.data)+1 3:ndims(ks2.data)]); 
        ks2.dim = [ks2.dim(1:2), 'slc', ks2.dim(3:end)]; % Add the 'slc' in to the dim vector
        ks = ks2.data; 
    else 
        DIM = 3; 
        ks = ks2.data; 
    end 
     
    %%%% Convert to Image space %%%%%
    [ims, invcov] = CreateImg3(ks, DIM, 3,[1 1], []);   close(gcf); 
    ims = rmOs(ims); %remove RO oversampling from the image
    % take care of slice PF zerofilling:
    if sPFfraction ~= 1
        ims = ZFI_Filter2(ims,1,1,1/sPFfraction,DIM,'N',1,1,'n',[]); 
    end 
     
    % Get acquisition resolution  
    AcqResolution = getImageResolution(header);
        
    %%%% ZFI %%%%
     if doZFI == 1
        KMatrixSize = [size(ims,1) size(ims,2), size(ims,3)];
        FOV = KMatrixSize.*AcqResolution;

        zfiFactor = FOV./(desR.*KMatrixSize); % calculate ZFI factor based on desR
        imsZFI = ZFI_Filter2(ims,zfiFactor(1),zfiFactor(2),zfiFactor(3),DIM,'N',1,1,'n',[]); 
        
        % OPTIONAL: reconned data will be saved to file with ZFI factor in
        % the name:
        zfistr = sprintf('%.2fx%.2fx%.2f',zfiFactor(1), zfiFactor(2), zfiFactor(3)); 
        savename = [savename '_ZFI_' zfistr '.mat']; 
        
        %%%% Calculate MR Coordinates %%%%
        [imsZFI, PosPCS, PosDCS, geomInfo] = calc_coords_RAW(header,imsZFI,zfiFactor);
        ims = imsZFI;
    else 
        zfiFactor = [1 1 1];
        savename= [savename '.mat']; 
        
        %%%% Calculate MR Coordinates %%%%
        [ims, PosPCS, PosDCS, geomInfo] = calc_coords_RAW(header, ims);
    end 
        
    %%% Calculate temperatures %%%%
    TE = header.alTE/1000000;
    B0 = header.sProtConsistencyInfo.flNominalB0;
    BLtemp = fileinfo(m).BLtemp;
    temps = CalcPhDiff(ims, TE(1),[1,4],B0) + fileinfo(m).BLtemp; 
     
    
    %%%% Save variables %%%%%
    save([workdir  savename],'ks','header','ims','temps','geomInfo','BLtemp');
    
end 
end 

%% Load an MRTI Dicom

dicomname = 's000026 segEPI_HIFU11c_C3M9_heat9_56pct_15s'; 
load_image_DICOM(dicomdatadir, dicomname, savedir_dicom,0,1); 
cd(workdir)   
%% Load Dicoms & make masks

dicomdata.fnames = {'s000004 t2_spc_1mm_iso_cor_TE300_2ave'; 
                 's000041 t2_spc_1mm_C6_2L';
                 's000078 t2_spc_1mm_C4_1L';
                 's000111 t2_spc_1mm_C4_1L'};

load([workdir 'mouse_file_info.mat']); 
savedir_dicom = workdir; 

for m = 1:4
 
    dicomname = fileinfo(m).T2wDicom;
    load_image_DICOM(dicomdatadir, dicomname, savedir_dicom,0,1); 
    
    cd(workdir)
    tempdata = dir(['temps_' fileinfo(m).anID '*.mat']);
    load(tempdata(2).name,'temps','geomInfo'); 
    targetImge = temps; 
    targetInfo = geomInfo; 
    
    load([dicomname '.mat'],'imge','geomInfo'); 
    sourceImge = imge; 
    sourceInfo = geomInfo; 
    
    [imge_MRTI] = coregister_images('nearest',targetImge,targetInfo,sourceImge,sourceInfo); 
    overlayVolume(imge_MRTI,temps(:,:,:,7),'title',['mouse ' num2str(m)],'overRange',[33 55],'overMap','jet',...
                               'maskRange',[2 Inf]); 
    
    save([dicomname '.mat'],'imge_MRTI','-append'); 
    
    %Rename files: 
    movefile([dicomname  '.mat'], [fileinfo(m).T2wDicom_save]); 
    
    %%% OPTIONAL: Visualize bounding boxes %%%%%%
    %Visualize bounding boxes 
    PosT2w = meshgrid_from_affine(targetInfo.AffineDCS, imge); 
    PosMRTI = meshgrid_from_affine(sourceInfo.AffineDCS, temps);
    figure; 
    plotPos(PosT2w,[],1,'b'); hold on; 
    plotPos(PosMRTI,[],1,'r'); 
    title({'DICOM: blue';'MRTI: red'})
end 

%% SPECIAL CASE: Correct shift in C4M13 MRTI data
% S. Johnson 3/8/2022
% MRTI data is shiftedd din the HF direction relative to other scans in
% mouse 14, for an unknown reason. This section applies a manual
% translation to the MRTI geomInfo affine matrices such that the agar hole
% in MRTI and T2w images align after reslicing. 

load([workdir 'mouse_file_info.mat']); 
m = 3; 
% load ims of MRTI, get magnitude image. Generate DCS meshgrid
tempdata = dir(['temps_' fileinfo(m).anID '*.mat']);
load(tempdata(3).name,'ims','temps','geomInfo');
geom_temps = geomInfo; 
DCS_temps = meshgrid_from_affine(geom_temps.AffineDCS, squeeze(ims(:,:,:,1)));
mag = abs(ims); 
% load T2w image 
load(fileinfo(m).T2wDicom_save, 'imge','geomInfo');
geom_t2 = geomInfo; 
DCS_t2 = meshgrid_from_affine(geom_t2.AffineDCS, imge); 

% use MRDataFig to locate the top of the agar hole (in HF direction)
%MRDataFig(1, 'mag'); 
%MRDataFig(2, 'imge'); 

% calculate the HF offset between MRTI and T2w image in DCS coords
MRloc_temps = squeeze(DCS_temps(138, 63, 4,:))';
MRloc_t2 = squeeze(DCS_t2(53, 110, 32,:))';
offset = MRloc_t2(3) - MRloc_temps(3) % HF offset in mm 


% translate the DCS and PCS matrices (same translation in HFP patient
% position)
geom_temps2 = geom_temps; 
T = zeros(size(geom_temps.AffineDCS));
T(3,4) = offset; 
geom_temps2.AffineDCS = geom_temps.AffineDCS + T; 
geom_temps2.AffinePCS = geom_temps.AffinePCS + T; 


% re-do the T2w reslice: 
targetImge = temps; 
targetInfo = geom_temps2; 

sourceImge = imge; 
sourceInfo = geom_t2; 

[imge_MRTI] = coregister_images('nearest',targetImge,targetInfo,sourceImge,sourceInfo); 
overlayVolume(imge_MRTI,temps(:,:,:,7),'title',['mouse ' num2str(m)],'overRange',[33 55],'overMap','jet',...
                           'maskRange',[2 Inf]); 

save([fileinfo(m).T2wDicom_save],'imge_MRTI','-append'); 

% save the new geomInfo variable to all "temps" in M13 
geomInfo = geom_temps2; 
for f = 1:length(tempdata)
    save(tempdata(f).name, 'geomInfo','-append'); 
end 

%% 










%%%% END RECONSTRUCTION %%%%%