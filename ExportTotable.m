clear 
workdir = '/v/raid10/users/sjohnson/Experiment Analysis/IACUC21-04003/FinalSubjectData/'
addpath(genpath('/v/raid10/users/sjohnson/Matlab Code/Packages/')); 

cd(workdir)
filelist = dir('C*B*.mat'); 

for f = 1:length(filelist)
    load(filelist(f).name)
    
    
    for m = 1:length(mousedata)
        if f == 1 && m == 1
            allmousedata = mousedata(m); 
        else
            allmousedata = cat(1, allmousedata, mousedata(m));  
        end 
    end 
end 


mrtiTable = struct2tableAll(allmousedata);

writetable(mrtiTable, 'MRTIstats.csv'); 