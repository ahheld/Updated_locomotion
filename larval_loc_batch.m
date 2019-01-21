function [out] = larval_loc_batch(genotype,sec,iteration)

files = glob(['Data\',genotype,'*']); %find genotypes named in GUI in data folder
distances = [];
for n = 1:numel(files)
    try %this try loop should prevent rogue videos from disrupting the whole information flow
    distance = larval_loc_func_2(files{n},sec);
    num = n*ones(numel(distance),1);
    distances = [distances; distance num];
    end
end

genotypename = genotype(1:end-6);

xlswrite('locomotion_out',{genotypename},'Summary',[char(iteration*2+64),'1']);
xlswrite('locomotion_out',{'videonumber'},'Summary',[char(iteration*2+65),'1']);
xlswrite('locomotion_out',distances,'Summary',[char(iteration*2+64),'2']);

out = distances(:,1);