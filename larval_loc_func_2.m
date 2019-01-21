%% This function analyzes distance traveled of adult flies over a selected number of frames. The function is meant to work with wmv files taken with the dinolite camera in the Wharton lab. Supply two inputs: the video name and number of frames to be analyzed.

% example: To analyze frames 1:1200 in LoxP_2.wmv, enter:
% adult_loc_func('LoxP_2.wmv',1200)
%
% The output will be a vector of distance traveled for each fly and a figure showing distance traveled in pixels 

function [real_dist_traveled] = larval_loc_func_2(video_in,sec)

sec = str2num(sec); %bug where sec wasn't being read as a number

%% this section lays the groundwork for opening frames for analysis
vid = VideoReader(video_in);
frame_count = 1:10:(sec*30+1); %this will be the limitation on frames read
allcenters = cell(numel(frame_count),1); %cell array storing object location information
bowlcenter = [320,240];


%disp(frame_count(end))

% find the right threshold for videos
positives = [];

for n=.2:.1:.6
    frame = read(vid,1);
    BW = im2bw(frame,n);
    bound = bwboundaries(BW,'noholes');
    %next part finds objects of correct size
    [row,column] = cellfun(@size,bound);
    bound_filtered = bound(row>5 & row<30); %between 10 and 25 seems optimal
    delete = zeros(numel(bound_filtered),1);
    for nn=1:length(delete) %this loop deletes objects outside the bowl using bowlcenter
        testcoord = bound_filtered{nn}(1,:);
        if ((bowlcenter(1)-testcoord(2)).^2+(bowlcenter(2)-testcoord(1)).^2).^.5>150
            delete(nn)=1;
        end
    end
    for nn=1:length(delete)
        if delete(nn)==1
            bound_filtered{nn}=[];
        end
    end
    bound_filtered = bound_filtered(~cellfun(@isempty,bound_filtered));
    positives = [positives; n,numel(bound_filtered)];
    
    %this part shows boundaries
    %imshow(BW);
    %hold on
    %for k = 1:length(bound_filtered);
    %    boundary = bound_filtered{k};
    %    plot(boundary(:,2),boundary(:,1),'g');
    %end
    %scatter(320,240,'b') %center
end

%need to delete thresholds with >5 objects detected
save = positives(:,2)<6;
positives = [positives(save,1),positives(save,2)];


[val,pos] = max(positives(:,2));
threshold = pos*.1+.1;

%% this section opens frames, applies a filter, and identifies points that are flies

for n = 1:numel(frame_count) %reads up to and including frame number in frame_count ----------- possibly use parfor instead? or use parfor when calling this function?
    frame = read(vid,frame_count(n));
    BW = im2bw(frame,threshold); %threshold calculated above
    bound = bwboundaries(BW,'noholes');
    %next part finds objects of correct size
    [row,column] = cellfun(@size,bound);
    bound_filtered = bound(row>5 & row<30); %between 10 and 25 seems optimal   
    centers = cell2mat(cellfun(@mean,bound_filtered,'UniformOutput',false)); %make centers into matrix
    allcenters{n}=centers; %put matrix into cell array
    
    %this part shows boundaries
    %imshow(BW);
    %hold on
    %for k = 1:length(bound_filtered);
    %    boundary = bound_filtered{k};
    %    plot(boundary(:,2),boundary(:,1),'g');
    %    scatter(centers(k,2),centers(k,1),'r');
    %end
end

vid = []; %the video object takes up a lot of CPU usage, so I delete it afterward

%% next section will be to connect the points between frames
% In order to identify the correct fly during collisions, I need to connect
% the assembled fly with the next frame's points in order of most confident
% interaction (lowest distance between points across frames). I made an interaction matrix between fly endpoints and potential next points, and assembly is done by ordered confidence. 

objnum = length(allcenters{1}(:,1)); %number of objects in first frame
%need to build in elimination of background objects

start = [];
for n=1:objnum
    score = ((320-allcenters{1}(n,2)).^2+(240-allcenters{1}(n,1)).^2).^.5; %find objects' distances from center
    if score<150 %if object is less than 150 from center at start, then it counts as a real object
        start = [start; allcenters{1}(n,:)];
    end
end

objnum = numel(start)/2;        

assembly = cell(objnum,1); %this is the cell matrix for aligning the positions between frames


for n = 1:objnum %this line assigns the first 
    assembly{n}=start(n,:);
end

for n=1:(length(allcenters)-1) %go through all frames
    interactions = []; %interaction matrix
    match_objects = allcenters{n}; %points to be matched
    for nn=1:objnum %for flies detected
        fly = assembly{nn}(end,:);
        for nnn=1:numel(match_objects(:,1)) %for match objects
            dist = ((fly(1)-match_objects(nnn,1)).^2+(fly(2)-match_objects(nnn,2)).^2).^.5; %calculate distance between objects
            interactions = [interactions; nn, nnn, dist, match_objects(nnn,1),match_objects(nnn,2)]; %make vector of fly#, match_object#, distance between, match_object_x, match_object_y
        end
    end
    %interaction matrix complete, start matching to existing assembly
    while numel(interactions)>0
        [val,loc] = min(interactions(:,3)); %find row # of shortest interaction
        match = interactions(loc,:); %values for shortest interaction
        object = assembly{match(1)}; %existing assembly for fly number (value is match(1)) 
        object = [object; match(4:5)]; %tag on next point (value is match(4:5))
        if val<50; assembly{match(1)} = object; end %place matched points into assembly if value of minimum interaction>30
        delete_fly = interactions(:,1)==match(1); %identify matched fly in interaction matrix
        interactions(delete_fly,:) = []; %delete fly from interaction matrix (do not want to match twice in one frame)
        delete_match_obj = interactions(:,2)==match(2);
        interactions(delete_match_obj,:) = [];
    end
end

    
%% now that points have been strung together, need to distinguish background noise from real movement
% looked at fly that did not move, and picked up max variation of 1.1 pixels in
% distance between points

trace =cell(objnum,1); %this will be the trace of the distance traveled by the fly

for n = 1:objnum
    trace{n}=assembly{n}(1,:); %load trace array with first point
end


objects_dist_travel = []; %this is the matrix for storing distance traveled
for n=1:objnum %this is number of detected flies
    for nn=1:length(assembly{n}(:,1)) %loops for comparing test points for distance traveled >1.5
        objtrace = trace{n}; %this is the individual object being tracked
        start = objtrace(end,:); %start point for comparisons
        test = assembly{n}(nn,:); %points to be tested
        if (((start(1,1)-test(1,1)).^2+(start(1,2)-test(1,2)).^2).^.5)>1.5 %1.5 is conservative movement above noise
            objtrace = [objtrace; test]; %tag successful test points onto trace
            trace{n}=objtrace; %insert objtrace into trace
        end
    end
    
    %this for loop calculates distance traveled for trace
    sum_distance = 0;
    if length(trace{n}(:,1))>1 %if more than one point in trace
        for nn=1:(length(trace{n})-1)
            start = trace{n}(nn,:); %start point
            next = trace{n}((nn+1),:); %next point
            distance = ((start(1,1)-next(1,1)).^2+(start(1,2)-next(1,2)).^2).^.5; %distance between points
            sum_distance = sum_distance + distance; %total distance traveled
        end
    end
    objects_dist_travel = [objects_dist_travel; sum_distance]; %save distance traveled
end

real_dist_traveled = objects_dist_travel*.3351; %adjustment for frame scale, output in millimeters

%% make plot of distances traveled for visual verification
figure %('Visible','off'); %make figure, but keep it invisible
colorloop = ['r','g','b','w','y','m','c'];%loops colors used for fly visualization
colorloop  = [colorloop colorloop colorloop colorloop colorloop]; %makes colors for up to 32 flies
%scatters distance traveled above noise
imshow(BW)
hold on
for n=1:objnum
    scatter(trace{n}(:,2),trace{n}(:,1),colorloop(n),'filled');
    tx= text((trace{n}(1,2)+5),(trace{n}(1,1)+5),strcat('larva',int2str(n),':',int2str(real_dist_traveled(n)))); %text is shifted 5 up and 5 right from first point, text displays movement
    tx.Color = colorloop(n);
end

tx2 = text(1,10,strcat(int2str(sec),'sec analyzed'));
tx2.Color = 'w';

%set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') %this line of code changes the property of the figure to make it visible when opened
saveas(gcf,[video_in(1:(end-4)),'_analysis.fig'],'fig') %this saves the figure in its current form
close %closes the figure