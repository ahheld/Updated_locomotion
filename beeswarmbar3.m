%This function plots data from excel spreadsheets. The first input must be
%the file name, the second is the data sheet, and the 3rd is an array of
%data to be read and compared. All the letters in the array must be
%strings. The 4th is how closely data should be binned. The 5th is the y
%axis limit. The 6th is the y axis numbers. The 7th should be how the y
%axis is numbered.
%
%NOTE: This script does not work well with categorical data. Try to use
%continuous data.
%
%Example: beeswarmbar3('4.1.16_Larval Movement
%Summary.xls','Combined',{'B2:B50';'L2:L50';'D2:D50';'J2:J50';'H2:H50';'F2:F50'],10,300,[0,100,200,300],{'0','100','200','300'})

function graph = beeswarmbar3(datacolumns,proximity,lim,points,pointlabels)

groupnum= numel(datacolumns);
figure('Color',[1 1 1],'Position',[100, 100,(75*groupnum)+300,500]);

%loads in data from excel spreadsheet, blocked in version 3
%for n=1:groupnum
%    data{n} = xlsread(file,sheet,strcat(datacolumns{n}));
%    data{n}(isnan(data{n}))=[]; %deletes empty positions
%end
data=datacolumns;

for n=1:groupnum
    %statistics
    avg= mean(data{n});
    counts = numel(data{n});
    SEM = std(data{n})/sqrt(counts); %mean, lower, upper CI
    
    nscale = (1.2*n)-.5; %controls how closely genotypes are grouped
    
    
    %format the data for plotting
    %scatter(nscale*ones(numel(data{n}),1),data{n}, 50, 'x','k') %plots x and y values, with x representing group number
    position = positionfinder2(data{n});
    
    %the if statement allows a specific column's data point color to be
    %changed
    if n ==4
    scatter(position,sort(data{n}),15,[.2 .2 .2],'filled');
    else
    scatter(position,sort(data{n}),15,[.2 .2 .2],'filled');
    end
    hold on
    
    %this next bit specifies bar colors and widths
    meanwidth = .5;
    SEMwidth = meanwidth/2; 
    
    
    if n == 4 %this first if statement can be used to change the color of a specific bar
    color = 'r';
    line([nscale-meanwidth, nscale+meanwidth],[avg,avg],'color',color,'linewidth',2.5); %mean
    line([nscale-SEMwidth, nscale+SEMwidth],[avg-SEM,avg-SEM],'color',color,'linewidth',2); %lower limit
    line([nscale-SEMwidth, nscale+SEMwidth],[avg+SEM,avg+SEM],'color',color,'linewidth',2); %upper limit
    line([nscale, nscale],[avg+SEM,avg-SEM],'color',color,'linewidth',2); %connects upper and lower limit
    else
    line([nscale-meanwidth, nscale+meanwidth],[avg,avg],'color','r','linewidth',2.5); %mean
    line([nscale-SEMwidth, nscale+SEMwidth],[avg-SEM,avg-SEM],'color','r','linewidth',1.5); %lower limit
    line([nscale-SEMwidth, nscale+SEMwidth],[avg+SEM,avg+SEM],'color','r','linewidth',1.5); %upper limit
    line([nscale, nscale],[avg+SEM,avg-SEM],'color','r','linewidth',1.5); %connects upper and lower limit
    end

    
end

ax.YAxis.LineWidth = 10000;
ax.XAxis.LineWidth = 10000;
xlim([0 nscale+(meanwidth*1.5)]);
%these next few lines control the y axis
ylim([0 lim])
%ylim([0 mean(data{1}*2)]);
set(gca,'Ytick',points,'YTickLabel',pointlabels);
%this next line normalizes the y label to 100% of 1st column
%set(gca,'YTickLabel',{'0','50','100','150','200'},'YTick',[0,mean(data{1})/2,mean(data{1}),mean(data{1})*1.5,mean(data{1}*2)]);%this line is specifically for locomotion
set(gca,'XTickLabel',{},'XTick',[]);
set(gca,'FontSize',20);
set(gca,'FontName','Arial');
set(gca,'Linewidth',2);

hold off


function position = positionfinder2(genotype)

%the point of this function is to scatter overlapping points in a horizontal manner
%without distorting their values (their vertical location)
    
genotype = sort(genotype);
position = [];
dataproximity = proximity; %controls horizontal spread bin size
dataspread = .1; %controls how far horizontally the data points are spread

a=1;
while a<=(numel(genotype))
    if genotype(end)==genotype(a); position = [position; nscale]; %if it's the last data point (because last point wasn't grouped)
    elseif genotype(a+1)-genotype(a)>dataproximity; %if data point is not near others, plot in middle
        position = [position; nscale];
    else group = find(genotype<(genotype(a)+dataproximity) & genotype >= genotype(a)); %if two points are within dataproximity, group them
        if mod(numel(group),2)==0; %if even
            for o=1:numel(group);
            if o==1; position = [position; nscale-(dataspread/2)];
            else position = [position; position(end)-(dataspread*(o-1)*(-1)^(o-1))]; %equation moves data left/right by dataspread
            end
            end
        else 
            for o=1:numel(group);
            if o==1; position = [position; nscale];
            else position = [position; position(end)+ (dataspread*(o-1)*(-1)^(o-1))]; %equation moves data to left/right by dataspread
            end
            end
        end
    end
    a=numel(position)+1;
end
end
end