function [] = linePlot(varargin)

if isa(varargin{1},'matlab.graphics.axis.Axes')
    ha = varargin{1};
    varargin = varargin(2:end);
else
    ha = gca;
end

if isa(varargin{1},'cell')
    dataCell = varargin{1};
else
    dataCell = {varargin{1}};
end

plotArgs = varargin(2:end);

for i_ = 1:numel(dataCell)
    
    data = dataCell{i_};

szM1 = size(data);
plotArgs = varargin(2:end);

if szM1(1)==2
    xdata = squeeze(data(1,:,:));
    ydata = squeeze(data(2,:,:));

    plot(ha,xdata,ydata,plotArgs{:});
elseif szM1(1)==3
    xdata = squeeze(data(1,:,:));
    ydata = squeeze(data(2,:,:));
    zdata = squeeze(data(3,:,:));

    plot3(ha,xdata,ydata,zdata,plotArgs{:});
end

hold(ha,'on');
end
hold(ha,'off');

end

