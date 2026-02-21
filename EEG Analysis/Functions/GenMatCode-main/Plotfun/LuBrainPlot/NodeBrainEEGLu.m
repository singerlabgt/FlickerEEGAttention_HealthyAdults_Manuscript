function NodeBrainEEGLu(ChanPos,ChanWeight,NodeTh,Color)

if isempty(ChanWeight)
   ChanWeight=ones(length(ChanPos.X),1); 
end
ChanWeight=ones(length(ChanPos.X),1); % taking this line out of the if statement above makes all nodes the same size
% x(:,1)=ChanPos.Y;
% x(:,2)=ChanPos.X;
% x(:,3)=ChanPos.Z;

x=ChanPos.ColinCoord;

xstep=(max(x)-min(x))/10;
xrange=[min(x)-xstep;max(x)+xstep];

radius=mean(nanmax(x)-nanmin(x))/10;  % make node bigger (divide by 2)

%%% Define a fixed radius for all nodes (adjust as needed) 2025-02-11 MKA
radius = 0.0025;

% Ensure color consistency
if size(Color,1)==1
   Color=repmat(Color,length(ChanWeight),1);
elseif size(Color,1)==0
   Color=repmat([0.6 0.6 0.6],length(ChanWeight),1);
else
    
end

% % if AppNode==1
% %    x(:,1)=x(:,1)-100; 
% % end
% % 
   
    
%%% Show only channel labels that exceed Node threshold NodeTh  (pre-2025-02-11)
% Index=find(abs(ChanWeight)>=NodeTh);
% for i=1:length(Index)
%     % text(x(Index(i),1)-1,x(Index(i),2),x(Index(i),3),ChanPos.labels{Index(i)}); % original. does not show up
%     text(x(Index(i),1)-0.01,x(Index(i),2),x(Index(i),3),ChanPos.labels{Index(i)}); % added 9/19/24 MKA
% end

%%% Show all channel labels (2025-02-11 MKA)
for iCh=1:length(ChanWeight)
        % Adjust vertical position for first two channels
    if iCh == 1 || iCh == 30
        yOffset = -0.002; % Slightly higher position
        xOffset = -0.005;
    elseif iCh == 14 || iCh == 18  % channel PO3 and PO4
        yOffset = -0.002; % slightly
        xOffset = -0.005;
    else
        yOffset = -0.01; % Default position
        xOffset = 0.001;
    end
    %raise all channels up
    zOffset = 0.1;
    % text(x(Index(i),1)-1,x(Index(i),2),x(Index(i),3),ChanPos.labels{Index(i)}); % original. does not show up
    text(x(iCh,1)+yOffset,x(iCh,2)+xOffset,x(iCh,3)+zOffset,ChanPos.labels{iCh}); % added 2/11/25 MKA
end

%%% Plot nodes
% Shift nodes slightly upward to appear above edges
% z_offset = 0.02; % Adjust this value based on your data scale
% x(:,3) = 0.12;
% x(:,3) = x(:,3) + z_offset;

SpherePlot(x,radius,Color,ChanWeight);