% Plotting
  if ~exist('plotRolArray','var')
      figure(56)
      imagesc(maskedMaps.localiser(:,:,SliceIdx));
      title('Select Rect Regions to Plot')
      h2plot = imrect;  wait(h2plot);
      plotMask = logical(createMask(h2plot));
      [rol,cow] = find(plotMask);
      plotRolArray = rol(1):rol(end);
      plotCowArray = cow(1):cow(end);     
  end
 % Other plots
 %%
figure(90)
set(gcf,'color','w','InvertHardcopy','off')
set(gcf,'units','centimeters','position',[4 4 40 20],'paperunits','centimeters','paperposition',[0 0 40 20])
clf
PlotFALim = [0 20]; FontSize = 9;   CLB = cell(7,1);   CLBFontSize = 8;
Figs = cell(4,3);

Figs{2,1} = subplot(3,2,3);
imagesc(ImgToPlot.b1CP(plotRolArray,plotCowArray,SliceIdx)'); title('DREAM B1 Map','FontSize',FontSize);axis off;
CLB{1} = colorbar('FontSize',CLBFontSize);
CLB{1}.YLabel.String = 'Hz/V';    CLB{1}.YLabel.Rotation = 0;     
CLB{1}.YLabel.Position = [0.5 4.2 0];
hold on;colormap(Figs{2,1},'gray')
visboundaries(maskedMaps.mask(plotRolArray,plotCowArray,SliceIdx)','Color','r','LineWidth',0.4,...
    'EnhanceVisibility', true,'LineStyle','-');axis off;
nudge(Figs{2,1},[0 0.05 0 0])

Figs{3,1} = subplot(3,2,5); 
imagesc(ImgToPlot.b0(plotRolArray,plotCowArray,SliceIdx)',[-200 200]); title('B0 Map','FontSize',FontSize);axis off;hold on;
CLB{2} = colorbar('FontSize',CLBFontSize);
CLB{2}.YLabel.String = 'Hz';    CLB{2}.YLabel.Rotation = 0;     
%CLB{2}.YLabel.Position = [0 235 0];
CLB{2}.YLabel.Position = [0.5 235 0];
colormap(Figs{3,1},'gray')
visboundaries(maskedMaps.mask(plotRolArray,plotCowArray,SliceIdx)','Color','r','LineWidth',0.4,...
    'EnhanceVisibility', true,'LineStyle','-');axis off;
nudge(Figs{3,1},[0 0.1 0 0])

Figs{1,1} = subplot(3,2,1); 
imagesc(maskedMaps.localiser(plotRolArray,plotCowArray,SliceIdx)'); title('TOF image for vessel masking','FontSize',FontSize);axis off;colormap(Figs{1,1},'gray')
hold on
visboundaries(maskedMaps.mask(plotRolArray,plotCowArray,SliceIdx)','Color','r','LineWidth',0.4,...
    'EnhanceVisibility', true,'LineStyle','-');axis off;
hp3 = get(subplot(3,2,1),'Position');
set(Figs{1,1},'Position',[hp3(1) hp3(2) Figs{3,1}.Position(3) hp3(4)]);



HoriOffset = 0.038;
Figs{1,2} = subplot(3,4,3);
imagesc(FullImage.CP(plotRolArray,plotCowArray,2)',PlotFALim); axis off;
title('Flip angle maps','FontSize',FontSize);
ylabel('CP Mode','FontSize',FontSize,'FontWeight','bold');
Figs{1,2}.YLabel.Visible = 'on';
set(Figs{1,2},'Position',[Figs{1,2}.Position(1) Figs{1,2}.Position(2) Figs{1,2}.Position(3) Figs{1,2}.Position(4)]);
colormap(Figs{1,2},'hot')

Figs{2,2} = subplot(3,4,7);
imagesc(FullImage.PhaseOnly(plotRolArray,plotCowArray,2)',PlotFALim);axis off;
ylabel('Phase Only','FontSize',FontSize,'FontWeight','bold');
Figs{2,2}.YLabel.Visible = 'on';
set(Figs{2,2},'Position',[Figs{2,2}.Position(1) Figs{2,2}.Position(2) Figs{2,2}.Position(3) Figs{2,2}.Position(4)]);
colormap(Figs{2,2},'hot')
nudge(Figs{2,2},[0 0.05 0 0]);

Figs{3,2} = subplot(3,4,11);
imagesc(FullImage.AS(plotRolArray,plotCowArray,2)',PlotFALim);axis off;
set(Figs{3,2},'Position',[Figs{3,2}.Position(1) Figs{3,2}.Position(2) Figs{3,2}.Position(3) Figs{3,2}.Position(4)]);
ylabel('Full RF Shimming','FontSize',FontSize,'FontWeight','bold');
Figs{3,2}.YLabel.Visible = 'on';
CLB{6} = colorbar('southoutside','FontSize',CLBFontSize);
CLB{6}.YLabel.String = 'Degrees (°)';    CLB{6}.YLabel.Rotation = 0;   
colormap(Figs{3,2},'hot')
nudge(Figs{3,2},[0 0.1 0 0])
nudge(CLB{6},[0 0.02 0 0])

PlotDiffLim = [-20 20];
Figs{1,3} = subplot(3,4,4);
imagesc(DiffImage.CP(plotRolArray,plotCowArray,2)',PlotDiffLim); 
title('Difference images','FontSize',FontSize);axis off;
set(Figs{1,3},'Position',[Figs{1,3}.Position(1)-HoriOffset Figs{1,3}.Position(2) ...
    Figs{1,3}.Position(3) Figs{1,3}.Position(4)]);
colormap(Figs{1,3},'parula')

Figs{2,3} = subplot(3,4,8);
imagesc(DiffImage.PhaseOnly(plotRolArray,plotCowArray,2)',PlotDiffLim); axis off;
set(Figs{2,3},'Position',[Figs{2,3}.Position(1)-HoriOffset Figs{2,3}.Position(2) ...
    Figs{2,3}.Position(3) Figs{2,3}.Position(4)]);
colormap(Figs{2,3},'parula')
nudge(Figs{2,3},[0 0.05 0 0])


Figs{3,3} = subplot(3,4,12);
imagesc(DiffImage.AS(plotRolArray,plotCowArray,2)',PlotDiffLim);axis off
set(Figs{3,3},'Position',[Figs{3,3}.Position(1)-HoriOffset Figs{3,3}.Position(2) ...
    Figs{3,3}.Position(3) Figs{3,3}.Position(4)]);
CLB{7} = colorbar('southoutside');
CLB{7} = colorbar('southoutside','FontSize',CLBFontSize);
colormap(Figs{3,3},'parula')
nudge(Figs{3,3},[0 0.1 0 0])
nudge(CLB{7},[0 0.02 0 0])

CLB{7}.YLabel.String = 'Percentage error (%)';    CLB{7}.YLabel.Rotation = 0;     
%CLB{7}.YLabel.Position = [0.5 700 0];

hp4 = get(Figs{3,3},'Position');
%colorbar('Position', [hp4(1)+hp4(3)+0.015  hp4(2)  0.02  0.815],'FontSize',FontSize);



%%
figure(92)
Histo = cell(3,1);
[f1,xi] = ksdensity(abs(FAFinal.AS),linspace(15,25,500));
plot(xi,f1,'b');hold on
[f1,xi] = ksdensity(abs(FAFinal.CP),linspace(15,25,500));
plot(xi,f1,'r');
[f1,xi] = ksdensity(abs(FAFinal.PhaseOnly),linspace(15,25,500));
plot(xi,f1,'m');
%%
figure(93)
subplot(1,2,1)
set(gcf,'color','w','InvertHardcopy','off')
set(gcf,'units','centimeters','position',[4 4 40 20],'paperunits','centimeters','paperposition',[0 0 40 20])
%histogram(abs(FAFinal.AS),20,'Normalization','probability');hold on;
%[f1,xi] = ksdensity(abs(FAFinal.AS),linspace(17,22,100));
%plot(xi,f1,'b');

% Histo{2,1} = histfit(abs(FAFinal.CP),10,'kernel'); delete(Histo{2,1}(1));  Histo{2,1}(2).Color = 'c';  hold on
% Histo{3,1} = histfit(abs(FAFinal.PhaseOnly),10,'kernel');  delete(Histo{3,1}(1));  Histo{3,1}(2).Color = 'b';
% Histo{1,1} = histfit(abs(FAFinal.AS),10,'kernel');  delete(Histo{1,1}(1));  Histo{1,1}(2).Color = 'r'; 

edges = 16:0.1:24;

Histo{1,1} = histogram(abs(FAFinal.CP),edges,'FaceColor','b');    Histo{1,1}.Normalization = 'probability'; hold on
Histo{2,1} = histogram(abs(FAFinal.PhaseOnly),edges,'FaceColor','g');    Histo{2,1}.Normalization = 'probability';
Histo{3,1} = histogram(abs(FAFinal.AS),edges,'FaceColor','y');    Histo{3,1}.Normalization = 'probability';
box off
xlabel('Flip angle (°)')
ylabel('Probability'); ylim([0 0.18])
lgd = legend('CP mode','Phase only shimming','Full RF shimming');legend('boxoff')   
lgd.FontSize = 13;set(gca, 'FontSize', 16)

  subplot(1,2,2)
  set(gcf,'color','w','InvertHardcopy','off')
  set(gcf,'units','centimeters','position',[4 4 40 20],'paperunits','centimeters','paperposition',[0 0 40 20])
  Mode = {'CP Mode','Phase only','Full RF shimming'};   VesselName = {'RICA','RVA','LICA','LVA','Total'};
  BarChart = bar(1:5,BarMtx);  lgdTmp = legend(Mode); box off;  ylim([0.6 1]);  ylabel('Labelling efficiency');  
  lgdTmp.FontSize = 13; set(gca,'XTickLabel',VesselName);set(gca, 'FontSize', 16);legend('boxoff')  
  %pause(0.1);
  hold on
  for iDx = 1:numel(BarChart)
    xData = BarChart(iDx).XData+BarChart(iDx).XOffset;
    errorbar(xData,BarMtx(:,iDx),StdMtx(:,iDx),'r.')
  end
  BarChart(1).FaceColor = 'b'; BarChart(2).FaceColor = 'g'; BarChart(3).FaceColor = 'y';
  hold off
  
   %%
 matched = dt.searchForSeriesInstanceNumber(3,80);
 imgObj = Spectro.dicomImage({matched.Filename});
 
  matched = dt.searchForSeriesInstanceNumber(21,1);
 ASLObj = Spectro.dicomImage({matched.Filename});
 %%
 figure(123)
set(gcf,'color','w','InvertHardcopy','off')
set(gcf,'units','centimeters','position',[4 4 40 20],'paperunits','centimeters','paperposition',[0 0 40 20])
subplot(1,2,1)
 clf
 himpatch = Spectro.PlotCsi.impatch2d(imgObj.image,imgObj.info{1}.PixelSpacing);
 h = gca;
 h.YDir = 'reverse';
 
 % Draw sat band on the localiser image

% Read RSat geometry from the protocol

rsatNum = ASLObj.getMrProtocolNumber('sRSatArray.lSize');

for rsatDx = 1:rsatNum
    rsat(rsatDx).position = [...
        ASLObj.getMrProtocolNumber(sprintf('sRSatArray.asElm[%d].sPosition.dSag',rsatDx-1));...
        ASLObj.getMrProtocolNumber(sprintf('sRSatArray.asElm[%d].sPosition.dCor',rsatDx-1));...
        ASLObj.getMrProtocolNumber(sprintf('sRSatArray.asElm[%d].sPosition.dTra',rsatDx-1))];
    
    rsat(rsatDx).normal = [...
        ASLObj.getMrProtocolNumber(sprintf('sRSatArray.asElm[%d].sNormal.dSag',rsatDx-1));...
        ASLObj.getMrProtocolNumber(sprintf('sRSatArray.asElm[%d].sNormal.dCor',rsatDx-1));...
        ASLObj.getMrProtocolNumber(sprintf('sRSatArray.asElm[%d].sNormal.dTra',rsatDx-1))];
    
    rsat(rsatDx).thickness = ASLObj.getMrProtocolNumber(sprintf('sRSatArray.asElm[%d].dThickness',rsatDx-1));
end

for slices = 1:20
    rsat(2+slices).position = ASLObj.imagePositionPatient;
    rsat(2+slices).position(3) = rsat(2+slices).position(3) +(slices-1)*ASLObj.info{1}.SliceThickness;
    rsat(2+slices).normal =ASLObj.sliceNormal;
    rsat(2+slices).thickness = ASLObj.info{1}.SliceThickness;
    rsatNum = rsatNum+1;
end
% Mark intersection of saturation band with displayed reference images
    axes(h);
    
    
    % Now calculate and plot new lines
    infoRef = imgObj.info{1};

    % Ref is the plane of the image on screen. Rectangle.
    ref_geom.imageOrientationPatient = ...
        reshape(infoRef.ImageOrientationPatient,3,2);
    ref_geom.normal = ...
        cross(ref_geom.imageOrientationPatient(:,1),...
              ref_geom.imageOrientationPatient(:,2));
    ref_geom.unitVecs = ...
        [ref_geom.imageOrientationPatient ...
        ref_geom.normal];
    
    % Store dimensions
    ref_geom.pixelSpacing = infoRef.PixelSpacing;
    ref_geom.rows = double(infoRef.Rows);
    ref_geom.columns = double(infoRef.Columns);

    % Mark corner positions
    ref_geom.imagePositionPatient = ...
        infoRef.ImagePositionPatient;
    
    ref_geom.imagePositionEndRow1 = ref_geom.imagePositionPatient...
    + ref_geom.imageOrientationPatient(:,1) * ref_geom.pixelSpacing(2) * ref_geom.columns;

    ref_geom.imagePositionEndColumn1 = ref_geom.imagePositionPatient...
    + ref_geom.imageOrientationPatient(:,2) * ref_geom.pixelSpacing(1) * ref_geom.rows;
  
    ref_geom.imagePositionFarCorner = ...
        ref_geom.imagePositionPatient + ...
        ref_geom.imageOrientationPatient(:,1) * ref_geom.pixelSpacing(2) * ref_geom.columns + ...
        ref_geom.imageOrientationPatient(:,2) * ref_geom.pixelSpacing(1) * ref_geom.rows;
    
    for rsatDx=1:rsatNum
        %% Calculate the intersection of two planes.
        % Plane 2 is the RSat. Infinite extent.
        rsat_geom = rsat(rsatDx);
        
        % Test for parallel planes
        if norm(cross(ref_geom.normal,rsat_geom.normal)) < 1e-14
            warning('Spectro:PlotCsi:calcPlaneIntersect','Parallel sat band. Cannot yet plot intersection.')
            
            hLine = NaN;
            pointOnLine = [];
            directionVector = [];
            
            continue
        end
        
        % unitVecs' is a ROTATION MATRIX from i, j, k PATIENT COORDS to the
        % coordinate system of the reference image
        % (and unitVecs goes the other way).
        
        ref_geom.RotMatrix = ref_geom.unitVecs';
        
        minusPRef = dot(ref_geom.normal, ref_geom.imagePositionPatient);
        minusP2 = dot(rsat_geom.normal, rsat_geom.position);

        minusP2 = minusP2 - rsat_geom.thickness / 2; % Show edges of sat band not centre
        
        % Calculate direction vector along line of intersection
        directionVector = null([ref_geom.normal rsat_geom.normal]');
        directionVectorProjected = ref_geom.RotMatrix*directionVector;
        
        % Calculate a point on the (extension of the) line of intersection
        pointOnLine = [ref_geom.normal rsat_geom.normal]' \ [minusPRef; minusP2];
        pointOnLineProjected = ref_geom.RotMatrix*(pointOnLine - infoRef.ImagePositionPatient);
        
        % % Check commonPoint is in both planes
        % maxdiff(dot(normalRef,pointOnLine),minusPRef)
        % maxdiff(dot(normal2,pointOnLine),minusP2)
        
        %% Find the start and end of this line, by seeing whether/where it
        %% intersects the four sides of the displayed reference image rectangle.
        
        % We will do this calculation in the projected plane to make it a 2D rather
        % than full 3D geometric problem.
        
        % TODO: I'm sure these don't need calculating like this... They are just
        % pixelspace .* [rows columns] or something like that.
        cornerOrigin_InPlane = ref_geom.unitVecs'*(ref_geom.imagePositionPatient - ref_geom.imagePositionPatient);
        cornerEndRow1_InPlane = ref_geom.unitVecs'*(ref_geom.imagePositionEndRow1 - ref_geom.imagePositionPatient);
        cornerEndColumn1_InPlane = ref_geom.unitVecs'*(ref_geom.imagePositionEndColumn1 - ref_geom.imagePositionPatient);
        cornerFar_InPlane = ref_geom.unitVecs'*(ref_geom.imagePositionFarCorner - ref_geom.imagePositionPatient);
        
        res{1} = Spectro.PlotCsi.intersectLineAndLineSegment(pointOnLineProjected,directionVectorProjected,cornerOrigin_InPlane,cornerEndRow1_InPlane);
        res{2} = Spectro.PlotCsi.intersectLineAndLineSegment(pointOnLineProjected,directionVectorProjected,cornerOrigin_InPlane,cornerEndColumn1_InPlane);
        res{3} = Spectro.PlotCsi.intersectLineAndLineSegment(pointOnLineProjected,directionVectorProjected,cornerEndRow1_InPlane,cornerFar_InPlane);
        res{4} = Spectro.PlotCsi.intersectLineAndLineSegment(pointOnLineProjected,directionVectorProjected,cornerEndColumn1_InPlane,cornerFar_InPlane);
        
        lineEnds = zeros(3,0);
        for idx=1:4
            lineEnds = [lineEnds res{idx}.points]; %#ok<AGROW>
        end
        
        %% Now calculate for the other side of the sat band
        minusP2_plus = dot(rsat_geom.normal, rsat_geom.position) + rsat_geom.thickness / 2;
        
        % Calculate a point on the (extension of the) line of intersection
        pointOnLine_plus = [ref_geom.normal rsat_geom.normal]' \ [minusPRef; minusP2_plus];
        pointOnLineProjected_plus = ref_geom.RotMatrix*(pointOnLine_plus - infoRef.ImagePositionPatient);

        res_plus{1} = Spectro.PlotCsi.intersectLineAndLineSegment(pointOnLineProjected_plus,directionVectorProjected,cornerOrigin_InPlane,cornerEndRow1_InPlane);
        res_plus{2} = Spectro.PlotCsi.intersectLineAndLineSegment(pointOnLineProjected_plus,directionVectorProjected,cornerOrigin_InPlane,cornerEndColumn1_InPlane);
        res_plus{3} = Spectro.PlotCsi.intersectLineAndLineSegment(pointOnLineProjected_plus,directionVectorProjected,cornerEndRow1_InPlane,cornerFar_InPlane);
        res_plus{4} = Spectro.PlotCsi.intersectLineAndLineSegment(pointOnLineProjected_plus,directionVectorProjected,cornerEndColumn1_InPlane,cornerFar_InPlane);
        
        lineEnds_plus = zeros(3,0);
        for idx=1:4
            lineEnds_plus = [lineEnds_plus res_plus{idx}.points]; %#ok<AGROW>
        end
        
        %% Check whether any vertices are included
        inRegion = @(x) sign(dot(x - pointOnLineProjected, cross(directionVectorProjected,[0 0 1]))) ~= sign(dot(x - pointOnLineProjected_plus, cross(directionVectorProjected,[0 0 1])));
        
        allPts = [lineEnds lineEnds_plus];
        if inRegion(cornerOrigin_InPlane)
            allPts(:,end+1) = cornerOrigin_InPlane;
        end
        if inRegion(cornerEndRow1_InPlane)
            allPts(:,end+1) = cornerEndRow1_InPlane;
        end
        if inRegion(cornerEndColumn1_InPlane)
            allPts(:,end+1) = cornerEndColumn1_InPlane;
        end
        if inRegion(cornerFar_InPlane)
            allPts(:,end+1) = cornerFar_InPlane;
        end
        
        allPts = allPts(1:2,:);
        allPtsDx = convhull(allPts.');
        
        allPts = allPts(:,allPtsDx);
        
        %% Special handling for Matlab R2014b to prevent clicking on the patch
        if verLessThan('matlab','8.4')
            hitTestProp = {'HitTest','off'};
        else
            hitTestProp = {'HitTest','off','PickableParts','none'};
        end
        
        %% Plot all
        hold on
        if rsatDx==2 
        hPatch = patch(allPts(1,:),allPts(2,:),[1 1 0],'FaceAlpha',0.2,'Tag','RSat',hitTestProp{:},'UserData',rsatDx);
        hLine = line(lineEnds(1,:),lineEnds(2,:));
        set(hLine,'xliminclude','off','yliminclude','off','Color',[1 1 0],'Tag','RSat',hitTestProp{:},'UserData',rsatDx)
        
        hLine_plus = line(lineEnds_plus(1,:),lineEnds_plus(2,:));
        set(hLine_plus,'xliminclude','off','yliminclude','off','Color',[1 1 0],'Tag','RSat',hitTestProp{:},'UserData',rsatDx)
        
        elseif rsatDx>2 
            hPatch = patch(allPts(1,2:3),allPts(2,2:3),[1 0 0],'edgecolor','r','FaceAlpha',0.2,'Tag','RSat',hitTestProp{:},'UserData',rsatDx);
        end
        
    end

 TOF = subplot(1,2,2);
 imagesc(maskedMaps.localiser(:,:,12)');colormap(TOF,'gray');axis off
% hold on
% visboundaries(maskedMaps.mask(:,:,12)','Color','r','LineWidth',0.4,...
%     'EnhanceVisibility', true);axis off;
% Anno = cell(4,1);
% Anno{1,1} = annotation('textarrow',[38 175]/10,[92 199]/10,'String','RICA');
%%
Spectro.PlotCsi.impatch2d(ASLObj.image,ASLObj.info{1}.PixelSpacing)
