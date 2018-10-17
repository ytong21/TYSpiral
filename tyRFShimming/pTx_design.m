classdef pTx_design < handle
  % pTx_design:  a pTx pulse design class
   properties (Access = public)
 
   end
   properties (Access = protected)
     maskedMaps;  % Masked localiser, positions, b1 map and b0 map
     VesselMask;
     SliceIdx;
     gaussian_struct;
     verse_struct;
   end
   
   methods
       function obj = pTx_design(ptxFMObj,Vessel,SliceIdx)
            % constructor:  load localiser, positions, b1 map and b0 map
            maskedMaps.b1SensMasked = ptxFMObj.getB1PerV('uT','delete');
            maskedMaps.b0MapMasked = ptxFMObj.getB0('Hz','delete');
            maskedMaps.b0MapMaskedRad = (2*pi)*maskedMaps.b0MapMasked;
            maskedMaps.posVox = ptxFMObj.getPositions('cm','delete').';
            maskedMaps.localiser = ptxFMObj.getLoc('none');
            maskedMaps.b1SensMaskedHz = ptxFMObj.getB1PerV('Hz','delete');
            maskedMaps.mask = ptxFMObj.getMask();
            maskedMaps.TargetMasked = Vessel.TargetMasked;
            obj.maskedMaps = maskedMaps;
            obj.SliceIdx = SliceIdx;
       end
       
       function obj = set_vessel_mask(obj,VesselMask_in)
           obj.VesselMask = VesselMask_in;
       end
           
       function draw_vessel(obj)
           figure(96)
           LocToPlot = obj.maskedMaps.localiser(:,:,obj.SliceIdx);
           LocMax = max(LocToPlot(:));
           for iDx = 1:numel(obj.VesselMask)
              LocTemp = LocToPlot;
              LocTemp(obj.VesselMask{iDx}) = LocMax;
              imagesc(LocTemp,[min(LocToPlot(:)) max(LocToPlot(:))]); axis off;drawnow;
              pause(0.5)
           end
           close(96)
       end
       
       function AFull = getAMatrixFull(obj,RF_pulse_shape)
           % Make all gradients zero. We are only concerned about
           % on-resonance spins at the moment. Phase ramps have been
           % applied to the rf pulse to account for off-centred slices
            disp('Calculating system matrix...');
            gr = zeros(numel(RF_pulse_shape),3);
            rfOn = true(size(gr,1),1);
            tp = 10E-6 * ones(size(rfOn)); % in seconds
            DA = genAMatFull(tp,rfOn,gr,obj.maskedMaps.b1SensMaskedHz,...
                obj.maskedMaps.b0MapMasked,obj.maskedMaps.posVox);
            % To take the RF pulse shape into account. Construct a pseudo-diagonal matrix
            RFDiag = zeros(size(DA,2),8);
            RFSize = numel(bv);
            for iDx = 1:8
                RFDiag(((1:RFSize)+(iDx-1)*RFSize),iDx) = RF_pulse_shape;
            end
            AFullMag = DA*RFDiag;  %    magnetization /V
            AFull = asin(AFullMag);   % rad/V  Dupas paper in rad/V
            disp('Complete');
       end
       
       function obj = make_RF_pulses(obj,gaussian_duration_in_us,time_step_in_us,gaussian_GradAmp_in_mT_per_m)
            obj.gaussian_struct.RF_pulse = gaussian_RF_dsv(time_step_in_us, gaussian_duration_in_us, 1);
            b = [0 0 obj.RFStruct_gaussian.RF_pulse 0 0]; % arbitrary unit, peak=1
            g = [0 0 gaussian_GradAmp_in_mT_per_m*ones(size(RFStruct_Example.RF_pulse)) 0 0];% in mT/m
            dt = 10e-3;                     % in ms
            bmax = 0.5;
            gmax = 22;                      % in mT/m
            smax = 200;                     % in mT/(m*ms)
            dtout = 10e-3;                  % in ms
            emax = sum(b.^2)*dt*0.5;        % in units of b*b*dt
            [bv,gv] = mintverse(b,g,dt,bmax,gmax,smax,dtout,emax);
            bv = bv/max(bv);
            % bv, gv are column vectors, and gaussian_RF_dsv gives row
            % vectors. Need to take this into account to avoid later
            % mistakes in pulse design
            
            
            obj.gaussian_struct.GradAmp_in_mT_per_m = gaussian_GradAmp_in_mT_per_m;
            obj.verse_struct.RF_pulse = bv';
            obj.verse_struct.Gradient = gv';
            
            obj.gaussian_struct.AFull = getAMatrixFull(obj.gaussian_struct.RF_pulse);
            obj.verse_struct.AFull = getAMatrixFull(obj.verse_struct.RF_pulse);
            
       end
       

   end
 end