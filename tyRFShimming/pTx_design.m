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
     param;
   end
   
   methods
       function obj = pTx_design(ptxFMObj,Vessel,SliceIdx,param)
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
            obj.param = param;
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
            RFSize = numel(RF_pulse_shape);
            for iDx = 1:8
                RFDiag(((1:RFSize)+(iDx-1)*RFSize),iDx) = RF_pulse_shape;
            end
            AFullMag = DA*RFDiag;  %    magnetization /V
            AFull = asin(AFullMag);   % rad/V  Dupas paper in rad/V
            disp('Complete');
       end
       
       function obj = make_RF_pulses(obj,gaussian_duration_in_us,time_step_in_us,gaussian_GradAmp_in_mT_per_m)
            RF_temp = gaussian_RF_dsv(time_step_in_us, gaussian_duration_in_us, 1);
            obj.gaussian_struct.RF_pulse = RF_temp(1,:);
            b = [0 0 obj.gaussian_struct.RF_pulse 0 0]; % arbitrary unit, peak=1
            g = [0 0 gaussian_GradAmp_in_mT_per_m*ones(size(obj.gaussian_struct.RF_pulse)) 0 0];% in mT/m
            dt = time_step_in_us*1e-3;                     % in ms
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
            
            obj.gaussian_struct.AFull = obj.getAMatrixFull(obj.gaussian_struct.RF_pulse);
            obj.verse_struct.AFull = obj.getAMatrixFull(obj.verse_struct.RF_pulse);
       end
       function [bOut,out_struct] = runVE(obj,AFullFunc,CGtikhonov)
            MaxVEIter = 50;
            cost = inf;
            xCurr = zeros(obj.param.numCh,1);
            Tol = 0.001;    
            %targetFAInRad = maskedMaps.TargetMasked*deg2rad(param.targetFlipAngle);
            targetFAInRad = ones(size(obj.maskedMaps.b0MapMasked))*deg2rad(obj.param.targetFlipAngle);

            phiTarget = angle(sum(obj.maskedMaps.b1SensMaskedHz,2));
            z = exp(1i*phiTarget);
            fullImage =  zeros(numel(obj.maskedMaps.mask),1);

            ALambda = pinv(AFullFunc'*AFullFunc + CGtikhonov*eye(numel(xCurr)))*AFullFunc';
            for iDx = 1:MaxVEIter

                targetFA = z .* targetFAInRad;
            %     [xs,~ ] = qpwls_pcg(xCurr, AFullFunc,1, targetM, 0, CGtikhonov, 1, cgIterations);
            %     xCurr = xs(:,end);
                xCurr = ALambda*targetFA;currFA = AFullFunc*xCurr;
                costNew = norm(currFA - targetFA) + CGtikhonov*(xCurr'*xCurr);
                costDiff = abs(cost-costNew)/cost;
                if costDiff < Tol
                    break
                end
                cost = costNew;

                % Smooth phase
                fullImage(logical(obj.maskedMaps.mask(:))) = angle(currFA);
                fullImage = reshape(fullImage,size(obj.maskedMaps.mask));
                fullImage = imgaussfilt3(fullImage,2,'FilterSize',5);

                phiTarget = fullImage(obj.maskedMaps.mask(:));
                z = exp(1i*phiTarget);    
            end

            bOut = xCurr;
            out_struct.finalPwr = norm(xCurr);
            %finalCost = costNew;
            %finalCost = cost;
            finalMag = zeros(numel(obj.maskedMaps.mask),1);
            finalMag(logical(obj.maskedMaps.mask(:))) = AFullFunc*xCurr;
            out_struct.finalMag = reshape(finalMag,size(obj.maskedMaps.mask));
            %finalRMSE = norm(currFA - targetFA);
            finalRMSE = norm(abs(currFA) - abs(targetFA))/norm(abs(targetFA));
            disp(finalRMSE)
            out_struct.finalRMSE = finalRMSE;
       end
       function [bOut,out_struct] = runAS(obj,bVE,RF_pulse,AFull)
                       %Setting upper and lower bounds
            maxV = 239; % Check where I got this.
            clear ub lb
            ub(1:obj.param.numCh) = maxV;
            ub((1+obj.param.numCh):(obj.param.numCh*2)) = inf;

            lb(1:obj.param.numCh) = 0;
            lb((1+obj.param.numCh):(obj.param.numCh*2)) = -inf;

            %Setting optimization options

            optionsFMin = optimoptions(@fmincon);
            optionsFMin.Algorithm = 'active-set';
            optionsFMin.Display = 'none';
            optionsFMin.MaxFunctionEvaluations = obj.param.MaxEvaluation;
            optionsFMin.SpecifyConstraintGradient = false;
            optionsFMin.OptimalityTolerance = obj.param.tol;
            optionsFMin.FiniteDifferenceType = 'central';
            optionsFMin.MaxIterations = 40000;
            protectedModeConstraints = CoilConstraints.novaCoil( true );
            nonlincon = @(x) TYpowerConstraints_AS_Global_OneSpoke(x,protectedModeConstraints,...
                obj.param.RFsep,RF_pulse);

            % Creating a function handle based on getAMatSimp
            xInitial = [abs(bVE);angle(bVE)];
            %TargetFA = maskedMaps.TargetMasked*deg2rad(param.targetFlipAngle);
            TargetFA = ones(size(obj.maskedMaps.b0MapMasked))*deg2rad(obj.param.targetFlipAngle);
            %FunHandle = @(x) norm(abs(AFull*(x(1:8).*exp(1i*x(9:16)))) - TargetFA);
            FunHandle = @(x) (norm(abs(AFull*(x(1:8).*exp(1i*x(9:16)))) - TargetFA)/norm(TargetFA));
            %FunHandle = @(x) sum(abs(abs(AFull*(x(1:8).*exp(1i*x(9:16))))-TargetFA))/sum(TargetFA);
            %FunHandle = @(x) -goodnessOfFit(abs(AFull*(x(1:8).*exp(1i*x(9:16)))),TargetFA,'NRMSE');
            %options = optimoptions('fmincon','MaxIterations',4000);
            [bOut,out_struct.fval,out_struct.exitflag,out_struct.output] = fmincon(FunHandle,xInitial,...
                [],[],[],[],lb.',ub.',nonlincon,optionsFMin);
       end
       function PulseStruct = runShimCalculation(obj,PulseStruct,tikhonovArray)
              PulseStruct.bVE = zeros(8,numel(tikhonovArray));
              PulseStruct.VE_output = cell(1,numel(tikhonovArray));
              for iDx = 1:numel(tikhonovArray)
                  [PulseStruct.bVE(:,iDx),PulseStruct.VE_output{iDx}] = ...
                       obj.runVE(PulseStruct.AFull,tikhonovArray(iDx));
                  disp(PulseStruct.VE_output{iDx}.finalRMSE)  
              end
              PulseStruct.bAS = zeros(16,numel(tikhonovArray));
              PulseStruct.AS_output = cell(1,numel(tikhonovArray));
              for iDx = 1:size(PulseStruct.bVE,2)
                    [PulseStruct.bAS(:,iDx),PulseStruct.AS_output{iDx}] = runAS(obj,PulseStruct.bVE(:,iDx),...
                        PulseStruct.RF_pulse,PulseStruct.AFull);
              end
              error_vec = zeros(1,numel(PulseStruct.AS_output));
              for iDx = 1:numel(PulseStruct.AS_output)
                  error_vec(iDx) = PulseStruct.AS_output{iDx}.fval;
              end
              [~, minIndex] = min(error_vec);
              PulseStruct.minIndex = minIndex;
              PulseStruct.bmin = PulseStruct.bAS(1:8,minIndex).*exp(1i*PulseStruct.bAS(9:16,minIndex));
       end
       function obj = run_gaussian(obj)
           obj.gaussian_struct = obj.runShimCalculation(obj.gaussian_struct,power(10,-7:0.5:-3));
       end
       function obj = run_verse(obj)
           obj.verse_struct = obj.runShimCalculation(obj.verse_struct,power(10,-7:0.5:-3));
       end
       function magnetization = runBlochSim(obj,PulseStruct)
            RFToSim = bsxfun(@times,repmat(PulseStruct.RF_pulse,1,8),b.');
            GToSim = zeros(numel(RFToSim)/8,3);
            %GToSim(:,3) = gv;
            %GToSim(:,3) = 425.8*RFStruct.G_amp*ones(numel(RFToSim)/8,1);     %from mT/m to Hz/cm
            B1ToSim = obj.maskedMaps.b1SensMaskedHz;
            magnetization = make_blochSim(RFToSim,B1ToSim,...
                obj.maskedMaps.b0MapMasked,GToSim,10E-6,obj.maskedMaps.posVox,obj.maskedMaps.mask);
            SlicesMask = obj.maskedMaps.mask(:,:,obj.SliceIdx);
            PulseStruct.FullImage = zeros(size(SlicesMask));
            PulseStruct.FullImage = asind(magnetization.mxy);
            
       end
   end
   methods (Access = public)
       function out = get_gaussian(obj)
           out = obj.gaussian_struct;
       end
       function out = get_verse(obj)
           out = obj.verse_struct;
       end       
   end
 end