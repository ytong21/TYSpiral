function out = SimSeq(param,PulsesFromDSV)

%   Calculating dead time based on seqence params;
tEPI_per_slice = 85930e-6;
tEPI = param.num_EPI_slices*tEPI_per_slice;   %   Based on poet simulation. partial Fourier 6/8, base resolution 64
tWET = (8680*4+5000*4)*(1e-6);          %   8680 of rf event (7680 rf pusle+500ramp up/down) + 5000 gr event
tPCASL = param.tagging_duration;
tPLD = param.PLD;                       %   PLD defined as the end of PCASL and the middle of first
                                        %   excitation pulse (90 deg)
                                        
%   Adjusting tPLD. Overlaps with tEPI
tFatSat = 9600e-6+780e-6;               %   This includes the 780us gradient event after fat sat.
%tSpoiler = 2000e-6;
tRFSliceSelect = 2560e-6;               %   Origninally, PLD defined towards the middle of this pulse;
tPLD = tPLD-tRFSliceSelect/2-tFatSat;   %   Now PLD ends at the beginning of the first FAT SAT.

tDead = param.TR - tEPI - tWET - tPCASL - tPLD;
EnergyPerTR_NonPCASL = calc_EnergyPerTR_NonPCAS(param,PulsesFromDSV);

%   Starting from the WET pulses, see how far the 10s window covers.
if  param.TR < 10
    num_full_repititions = floor(10/param.TR);
    tRemaining = 10-num_full_repititions*param.TR;
elseif param.TR == 10
    num_full_repititions = 1;
    tRemaining = 0;
else
    num_full_repititions = 0;
    tRemaining = 10;
end


[sEndPoint,tModule] = determin_endpoint(tRemaining,tWET,tPCASL,tEPI,tDead);

%   Excluding PCASL and tRemaining
[EnergyLastTR_NonPCASL,tAdditionalPCASL] = calc_EnergyLastTR(sEndPoint,tModule,param,...
    EnergyPerTR_NonPCASL,PulsesFromDSV,tEPI_per_slice);

energyTot_NonPCASL_10Vref = num_full_repititions*EnergyPerTR_NonPCASL + EnergyLastTR_NonPCASL;
energyTot_NonPCASL = energyTot_NonPCASL_10Vref*(param.RefVol/10)^2;
numGaussian = num_full_repititions*round(param.tagging_duration/param.RFsep) + ...
    round(tAdditionalPCASL/param.RFsep);
out = struct('energyTot_NonPCASL',energyTot_NonPCASL,'numGaussian',...
    numGaussian,'sEndPoint',sEndPoint);
end


%   Local functions
function [sEndPoint,tModule] = determin_endpoint(tRemaining,tWET,tPCASL,tEPI,tDead)
    %   Consider WET pulses as the start, instead of the dead time.
    if tRemaining==0
        sEndPoint = 'divisible';
        tModule = 0;        
    elseif tRemaining>0 && tRemaining<=tWET
        sEndPoint = 'WET';
        tModule = tRemaining;
    elseif tRemaining>=tWET && tRemaining<tWET+tPCASL
        sEndPoint = 'PCASL';
        tModule = tRemaining - tWET;
    elseif tRemaining>=tWET+tPCASL && tRemaining<tWET+tPCASL+tEPI
        sEndPoint = 'EPI';
        tModule = tRemaining - (tWET+tPCASL);
    elseif tRemaining>=tWET+tPCASL+tEPI && tRemaining<tWET+tPCASL+tEPI+tDead;
        sEndPoint = 'Dead';
        tModule = tRemaining - (tWET+tPCASL+tEPI);
    end
end
function EnergyPerTR_NonPCASL = calc_EnergyPerTR_NonPCAS(param,PulsesFromDSV)
    EnergyPerTR_NonPCASL = PulsesFromDSV{5}.energy+...
        (PulsesFromDSV{7}.energy+PulsesFromDSV{8}.energy)*param.num_EPI_slices;
    %   WET, FatSat and Imaging pulses. 1 pair of FatSat & Imaging pulses per slice
    %   PulsesFromDSV energy based on 10V ref voltage
    %   Converting to real power based on ref vol in the main function.
end
function [EnergyLastTR_NonPCASL,tAdditionalPCASL] = calc_EnergyLastTR(sEndPoint,tModule,param,...
    EnergyPerTR_NonPCASL,PulsesFromDSV,tEPI_per_slice)
    switch sEndPoint
        case 'divisible'
            tAdditionalPCASL = 0;
            EnergyLastTR_NonPCASL = 0;
        case 'WET'
            tAdditionalPCASL = 0;
            EnergyLastTR_NonPCASL = calc_WET_energy(tModule,PulsesFromDSV);
        case 'PCASL'
            tAdditionalPCASL = tModule;
            EnergyLastTR_NonPCASL = PulsesFromDSV{5}.energy;
            % The fifth cell is all the WET pulses
        case 'EPI'
            tAdditionalPCASL = param.tagging_duration;
            EnergyLastTR_NonPCASL = calc_EPI_energy(tEPI_per_slice,tModule,PulsesFromDSV);
        case 'Dead'
            tAdditionalPCASL = param.tagging_duration;
            EnergyLastTR_NonPCASL = EnergyPerTR_NonPCASL;
    end
end
function EnergyLastTR_NonPCASL = calc_WET_energy(tModule,PulsesFromDSV)
    tFullSingleWET = (8680+5000)*1e-4;
    tRemWithinModule = rem(tModule,tFullSingleWET);
    nFullWET = (tModule-tRemWithinModule)/tFullSingleWET;
    energy_full_WETs = 0;
    if nFullWET >= 1
        for iDx = 1:nFullWET
            energy_full_WETs = energy_full_WETs+PulsesFromDSV{iDx}.energy;
            %   Note that this is still 10V ref voltage. 
        end
    end
    
    if tRemWithinModule > 0 && tRemWithinModule <= 500e-6
        energy_last_WET = 0;
    elseif tRemWithinModule > 500e-6 && tRemWithinModule <= (500+7650)*1e-6   
        % For some reason the pulses from dsv is only 7650us long
        tRemWithinModule_10usRaster = round((tRemWithinModule-500e-6)/10e-6);
        energy_last_WET = PulsesFromDSV{nFullWET+1}.pow_cumsum(tRemWithinModule_10usRaster);
    elseif tRemWithinModule > (500+7650)*1e-6
        energy_last_WET = PulsesFromDSV{nFullWET+1}.energy;
    end
    
    EnergyLastTR_NonPCASL = energy_full_WETs+energy_last_WET;
end
function EnergyLastTR_NonPCASL = calc_EPI_energy(tEPI_per_slice,tModule,PulsesFromDSV)
    tRemWithinModule = rem(tModule,tEPI_per_slice);
    nFullSlices = (tModule-tRemWithinModule)/tEPI_per_slice;
    energy_full_slices = nFullSlices*(PulsesFromDSV{7}.energy+...
        PulsesFromDSV{8}.energy);
    
    if tRemWithinModule > 0 && tRemWithinModule <= 3520e-6
        energy_last_slice = 0;
    elseif tRemWithinModule > 3520e-6 && tRemWithinModule <= (3520+2550)*1e-6
        % For some reason the pulses from dsv is only 2550us long
        tRemWithinModule_10usRaster = round((tRemWithinModule-3520e-6)/10e-6);
        energy_last_slice = PulsesFromDSV{7}.pow_cumsum(tRemWithinModule_10usRaster);
    elseif tRemWithinModule > (3520+2550)*1e-6 && tRemWithinModule <= (9600+780)*1e-6
        % FATSAT event is 9600us + additional gradients before imaging
        % pulse 780 us.
        energy_last_slice = PulsesFromDSV{7}.energy;
    elseif tRemWithinModule > (9600+780)*1e-6 && tRemWithinModule <= (9600+780+2560)*1e-6
        % Imaging pulse 2560us long
        tRemWithinModule_10usRaster = round((tRemWithinModule-(9600+780)*1e-6)/10e-6);
        energy_last_slice = PulsesFromDSV{7}.energy+...
            PulsesFromDSV{8}.pow_cumsum(tRemWithinModule_10usRaster);
    elseif tRemWithinModule > (9600+780+2560)*1e-6
        energy_last_slice = PulsesFromDSV{7}.energy + PulsesFromDSV{8}.energy;
    end    
    EnergyLastTR_NonPCASL = PulsesFromDSV{5}.energy+energy_full_slices...
        +energy_last_slice;
    
end