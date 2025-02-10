function [STC, modulation] = assemble_STC(modulation_label, modulation_frequency, components)

arguments
    modulation_label {mustBeMember(modulation_label,['OFF-OFF-OFF','ON-ON-ON','ON-OFF-OFF','Sinusoidal (continuos)','Sinusoidal (discrete)'])} = 'Sinusoidal (discrete)'
    modulation_frequency {mustBeNumeric} = 0
    components.beam  Beam  = Beam()
    components.piezo Piezo = Piezo()
    components.shunt Shunt = Shunt('C-')
end

beams(1:3)  = deal(components.beam);
piezos(1:3) = deal(components.piezo);
E_piezo_OFF = components.piezo.bindShunt(Shunt('OFF'), 2*pi * modulation_frequency).E(0);
E_piezo_ON  = components.piezo.bindShunt(components.shunt, 1e6).E(0);

% Modulation data
modulation.label      = modulation_label;
modulation.lambda     = sum([beams.L]);
modulation.wavenumber = 2*pi / modulation.lambda;
modulation.omega      = 2*pi * modulation_frequency;
modulation.period     = min(abs(1 / modulation_frequency), 1);
% modulation.period     = 1 / modulation_frequency;
modulation.mean       = (E_piezo_OFF + E_piezo_ON) / 2;
modulation.amplitude  = (E_piezo_OFF - E_piezo_ON) / 2;

switch modulation.label

    case 'OFF-OFF-OFF'
        piezos(1) = piezos(1).bindShunt(Shunt('OFF'), 1);
        piezos(2) = piezos(2).bindShunt(Shunt('OFF'), 1);
        piezos(3) = piezos(3).bindShunt(Shunt('OFF'), 1);

    case 'ON-ON-ON'
        piezos(1) = piezos(1).bindShunt(components.shunt, 1e6);
        piezos(2) = piezos(2).bindShunt(components.shunt, 1e6);
        piezos(3) = piezos(3).bindShunt(components.shunt, 1e6);

    case 'ON-OFF-OFF'
        piezos(1) = piezos(1).bindShunt(components.shunt, 1e6);
        piezos(2) = piezos(2).bindShunt(Shunt('OFF'), 1);
        piezos(3) = piezos(3).bindShunt(Shunt('OFF'), 1);
        warning('You should use Sinusoidal (discrete) @0Hz instead')

    case 'Sinusoidal (continuos)'
        piezos(1).E = @(t) modulation.mean + modulation.amplitude * cos(modulation.omega*t + (1-1)*2*pi/3);
        piezos(2).E = @(t) modulation.mean + modulation.amplitude * cos(modulation.omega*t + (2-1)*2*pi/3);
        piezos(3).E = @(t) modulation.mean + modulation.amplitude * cos(modulation.omega*t + (3-1)*2*pi/3);

    case 'Sinusoidal (discrete)'
        piezos(1).E = @(t) modulation.mean + modulation.amplitude * sign( cos(modulation.omega*t + (1-1)*2*pi/3) );
        piezos(2).E = @(t) modulation.mean + modulation.amplitude * sign( cos(modulation.omega*t + (2-1)*2*pi/3) );
        piezos(3).E = @(t) modulation.mean + modulation.amplitude * sign( cos(modulation.omega*t + (3-1)*2*pi/3) );

    otherwise
        error('Unknown modulation model: %s', modulation.label)

end

STC = [
    SpatioTemporalCell(beams(1), piezos(1))
    SpatioTemporalCell(beams(2), piezos(2))
    SpatioTemporalCell(beams(3), piezos(3))
    ];

end

