function [STC, modulation] = assemble_STC(modulation_label, modulation_frequency, components)

arguments
    modulation_label {mustBeMember(modulation_label,['OFF-OFF-OFF','ON-ON-ON','ON-ON-OFF','Sinusoidal (continuos)','Sinusoidal (discrete)'])} = 'Sinusoidal (discrete)'
    modulation_frequency {mustBeNumeric} = 0
    components.beam  Beam  = Beam()
    components.piezo Piezo = Piezo()
    components.shunt Shunt = Shunt('C-')
end

beams(1:3)  = deal(components.beam);
piezos(1:3) = deal(components.piezo);
shunts(1:3) = deal(components.shunt);

% Modulation data
modulation.label      = modulation_label;
modulation.lambda     = sum([beams.L]);
modulation.wavenumber = 2*pi / modulation.lambda;
modulation.omega      = 2*pi * modulation_frequency;
% modulation.period     = min(1 / modulation_frequency, 1);
modulation.period     = 1 / modulation_frequency;

E_sandwich_ON  = abs(SpatioTemporalCell(components.beam, components.piezo.bindShunt(components.shunt, modulation.omega)).E{1}(0));
E_sandwich_OFF = abs(SpatioTemporalCell(components.beam, components.piezo).E{1}(0));

modulation.mean = (E_sandwich_ON + E_sandwich_OFF) / 2;
modulation.amplitude = - (E_sandwich_ON - E_sandwich_OFF) / (E_sandwich_ON + E_sandwich_OFF);

% modulation.mean = components.beam.E;
% modulation.amplitude = 0.6;

switch modulation.label

    case 'OFF-OFF-OFF'
        % piezos(1) = piezos(1).bindShunt(shunts(1), modulation.omega);
        % piezos(2) = piezos(2).bindShunt(shunts(2), modulation.omega);
        % piezos(3) = piezos(3).bindShunt(shunts(3), modulation.omega);

    case 'ON-ON-ON'
        piezos(1) = piezos(1).bindShunt(shunts(1), modulation.omega);
        piezos(2) = piezos(2).bindShunt(shunts(2), modulation.omega);
        piezos(3) = piezos(3).bindShunt(shunts(3), modulation.omega);

    case 'ON-ON-OFF'
        piezos(1) = piezos(1).bindShunt(shunts(1), modulation.omega);
        piezos(2) = piezos(2).bindShunt(shunts(2), modulation.omega);
        % piezos(3) = piezos(3).bindShunt(shunts(3), modulation.omega);

    case 'Sinusoidal (continuos)'
        piezos(1).E = @(t) (modulation.mean * (1 + modulation.amplitude * cos(modulation.omega*t + (1-1)*2*pi/3 )));
        piezos(2).E = @(t) (modulation.mean * (1 + modulation.amplitude * cos(modulation.omega*t + (2-1)*2*pi/3 )));
        piezos(3).E = @(t) (modulation.mean * (1 + modulation.amplitude * cos(modulation.omega*t + (3-1)*2*pi/3 )));

    case 'Sinusoidal (discrete)'        
        piezos(1).E = @(t) (modulation.mean * (1 + modulation.amplitude * sign( cos(modulation.omega*t + (1-1)*2*pi/3 ))));
        piezos(2).E = @(t) (modulation.mean * (1 + modulation.amplitude * sign( cos(modulation.omega*t + (2-1)*2*pi/3 ))));
        piezos(3).E = @(t) (modulation.mean * (1 + modulation.amplitude * sign( cos(modulation.omega*t + (3-1)*2*pi/3 ))));

    otherwise
        error('Unknown modulation model: %s', modulation.label)

end

STC = [
    SpatioTemporalCell(beams(1), piezos(1))
    SpatioTemporalCell(beams(2), piezos(2))
    SpatioTemporalCell(beams(3), piezos(3))
    ];

end

