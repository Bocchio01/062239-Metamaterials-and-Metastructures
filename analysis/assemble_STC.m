function [STC, modulation] = assemble_STC(STC_modulation_label, Z_model_label)

shunts(1:3) = deal(Shunt(Z_model_label));
beams(1:3)  = deal(Beam());
piezos(1:3) = deal(Piezo());

% Modulation data
modulation.label      = STC_modulation_label;
modulation.lambda     = sum([beams.L]);
modulation.wavenumber = 2*pi / modulation.lambda;

switch modulation.label
    
    case {'OFF-OFF-OFF', 'ON-ON-ON', 'ON-ON-OFF'}
        modulation.omega     = 0;
        modulation.period    = 1;
        modulation.amplitude = 0;

    case {'Sinusoidal (continuos)', 'Sinusoidal (discrete)'}
        modulation.omega  = 0.2 * sqrt(beams(1).E / beams(1).rho) * modulation.wavenumber;
        modulation.period = 2*pi / modulation.omega;

        E_sandwich_ON  = abs(SpatioTemporalCell(Beam(), Piezo().bindShunt(Shunt(Z_model_label), modulation.omega)).E{1}(0));
        E_sandwich_OFF = abs(SpatioTemporalCell(Beam(), Piezo()).E{1}(0));
   
        modulation.amplitude = (E_sandwich_ON - E_sandwich_OFF) / (E_sandwich_ON + E_sandwich_OFF);
        modulation.amplitude = 27.5 / 100;

        clear E_sandwich_OFF E_sandwich_ON

    otherwise
        error('Unknown modulation model: %s', modulation.label)

end


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
        piezos(1).E = @(t) (beams(1).E * (1 + modulation.amplitude * cos(2*pi*1/modulation.period*t + (1-1)*2*pi/3 )));
        piezos(2).E = @(t) (beams(2).E * (1 + modulation.amplitude * cos(2*pi*1/modulation.period*t + (2-1)*2*pi/3 )));
        piezos(3).E = @(t) (beams(3).E * (1 + modulation.amplitude * cos(2*pi*1/modulation.period*t + (3-1)*2*pi/3 )));

    case 'Sinusoidal (discrete)'
        piezos(1).E = @(t) (beams(1).E * (1 + modulation.amplitude * sign( cos(2*pi*1/modulation.period*t + (1-1)*2*pi/3 ))));
        piezos(2).E = @(t) (beams(2).E * (1 + modulation.amplitude * sign( cos(2*pi*1/modulation.period*t + (2-1)*2*pi/3 ))));
        piezos(3).E = @(t) (beams(3).E * (1 + modulation.amplitude * sign( cos(2*pi*1/modulation.period*t + (3-1)*2*pi/3 ))));

    otherwise
        error('Unknown modulation model: %s', modulation.label)

end

STC = [
    SpatioTemporalCell(beams(1), piezos(1))
    SpatioTemporalCell(beams(2), piezos(2))
    SpatioTemporalCell(beams(3), piezos(3))
    ];

end

