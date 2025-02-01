classdef SpatioTemporalCell

    properties
        Beam
        Piezo
    end

    properties
        E
        J
        A
        rho
    end
    
    methods
        function obj = SpatioTemporalCell(beam, piezo)
            obj.Beam = beam;
            obj.Piezo = piezo;

            obj = obj.computeAveragedProps();

        end

        function obj = computeAveragedProps(obj)

            beam = obj.Beam;
            piezo = obj.Piezo;

            obj.A = {
                beam.A + 2 * piezo.A
                beam.A
                };

            obj.J = {
                beam.J + 2 * (piezo.J + piezo.A*(beam.h/2 + piezo.h/2)^2)
                beam.J
                };
            
            obj.rho = {
                (beam.rho*beam.A + 2*piezo.rho*piezo.A) / obj.A{1}
                beam.rho
                };
            
            obj.E = {
                @(t) (beam.E*beam.J + 2*piezo.E(t)*(piezo.J + piezo.A*(beam.h/2 + piezo.h/2)^2)) / obj.J{1}
                @(t) beam.E
                };

        end
        


    end
end

