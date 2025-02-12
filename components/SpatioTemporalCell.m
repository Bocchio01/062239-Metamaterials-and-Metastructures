classdef SpatioTemporalCell
    % SPATIOTEMPORALCELL Class representing a spatio-temporal unit cell.
    % This class models a composite structure composed of a beam and
    % piezoelectric elements, and computes averaged mechanical properties.

    properties
        Beam  % Instance of Beam class representing the beam component
        Piezo % Instance of Piezo class representing the piezoelectric component
    end

    properties
        E   % Effective Young's modulus function handle
        J   % Effective moment of inertia [m^4]
        A   % Effective cross-sectional area [m^2]
        rho % Effective mass density [kg/m^3]
    end

    methods
        function obj = SpatioTemporalCell(beam, piezo)
            % SPATIOTEMPORALCELL Constructor for the SpatioTemporalCell class.
            % Initializes the unit cell with a beam and piezoelectric material.
            %
            % Syntax:
            %   obj = SpatioTemporalCell(beam, piezo)
            %
            % Inputs:
            %   beam  - Instance of Beam class
            %   piezo - Instance of Piezo class
            %
            % Output:
            %   obj - Instance of SpatioTemporalCell with computed properties.

            obj.Beam = beam;  % Assign beam object
            obj.Piezo = piezo; % Assign piezoelectric object

            % Compute averaged properties of the unit cell
            obj = obj.computeAveragedProps();
        end

        function obj = computeAveragedProps(obj)
            % COMPUTEAVERAGEDPROPS Computes the effective mechanical properties.
            % This function calculates the averaged area, inertia, density,
            % and Young's modulus of the spatio-temporal unit cell.
            %
            % Output:
            %   obj - Updated instance of SpatioTemporalCell with computed properties.

            beam = obj.Beam;
            piezo = obj.Piezo;

            % Compute effective cross-sectional areas
            obj.A = {
                beam.A + 2 * piezo.A; % Total area with two piezo layers
                beam.A;               % Beam-only area
                };

            % Compute effective moments of inertia
            obj.J = {
                beam.J + 2 * (piezo.J + piezo.A * (beam.h/2 + piezo.h/2)^2); % Composite inertia
                beam.J; % Beam-only inertia
                };

            % Compute effective mass densities
            obj.rho = {
                (beam.rho * beam.A + 2 * piezo.rho * piezo.A) / obj.A{1}; % Composite density
                beam.rho; % Beam-only density
                };

            % Compute effective Young's modulus as a function handle
            obj.E = {
                @(t) real((beam.E * beam.J + 2 * piezo.E(t) * (piezo.J + piezo.A * (beam.h/2 + piezo.h/2)^2)) / obj.J{1}); % Composite modulus
                @(t) beam.E; % Beam-only modulus
                };
        end
    end
end
