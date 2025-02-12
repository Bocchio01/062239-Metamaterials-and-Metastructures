classdef Beam
    % BEAM Class representing a beam with mechanical properties.
    % This class defines a beam with properties such as Young's modulus,
    % cross-sectional dimensions, density, and inertia moment.

    properties
        E   % Young's modulus [Pa] (Material stiffness)
        h   % Beam height [m] (Thickness of the beam section)
        b   % Beam width [m]
        L   % Beam length [m]
        rho % Density [kg/m^3] (Material mass density)
        A   % Cross-sectional area [m^2] (Calculated property)
        J   % Cross-sectional moment of inertia [m^4] (Calculated property)
    end

    methods
        function obj = Beam(varargin)
            % BEAM Constructor for the Beam class.
            % Initializes beam properties and computes area and inertia.
            %
            % Syntax:
            %   obj = Beam('E', value, 'h', value, ...)
            %
            % Name-Value Pairs:
            %   'E'   - Young's modulus [Pa] (default: 69e9)
            %   'h'   - Beam height [m] (default: 1e-3)
            %   'b'   - Beam width [m] (default: 20e-3)
            %   'L'   - Beam length [m] (default: 24e-3)
            %   'rho' - Density [kg/m^3] (default: 2700)
            %
            % Output:
            %   obj - Instance of the Beam class.

            % Create input parser
            parser = inputParser;
            addParameter(parser, 'E', 69e9);    % Default Young's modulus (69 GPa)
            addParameter(parser, 'h', 1e-3);    % Default height (1 mm)
            addParameter(parser, 'b', 20e-3);   % Default width (20 mm)
            addParameter(parser, 'L', 24e-3);   % Default length (24 mm)
            addParameter(parser, 'rho', 2700);  % Default density (2700 kg/m^3)

            % Parse the inputs
            parse(parser, varargin{:});

            % Assign parsed values to properties
            obj.E = parser.Results.E;
            obj.h = parser.Results.h;
            obj.b = parser.Results.b;
            obj.L = parser.Results.L;
            obj.rho = parser.Results.rho;

            % Compute additional beam properties
            obj.A = obj.b * obj.h;          % Cross-sectional area
            obj.J = obj.b * (obj.h^3) / 12; % Cross-sectional moment of inertia
        end
    end
end
