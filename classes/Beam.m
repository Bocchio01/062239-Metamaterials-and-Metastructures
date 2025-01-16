classdef Beam

    properties
        E % Young's modulus [Pa]
        h % Height [m]
        b % Width [m]
        L % Length [m]
        rho % Density [Kg/dm^3]
        A % Cross-sectional area [m^2]
        J % Cross-sectional inertia moment [m^4]
    end
    
    methods

        function obj = Beam(varargin)

            parser = inputParser;
            addParameter(parser, 'E', 69e9);
            addParameter(parser, 'h', 1e-3);
            addParameter(parser, 'b', 20e-3);
            addParameter(parser, 'L', 24e-3);
            addParameter(parser, 'rho', 2700);
            
            % Parse the inputs
            parse(parser, varargin{:});
            
            % Assign properties
            obj.E = parser.Results.E;
            obj.h = parser.Results.h;
            obj.b = parser.Results.b;
            obj.L = parser.Results.L;
            obj.rho = parser.Results.rho;

            obj.A = obj.b .* obj.h;
            obj.J = obj.b .* (obj.h.^3) / 12;
        end

    end
end
