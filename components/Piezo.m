classdef Piezo

    properties
        Shunt
    end

    properties
        Y11_E % Piezoelectric Young modulus in short circuit [Pa]
        Y11_D % Piezoelectric Young modulus in open circuit [Pa]
        C_T   % Capacitance at constant stress [F]
        k31   % Electromechanical coupling coefficient
        h     % Thickness [m]
        b     % Width [m]
        L     % Length [m]
        rho   % Density [Kg/m^3]
        C_S   % Stiffness-adjusted capacitance [F]
        E
        A
        J
    end
    
    methods
        function obj = Piezo(varargin)
            
            % Input parser
            parser = inputParser;
            addParameter(parser, 'Y11_E', 62e9);
            addParameter(parser, 'C_T', 7.0e-9);
            addParameter(parser, 'k31', 0.351);
            addParameter(parser, 'h', 1e-3);
            addParameter(parser, 'b', 20e-3);
            addParameter(parser, 'L', 22e-3);
            addParameter(parser, 'rho', 7900);
            
            % Parse inputs
            parse(parser, varargin{:});
            
            % Assign properties
            obj.Y11_E = parser.Results.Y11_E;
            obj.C_T = parser.Results.C_T;
            obj.k31 = parser.Results.k31;
            obj.h = parser.Results.h;
            obj.b = parser.Results.b;
            obj.L = parser.Results.L;
            obj.rho = parser.Results.rho;
            
            % Derived properties
            obj.C_S = obj.C_T * (1 - obj.k31^2);
            obj.Y11_D = obj.Y11_E / (1 - obj.k31^2); % Series connection
            obj.E = @(t) obj.Y11_E;
            obj.A = obj.b * obj.h;
            obj.J = obj.b * (obj.h^3) / 12;
        
        end

        function obj = bindShunt(obj, shunt, w)

            obj.Shunt = shunt;
            % obj.E = @(t) obj.Y11_E * (1 + obj.k31^2 * 1i*w*obj.C_T*obj.Shunt.Z(w) / (1 + 1i*w*obj.C_S*obj.Shunt.Z(w)));
            obj.E = @(t) obj.Y11_D * (1 - obj.k31^2 / (1 + 1i*w*obj.C_S*obj.Shunt.Z(w)));

        end
        
    end
end
