classdef Piezo
    % PIEZO Class representing a piezoelectric material with electromechanical properties.
    % This class models a piezoelectric material with its physical and
    % electrical properties and allows binding to a shunt circuit.

    properties
        Shunt % Instance of a shunt circuit (Shunt class)
    end

    properties
        Y11_E % Piezoelectric Young's modulus in short circuit [Pa]
        Y11_D % Piezoelectric Young's modulus in open circuit [Pa]
        C_T   % Capacitance at constant stress [F]
        k31   % Electromechanical coupling coefficient
        h     % Thickness [m]
        b     % Width [m]
        L     % Length [m]
        rho   % Density [kg/m^3]
        C_S   % Stiffness-adjusted capacitance [F]
        E     % Young's modulus function handle (depends on time)
        A     % Cross-sectional area [m^2]
        J     % Cross-sectional moment of inertia [m^4]
    end
    
    methods
        function obj = Piezo(varargin)
            % PIEZO Constructor for the Piezo class.
            % Initializes piezoelectric properties and computes derived properties.
            %
            % Syntax:
            %   obj = Piezo('Y11_E', value, 'C_T', value, ...)
            %
            % Name-Value Pairs:
            %   'Y11_E' - Piezoelectric Young's modulus (short circuit) [Pa] (default: 62e9)
            %   'C_T'   - Capacitance at constant stress [F] (default: 7.0e-9)
            %   'k31'   - Electromechanical coupling coefficient (default: 0.351)
            %   'h'     - Thickness [m] (default: 1e-3)
            %   'b'     - Width [m] (default: 20e-3)
            %   'L'     - Length [m] (default: 22e-3)
            %   'rho'   - Density [kg/m^3] (default: 7900)
            %
            % Output:
            %   obj - Instance of the Piezo class.

            % Input parser for optional parameters
            parser = inputParser;
            addParameter(parser, 'Y11_E', 62e9);    % Default Young's modulus (62 GPa)
            addParameter(parser, 'C_T', 7.0e-9);    % Default capacitance (7 nF)
            addParameter(parser, 'k31', 0.351);     % Default electromechanical coupling coefficient
            addParameter(parser, 'h', 1e-3);        % Default thickness (1 mm)
            addParameter(parser, 'b', 20e-3);       % Default width (20 mm)
            addParameter(parser, 'L', 22e-3);       % Default length (22 mm)
            addParameter(parser, 'rho', 7900);      % Default density (7900 kg/m^3)

            % Parse inputs
            parse(parser, varargin{:});
            
            % Assign parsed values to object properties
            obj.Y11_E = parser.Results.Y11_E;
            obj.C_T = parser.Results.C_T;
            obj.k31 = parser.Results.k31;
            obj.h = parser.Results.h;
            obj.b = parser.Results.b;
            obj.L = parser.Results.L;
            obj.rho = parser.Results.rho;
            
            % Compute derived properties
            obj.C_S = obj.C_T * (1 - obj.k31^2);  % Stiffness-adjusted capacitance
            obj.Y11_D = obj.Y11_E / (1 - obj.k31^2); % Young's modulus (open circuit)
            obj.E = @(t) obj.Y11_E;  % Default Young's modulus function
            obj.A = obj.b * obj.h;    % Cross-sectional area
            obj.J = obj.b * (obj.h^3) / 12; % Cross-sectional moment of inertia
        
        end

        function obj = bindShunt(obj, shunt, w)
            % BINDSHUNT Binds a shunt circuit to the piezoelectric material.
            % This function modifies the effective Young's modulus based on
            % the electrical impedance of the shunt circuit.
            %
            % Syntax:
            %   obj = obj.bindShunt(shunt, w)
            %
            % Inputs:
            %   shunt - Shunt object (instance of the Shunt class)
            %   w     - Frequency [rad/s] at which to evaluate the shunt impedance
            %
            % Output:
            %   obj - Updated Piezo object with bound shunt.

            obj.Shunt = shunt; % Store the shunt circuit

            % Update the Young's modulus function based on the shunt impedance
            obj.E = @(t) obj.Y11_D * (1 - obj.k31^2 / (1 + 1i*w*obj.C_S*obj.Shunt.Z(w)));

        end
        
    end
end
