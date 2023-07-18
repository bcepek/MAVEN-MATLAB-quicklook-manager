function [amu, anu, c, q, m_p, m_O, m_O2, m_CO2, Rm] = constants
%---BEGIN---Constants
amu = 1.660539 * 10^-24; % [g]
anu = 1;
c = 3 * 10^10; %#ok<NASGU> %[cm/s]
q = 4.8 * 10^(-10); % SGSE charge unit
% amu = 1
m_p = 1 * amu;
m_O = 16 * amu;
m_O2 = 32 * amu;
m_CO2 = 44 * amu;
Rm = 3390;
%---END---Constants
end