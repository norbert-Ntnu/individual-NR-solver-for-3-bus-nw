function [angles,magnitudes] = my_threebus_NR_xxxxx(R12,X12,R13,X13,...
    P2_spec,V2_spec,P3_spec,Q3_spec)

    angles = [0,0,0];
    magnitudes = [1,1,1];
    V2 = abs(V2_spec); V1 = 1; delta1 = 0;
    
%from impedances to admittances
z12 = complex(R12,X12);
z13 = complex(R13,X13);
y12 = 1/z12;
y13 = 1/z13;
Y11 = y12 + y13; Y12 = -y12; Y13 = -y13;
Y21 = -y12; Y22 = y12; Y23 = 0;
Y31 = -y13; Y32 = 0; Y33 = y13;

%composing the Ybus
Ybus = [Y11 Y12 Y13; Y21 Y22 Y23; Y31 Y32 Y33];
[angle_Ybus, rho_Ybus] = cart2pol(real(Ybus), imag(Ybus));

Result=[rho_Ybus ; angle_Ybus];

%admittance magnitudes and angles:

Y21_mag = Result(2,1); Y21_rad = Result(5,1);
Y22_mag = Result(2,2); Y22_rad = Result(5,2);

Y31_mag = Result(3,1); Y31_rad = Result(6,1);
Y33_mag = Result(3,3); Y33_rad = Result(6,3);
 

%the state variables column vector [delta2, delta3, |V3|]:
delta_sv = [10; 10; 10];
sv = [0; 0; 1];

%the specified powers:
spec =[P2_spec; P3_spec; Q3_spec];

%the WHILE LOOP:
iter = 0;
while max(abs(delta_sv))>= 0.0001% && iter < 10
    iter = iter + 1;

    %the LF equations giving P2_calc, P3_calc and Q3_calc:
    calc =[abs(V2) * abs(V1) * Y21_mag * cos(Y21_rad - sv(1) + delta1) + (abs(V2))^2 * Y22_mag * cos(Y22_rad); 
        abs(sv(3)) * abs(V1) * Y31_mag * cos(Y31_rad - sv(2) + delta1) + abs(sv(3))^2 * Y33_mag * cos(Y33_rad);
    -abs(sv(3)) * abs(V1) * Y31_mag * sin(Y31_rad - sv(2) + delta1) - abs(sv(3))^2 * Y33_mag * sin(Y33_rad)];
    %calc = [p2_c; p3_c; q3_c];
    
    %column vector of mismatch powers:
    d_pow = spec - calc;
    
    %the Jacobian:
    J =[abs(V2) * abs(V1) * Y21_mag * sin(Y21_rad - sv(1) + delta1), 0, 0;
        0, abs(sv(3)) * abs(V1) * Y31_mag * sin(Y31_rad - sv(2) + delta1), abs(V1) * Y31_mag * cos(Y31_rad - sv(2) + delta1) + abs(sv(3)) * 2 * Y33_mag * cos(Y33_rad);
        0, abs(sv(3)) * abs(V1) * Y31_mag * cos(Y31_rad - sv(2) + delta1), -abs(V1) * Y31_mag * sin(Y31_rad - sv(2) + delta1) - abs(sv(3)) * 2 * Y33_mag * sin(Y33_rad)];

    
    % this is where the magic happens. It is basically the equivalent of 
    % [delta_Pi, delta_Qi] = J * [delta_angles, delta_voltages] but here
    % instead of obtaining delta_x by multiplying the power vector with the
    % inverse of the Jacobian, we use smth called "matrix left division"
    delta_sv = J\d_pow; % this implies that DC = J * Dx 
    sv = sv + delta_sv;
    angles = [1,sv(1),sv(2)];
    magnitudes = [1,V2,sv(3)];
end