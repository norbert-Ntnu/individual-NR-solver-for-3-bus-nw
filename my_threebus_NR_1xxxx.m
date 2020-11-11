function [angles,magnitudes] = my_threebus_NR_1xxxx(R12,X12,R13,X13,...
    P2_spec,V2_spec,P3_spec,Q3_spec)

    angles = [0,0,0];
    magnitudes = [1,1,1];
    V2 = V2_spec; V1 = 1; delta1 = 0;
    
%from impedances to admittances
z12 = complex(R12,X12);
z13 = complex(R13,X13);
y12 = 1/z12;
y13 = 1/z13;
Y11 = y12 + y13;
Y12 = -y12;
Y13 = -y13;
Y21 = -y12;
Y22 = y12;
Y23 = 0;
Y31 = -y13;
Y32 = 0;
Y33 = y13;

%admittance magnitudes and angles:
%Y11_mag = abs(Y11); Y11_rad = atan(imag(Y11)/real(Y11));
Y21_mag = abs(Y21); Y21_rad = atan(imag(Y21)/real(Y21));
Y22_mag = abs(Y22); Y22_rad = atan(imag(Y22)/real(Y22));
Y31_mag = abs(Y31); Y31_rad = atan(imag(Y31)/real(Y31));
Y33_mag = abs(Y33); Y33_rad = atan(imag(Y33)/real(Y33));

%composing the Ybus
Ybus = [Y11 Y12 Y13; Y21 Y22 Y23; Y31 Y32 Y33];

%flat start values:
delta2 = 0; delta3 = 0; V3 = 1;

%the state variables column vector:
sv = [delta2; delta3; V3];

%flat start values:
%delta2 = 0; delta3 = 0; V3 = 1;

%the specified powers:
p2_sp = P2_spec; p3_sp = P3_spec; q3_sp = Q3_spec;

%this is where the cookie crumbles
%the WHILE LOOP:
delta_sv = [1; 1; 1];
iter = 0;
while max(abs(delta_sv))>= 0.0001% && iter < 10
    iter = iter + 1;
    
    %the state variables column vector:
    %sv = [delta2; delta3; V3];
    


    %the LF equations giving P2_calc, P3_calc and Q3_calc:
    p2_c = abs(V2)*abs(V1)*Y21_mag*cos(Y21_rad-sv(1)+delta1)+(abs(V2))^2*Y22_mag*cos(Y22_rad);
    p3_c = abs(sv(3))*abs(V1)*Y31_mag*cos(Y31_rad-sv(2)+delta1)+(abs(sv(3)))^2*Y33_mag*cos(Y33_rad);
    q3_c = -(abs(sv(3))*abs(V1)*Y31_mag*sin(Y31_rad-sv(2)+delta1))-(((abs(sv(3)))^2*sin(Y33_rad)));
    
    %column vector of mismatch powers:
    d_pow = [p2_sp-p2_c; p3_sp-p3_c; q3_sp-q3_c];
    
    %the Jacobian:
    J11 = abs(V2)*abs(V1)*Y21_mag*sin(Y21_rad-sv(1)+delta1);
    J12 = 0;
    J13 = 0;
    J21 = 0;
    J22 = abs(sv(3))*abs(V1)*Y31_mag*sin(Y31_rad-sv(2)+delta1);
    J23 = abs(V1)*Y31_mag*cos(Y31_rad-sv(2)+delta1)+abs(sv(3))*2*Y33_mag*cos(Y33_rad);
    J31 = 0;
    J32 = abs(sv(3))*abs(V1)*Y31_mag*cos(Y31_rad-sv(2)+delta1);
    J33 = -(abs(V1)*Y31_mag*sin(Y31_rad-sv(2)+delta1))-(abs(sv(3))*2*sin(Y33_rad));

    
    % this is the Jacobian
    J = [J11 J12 J13; J21 J22 J23; J31 J32 J33];
    
    % this is where the magic happens. It is basically the equivalent of 
    % [delta_Pi, delta_Qi] = J * [delta_angles, delta_voltages] but here
    % instead of obtaining delta_x by multiplying the power vector with the
    % inverse of the Jacobian, we use smth called "matrix left division"
    delta_sv = J\d_pow; % this implies that DC = J * Dx 
    sv = sv+delta_sv;
    angles = [1,sv(1),sv(2)];
    magnitudes = [1,V2,sv(3)];
end