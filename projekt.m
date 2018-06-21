clf
%Input data and parameters
number_of_elements=10;
timesteps=62000;
material=2; % 1-Aluminium 2024; 2-Alloy Steel; 3-Stainless Steel
bakelit_ON=0;
chamber_wall_thickness=0.003; %w m
insulation_thickness=0.001; %w m
Di=0.063; %chamber inside diameter (meters)
massflow_n2o=0.013448728016102 ; %w kg/s
massflow_kerosene=0.002359425967737 ; %w kg/s
L=0.24; %chamber length (meters)
time=30; %in seconds
initial_temperature=300; %Kelvin
temperature_in_combastion_chamber=3000; %Kelvin
%cooling
water_temperature=300; %w K
waterflow_speed=0.1; %m/s (if no cooling put 0)
equivalent_radius_of_cooling_channel=0.015; %w metrach

%DO NOT CHANGE ANYTHING BELOW!

if bakelit_ON==1
    grubosc_calkowita=chamber_wall_thickness+insulation_thickness;
else
    grubosc_calkowita=chamber_wall_thickness;
end

%wsp_przejmowania_ciepla=1800; %podane w W/(M^2*K) dobrane z dupy, trzeba to policzyæ
Cp_N2O=0.21;  %w Btu/(lbm*R) www.engineeringtoolbox.com
Cp_kerosene=0.48; %w Btu/(lbm*R) www.engineeringtoolbox.com
Di2=Di*39.37; % inches
S=pi()/4*Di2^2*0.0069444; %w ft^2
massflow_n2o2=massflow_n2o*7936.64; %w lb/hr
massflow_kerosene2=massflow_kerosene*7936.64; %w lb/hr
Cp_average=(Cp_N2O*massflow_n2o+Cp_kerosene*massflow_kerosene)/(massflow_n2o+massflow_kerosene); %w Btu/(lbm*R)
G=(massflow_kerosene2+massflow_n2o2)/S; %w lb/(hr*ft^2)
L2=L*39.37; %inches
wsp_przejmowania_ciepla=0.024*Cp_average*G^0.8/Di2^0.2*(1+(Di2/L2)^0.7)*1269/223; %w W/(m^2*K) Mark's Hdbk. for Mechanical Engineers

if material==1 %Aluminium 2024
    case_density=2770;%kg/m^3
end
if material==2 %Alloy Steel 4130
    case_density=7840;%kg/m^3
end
if material==3 %Stainless Steel
    case_density=7900;%kg/m^3
end


h=chamber_wall_thickness/number_of_elements; %dimensional step
dt=time/timesteps; %time step
h_bakelit=insulation_thickness/1.5;
T=ones(number_of_elements+1,timesteps);
T(:,1)=initial_temperature;
T_iz=ones(3,timesteps);
T_iz(:,1)=initial_temperature;

%Obudowa komory
p_sh=specific_heat(material);
p_c=conductivity(material);
wsp_przewodzenia_ciepla=polyval(p_c,T(1,1)); %W/(M*K) vary with temperature
cieplo_wlasciwe=polyval(p_sh,T(1,1)); % J/(kg*K) vary with temperature
dyfuzyjnosc_termiczna=wsp_przewodzenia_ciepla/(cieplo_wlasciwe*case_density); %m^2/s
M=h^2/(dyfuzyjnosc_termiczna*dt);
C=(wsp_przejmowania_ciepla*h)/wsp_przewodzenia_ciepla;
if M<2
    error('Too few timesteps')
end

%Izolacja bakielitowa
wsp_przewodzenia_ciepla_bakelit=0.2; %W/(M*K) constant
gestosc_bakelit=1340; %kg/m^3
cieplo_wlasciwe_bakelit=920; % J/(kg*K) constant
dyfuzyjnosc_termiczna_bakelit=wsp_przewodzenia_ciepla_bakelit/(cieplo_wlasciwe_bakelit*gestosc_bakelit); %w m^2/s
M_bakelit=h_bakelit^2/(dyfuzyjnosc_termiczna_bakelit*dt);
C_bakelit=(wsp_przejmowania_ciepla*h_bakelit)/wsp_przewodzenia_ciepla_bakelit;
if M_bakelit<2
    error('Too few timesteps (bakelit)')
end

%cooling parameters
lepkosc_kinematyczna=1*10^(-6); % m^2/s for 293K
Re_lp=waterflow_speed*equivalent_radius_of_cooling_channel/lepkosc_kinematyczna;
Pr_p=7; %for 293K
Pr_s=2; %for 363K
Nu_p=0.037*(Re_lp^0.8)*(Pr_p^0.43)*((Pr_p/Pr_s)^0.25); % 5.170
wsp_przewodzenia_ciepla_wody=0.65; % W/(m*K) for 323K
wsp_przejmowania_ciepla_cool=wsp_przewodzenia_ciepla_wody*Nu_p/equivalent_radius_of_cooling_channel; %podane w W/(M^2*K)

if bakelit_ON==1
    korekcja=3;
else
    korekcja=0;
end

wsp_przewodzenia_ciepla=polyval(p_c,T(1,1)); %W/(M*K) vary with temperature
cieplo_wlasciwe=polyval(p_sh,T(1,1)); % J/(kg*K) ary with temperature
dyfuzyjnosc_termiczna=wsp_przewodzenia_ciepla/(cieplo_wlasciwe*case_density); % m^2/s

for k=2:timesteps+1
    if bakelit_ON==1
        T_iz(1,k)=1/M_bakelit*(2*C_bakelit*temperature_in_combastion_chamber+(M_bakelit-2*(C_bakelit+1))*T_iz(1,k-1)+2*T_iz(2,k-1));
        if M_bakelit<2*(C_bakelit+1)
        error('Too few timesteps (bakelit)')
        end
    else
        T(1,k)=1/M*(2*C*temperature_in_combastion_chamber+(M-2*(C+1))*T(1,k-1)+2*T(2,k-1));
        if M<2*(C+1)
        error('Too few timesteps (bakelit)')
        end
    end
    
    for i=2:(number_of_elements+korekcja)
        
        if (bakelit_ON==1) && (i==2) 
            T_iz(i,k)=T_iz(i-1,k-1)/M_bakelit+(M_bakelit-2)*T_iz(i,k-1)/M_bakelit+T_iz(i+1,k-1)/M_bakelit; %na razie wzor 3.198
        elseif (bakelit_ON==1) && (i==3)
            p=wsp_przewodzenia_ciepla*h_bakelit/(wsp_przewodzenia_ciepla_bakelit*h);
            T_iz(i,k)=1/(p+1)*(T_iz(i-1,k)+p*T(1,k-1)); %idealny styk temperatura Ts na podstawie wzoru 3.223
        elseif (bakelit_ON==1) && (i==4)
            T((i-korekcja),k)=(T_iz(i-1,k-1)+(M-2)*T(i-korekcja,k-1)+T(i+1-korekcja,k-1))/M;
        else
            %Obudowa komory
            wsp_przewodzenia_ciepla=polyval(p_c,T(i-korekcja,k-1)); %podane w W/(M*K) uwaga jest zmienne
            cieplo_wlasciwe=polyval(p_sh,T(i-korekcja,k-1)); % J/(kg*K) uwaga zmienne z temperatura
            dyfuzyjnosc_termiczna=wsp_przewodzenia_ciepla/(cieplo_wlasciwe*case_density); %w m^2/s
            M=h^2/(dyfuzyjnosc_termiczna*dt);
            C=(wsp_przejmowania_ciepla*h)/wsp_przewodzenia_ciepla;
            if M<2
                error('Too few timesteps')
            end 
            T((i-korekcja),k)=(T(i-1-korekcja,k-1)+(M-2)*T(i-korekcja,k-1)+T(i+1-korekcja,k-1))/M;
        end
    end
    
    %T(number_of_elements,k)=T(number_of_elements-1,k);
    Ccool=(wsp_przejmowania_ciepla_cool*h)/wsp_przewodzenia_ciepla;
    T((number_of_elements+1),k)=1/M*(2*Ccool*water_temperature+(M-2*(Ccool+1))*T(number_of_elements+1,k-1)+2*T(number_of_elements,k-1));
    
end
%Temperatura obliczana jest na brzegach podzia³ów

time_m=dt:dt:time;
wymiar=0:h:chamber_wall_thickness;
hold on
xlabel('X dimension [meters]')
ylabel('Temperature [Kelvins]')
for j=100:timesteps/10:timesteps
    x=T(1:end,j);
    plot(wymiar,x)
end
hold off

dlmwrite('Full_Temperature_Data.csv',T);

% Function matching value of material specific heat to its temperature
function sf=specific_heat(material)
if material==1 %Aluminium 2024
    sf_table=[875 925 1042]; %w J/(kg*K)
    temp_table=[300 400 600]; %w K
    T=table(temp_table, sf_table);
    sf=polyfit(T.temp_table, T.sf_table,1);
end
if material==2 %Alloy Steel 4130
    sf_table=[430 490 570 670 1260]; %w J/(kg*K)
    temp_table=[300 400 600 800 1000]; %w K
    T=table(temp_table, sf_table);
    sf=polyfit(T.temp_table, T.sf_table,1);
end
if material==3 %Stainless Steel
    sf_table=[477 515 557 582 611 640]; %w J/(kg*K)
    temp_table=[300 400 600 800 1000 1200]; %w K
    T=table(temp_table, sf_table);
    sf=polyfit(T.temp_table, T.sf_table,1);
end
end

% Function matching value of material conductivity to its temperature
function con=conductivity(material)
if material==1 %Aluminium 2024
    con_table=[177 186 186]; %w W/(m*K)
    temp_table=[300 400 600]; %w K
    T=table(temp_table, con_table);
    con=polyfit(T.temp_table, T.con_table,1);
end
if material==2 %Alloy Steel 4130
    con_table=[41.5 43.1 40.3 35 28.5 27.7]; %%w W/(m*K)
    temp_table=[300 400 600 800 1000 1200]; %w K
    T=table(temp_table, con_table);
    con=polyfit(T.temp_table, T.con_table,1);
end
if material==3 %Stainless Steel
    con_table=[14.9 16.6 19.8 22.6 25.4 28]; %%w W/(m*K)
    temp_table=[300 400 600 800 1000 1200]; %w K
    T=table(temp_table, con_table);
    con=polyfit(T.temp_table, T.con_table,1);
end
end