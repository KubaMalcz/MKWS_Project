clf
wsp_przewodzenia_ciepla=230; %podane w W/(M*K) uwaga jest zmienne, trzeba bedzie to zaimplementowac
wsp_przejmowania_ciepla=400; %podane w W/(M^2*K)
grubosc_scianki=0.003; %w m
ilosc_podzialow=12;
ilosc_krokow_czasowych=10000;
czas=1; %w sekundach
T=ones(ilosc_podzialow,ilosc_krokow_czasowych);
Temperatura_poczatkowa=293; %w Kelvinach
T(:,1)=Temperatura_poczatkowa;
gestosc=2700; %kg/m^3
cieplo_wlasciwe=934; % J/(kg*K) uwaga zmienne z temperatura
dyfuzyjnosc_termiczna=wsp_przewodzenia_ciepla/(cieplo_wlasciwe*gestosc); %w m^2/s
Temperatura_w_komorze_spalania=3000; %w Kelvinach

M=(grubosc_scianki/ilosc_podzialow)^2/(dyfuzyjnosc_termiczna*(czas/ilosc_krokow_czasowych));
C=(wsp_przejmowania_ciepla*(grubosc_scianki/ilosc_podzialow))/wsp_przewodzenia_ciepla;

temperatura_wody=293; %w K
lepkosc_kinematyczna=1*10^(-6); % w m^2/s dla 293K
szybkosc_przeplywu=1; %w m/s
srednica_rownowazna_kanalu_chlodzacego=0.015 %w metrach
Re_lp=szybkosc_przeplywu*srednica_rownowazna_kanalu_chlodzacego/lepkosc_kinematyczna
Pr_p=7; %dla 293K
Pr_s=2; %dla 363K
Nu_p=0.037*(Re_lp^0.8)*(Pr_p^0.43)*((Pr_p/Pr_s)^0.25) %wzor 5.170
wsp_przewodzenia_ciepla_wody=0.65; %w W/(m*K) dla 323K
dyfuzyjnosc_termiczna_cool=wsp_przewodzenia_ciepla/(cieplo_wlasciwe*gestosc); %w m^2/s
wsp_przejmowania_ciepla_cool=wsp_przewodzenia_ciepla_wody*Nu_p/srednica_rownowazna_kanalu_chlodzacego %podane w W/(M^2*K)
Mcool=(grubosc_scianki/ilosc_podzialow)^2/(dyfuzyjnosc_termiczna_cool*(czas/ilosc_krokow_czasowych));
Ccool=(wsp_przejmowania_ciepla_cool*(grubosc_scianki/ilosc_podzialow))/wsp_przewodzenia_ciepla;


for k=2:ilosc_krokow_czasowych
    T(1,k)=1/M*(2*C*Temperatura_w_komorze_spalania+(M-2*(C+1))*T(1,k-1)+2*T(2,k-1));
    for i=2:ilosc_podzialow-1
        T(i,k)=(T(i-1,k-1)+(M-2)*T(i,k-1)+T(i+1,k-1))/M; %na razie wzor 3.198
    end
    %T(ilosc_podzialow,k)=T(ilosc_podzialow-1,k);
    T(ilosc_podzialow,k)=1/Mcool*(2*Ccool*temperatura_wody+(Mcool-2*(Ccool+1))*T(ilosc_podzialow,k-1)+2*T(ilosc_podzialow-1,k-1));
end

%Temperatura obliczana jest na brzegach podzia³ów
czas_m=czas/ilosc_krokow_czasowych:czas/ilosc_krokow_czasowych:czas;
wymiar=0:grubosc_scianki/(ilosc_podzialow-1):grubosc_scianki;
hold on
for j=100:ilosc_krokow_czasowych/10:ilosc_krokow_czasowych
    x=T(:,j);
    plot(wymiar,x)
end
hold off
