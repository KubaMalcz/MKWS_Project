wsp_przewodzenia_ciepla=230; %podane w W/(M*K) uaga jest zmienne, trzeba bedzie to zaimplementowac
wsp_przejmowania_ciepla=4000; %podane w W/(M^2*K)
grubosc_scianki=0.003; %w m
ilosc_podzialow=10;
ilosc_krokow_czasowych=8000;
czas=3; %w sekundach
T=ones(ilosc_podzialow,ilosc_krokow_czasowych);
Temperatura_poczatkowa=293; %w Kelvinach
T(:,1)=Temperatura_poczatkowa;
gestosc=2700; %kg/m^3
cieplo_wlasciwe=934; % J/(kg*K) uwaga zmienne z temperatura
dyfuzyjnosc_termiczna=wsp_przewodzenia_ciepla/(cieplo_wlasciwe*gestosc); %w m^2/s
Temperatura_w_komorze_spalania=3000; %w Kelvinach

M=(grubosc_scianki/ilosc_podzialow)^2/(dyfuzyjnosc_termiczna*(czas/ilosc_krokow_czasowych));
C=(wsp_przejmowania_ciepla*(grubosc_scianki/ilosc_podzialow))/wsp_przewodzenia_ciepla;

for k=2:ilosc_krokow_czasowych
    T(1,k)=1/M*(2*C*Temperatura_w_komorze_spalania+(M-2*(C+1))*T(1,k-1)+2*T(2,k-1));
    for i=2:ilosc_podzialow-1
        T(i,k)=(T(i-1,k-1)+(M-2)*T(i,k-1)+T(i+1,k-1))/M; %na razie wzor 3.198
    end
    T(ilosc_podzialow,k)=T(ilosc_podzialow-1,k);
end

%Temperatura obliczana jest na brzegach podzia³ów
disp(T)
czas=[czas/ilosc_krokow_czasowych:czas/ilosc_krokow_czasowych:czas];
x=T(5,:);
plot(czas,x)