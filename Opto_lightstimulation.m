%optical stimulation
%3 x 30 repetitions of 100ms every 5s
clear all

a = arduino();
 configurePin(a,'D11', 'unset'); %ttl pulse
 configurePin(a,'D10', 'unset'); %light stimulation
for k = 1:3;
for j = 1:3; 
for i = 1:30; 
writePWMVoltage(a,'D10',2.5);
writePWMVoltage(a,'D11',3);
writePWMVoltage(a,'D11',0);
pause(0.1);
writePWMVoltage(a,'D10',0);
pause(0.2);
end
pause(5)
end
pause(30)
end