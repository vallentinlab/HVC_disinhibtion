%opening and plotting .rhd files generated using intan RHD2000 software
%with Neuronexus probes

read_Intan_RHD2000_file
subplot(2,1,1);
plot(t_amplifier, amplifier_data(1,:));
subplot(2,1,2);
plot(t_board_adc, board_adc_data(1,:));
