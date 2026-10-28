classdef (SharedTestFixtures = {ProteusTestEnvironment}) ...
        TransmitSignalTest < matlab.unittest.TestCase
%TRANSMITSIGNALTEST Custom transmit signals (issue #3: a mostly negative
%custom voltage waveform gave a zero amplitude and NaN signals).

    methods (Test, TestTags = {'Fast'})
        function negativeCustomVoltageIsScaledToAmplitude(testCase)
            Transducer = reference_transducer(reset_transducer());
            Transmit = reset_transmit(Transducer);
            Transmit.Type = 'Custom (voltage)';
            Transmit.AmplitudeMode = 'Voltage';
            Transmit.VoltageSignal = -[0 0.5 1 0.5 0 0.2 0];
            Transmit.VoltageAmplitude = 10;      % [V]

            Transmit = get_voltage_signal(Transmit);
            Transmit = get_pressure_signal(Transmit, Transducer);

            V = Transmit.VoltageSignal;
            testCase.verifyTrue(all(isfinite(V)));
            testCase.verifyEqual(max(abs(V)), 10, 'RelTol', 1e-12);
            testCase.verifyEqual(V(3), -10, 'RelTol', 1e-12);   % shape kept
            testCase.verifyTrue(isfinite(Transmit.MechanicalIndex));
            testCase.verifyGreaterThan(Transmit.AcousticPressure, 0);
        end

        function exampleCustomFilesHoldVariables(testCase)
            % The custom-data files are mat files with plain variables (as
            % documented in documentation/TransmitGUI.md, issue #8).
            env = testCase.getSharedTestFixtures('ProteusTestEnvironment');
            folder = fullfile(env.RepoRoot, 'example_custom_data');
            expected = {'custom_voltage_signal.mat','V'; ...
                'custom_pressure_signal.mat','p'; ...
                'custom_delays.mat','delays'; ...
                'custom_apodization.mat','apodization'; ...
                'custom_transmit_impulse_response.mat','IR'; ...
                'custom_receive_impulse_response.mat','IR'};
            for k = 1:size(expected,1)
                w = whos('-file', fullfile(folder, expected{k,1}));
                testCase.verifyEqual({w.name}, expected(k,2), expected{k,1});
                testCase.verifyTrue(strcmp(w.class,'double'), expected{k,1});
            end
        end
    end
end
