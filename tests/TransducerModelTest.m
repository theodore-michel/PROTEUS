classdef (SharedTestFixtures = {ProteusTestEnvironment}) ...
        TransducerModelTest < matlab.unittest.TestCase
%TRANSDUCERMODELTEST Calibrated transducer models (Transducer-models/*.mat)
%and the transducer catalogue.

    properties
        Root    % Repository root
    end

    methods (TestMethodSetup)
        function getRoot(testCase)
            env = testCase.getSharedTestFixtures('ProteusTestEnvironment');
            testCase.Root = env.RepoRoot;
        end
    end

    methods (Test, TestTags = {'Fast'})
        function p41KeepsCatalogueValuesWithCalibratedResponse(testCase)
            T.Type = 'P4-1';
            T = get_transducer_properties(T);
            % Original catalogue values (reference code):
            testCase.verifyEqual([T.BandwidthLow T.BandwidthHigh ...
                T.CenterFrequency], [1.5e6 3.5e6 2.5e6]);
            testCase.verifyEqual([T.NumberOfElements T.Pitch ...
                T.ElementWidth T.ElementHeight T.ElevationFocus], ...
                [96 0.295e-3 0.245e-3 16e-3 110e-3]);
            testCase.verifyEqual(T.Radius, 0);
            % Calibrated one-way impulse response (Verasonics format):
            M = load(fullfile(testCase.Root,'Transducer-models','P4-1.mat'));
            testCase.verifyEqual(T.ImpulseResponseType, 'Calibrated');
            testCase.verifyEqual(T.calibrated, 1);
            testCase.verifyEqual(T.SamplingRate, 250e6);
            testCase.verifyEqual(T.TransmitImpulseResponse, M.Trans.IR1wy(:)'*250e6);
            testCase.verifyEqual(T.ReceiveImpulseResponse,  M.Trans.IR1wy(:)'*250e6);
            % Unit gain, as the estimated impulse response:
            gain = @(T) max(abs(fft(T.TransmitImpulseResponse,2^16)))/T.SamplingRate;
            testCase.verifyEqual(gain(T), 1, 'RelTol', 0.01);
            testCase.verifyEqual(gain(reference_transducer(T)), 1, 'RelTol', 0.01);
        end

        function c16dGeometryComesFromModelFile(testCase)
            T.Type = 'C1-6D';
            T = get_transducer_properties(T);
            M = load(fullfile(testCase.Root,'Transducer-models','C1-6D.mat'));
            testCase.verifyEqual(T.Configuration, 'Curved');
            testCase.verifyEqual([T.NumberOfElements T.Pitch T.ElementWidth ...
                T.ElementHeight T.ElevationFocus T.Radius], ...
                [M.Transducer.NumberOfElements M.Transducer.Pitch ...
                M.Transducer.ElementWidth M.Transducer.ElementHeight ...
                M.Transducer.ElevationFocus M.Transducer.Radius]);
            testCase.verifyEqual(T.TransmitImpulseResponse, M.Transducer.IR_trans);
            testCase.verifyEqual(T.ReceiveImpulseResponse,  M.Transducer.IR_rec);
            testCase.verifyEqual(T.SamplingRate, M.Transducer.Fs_trans);
            % Grid-dependent data in the model file is not used:
            testCase.verifyFalse(isfield(T,'integration_points'));
        end

        function modelIsFoundFromAnyFolder(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            import matlab.unittest.fixtures.CurrentFolderFixture
            tmp = testCase.applyFixture(TemporaryFolderFixture);
            testCase.applyFixture(CurrentFolderFixture(tmp.Folder));
            T.Type = 'C1-6D';
            T = get_transducer_properties(T);
            testCase.verifyEqual(T.ImpulseResponseType, 'Calibrated');
        end

        function catalogueWithoutModelIsEstimated(testCase)
            for type = {'L22-14v','9L-D','RC6gV'}
                T.Type = type{1};
                T = get_transducer_properties(T);
                testCase.verifyEqual(T.ImpulseResponseType, 'Estimate', type{1});
                testCase.verifyEqual(T.calibrated, 0, type{1});
                testCase.verifyEqual(T.Radius, 0, type{1});
                clear T
            end
        end

        function guiDefaultIsCalibratedP41(testCase)
            T = reset_transducer();
            testCase.verifyEqual(T.Type, 'P4-1');
            testCase.verifyEqual(T.ImpulseResponseType, 'Calibrated');
            Tx = reset_transmit(T);
            testCase.verifyEqual(Tx.CenterFrequency, 2.5e6);
            testCase.verifyTrue(all(isfinite(Tx.PressureSignal)));
            testCase.verifyEqual(min(Tx.PressureSignal), -Tx.AcousticPressure, ...
                'RelTol', 1e-12);
        end

        function proteusFormatWithOwnImpulseResponses(testCase)
            folder = testCase.modelFolder();
            Transducer = struct('Pitch', 0.4e-3, ...
                'TransmitImpulseResponse', [0 1 -1 0]', ...
                'ReceiveImpulseResponse', [0 2 -2 0], ...
                'SamplingRate', 100e6, 'integration_points', zeros(2,2,3));
            save(fullfile(folder,'L22-14v.mat'),'Transducer');
            T.Type = 'L22-14v';
            T = get_transducer_properties(T);
            T = load_transducer_model(T, folder);
            testCase.verifyEqual(T.Pitch, 0.4e-3);
            testCase.verifyEqual(T.NumberOfElements, 128);   % catalogue
            testCase.verifyEqual(T.TransmitImpulseResponse, [0 1 -1 0]);
            testCase.verifyEqual(T.ReceiveImpulseResponse,  [0 2 -2 0]);
            testCase.verifyEqual(T.SamplingRate, 100e6);
            testCase.verifyFalse(isfield(T,'integration_points'));
        end

        function customLinearDropsCurvedConfiguration(testCase)
            T.Type = 'C1-6D';
            T = get_transducer_properties(T);
            T.Type = 'Custom linear';
            T = get_transducer_properties(T);
            testCase.verifyEqual(T.Configuration, 'Linear');
            testCase.verifyEqual(T.Radius, 0);
            testCase.verifyEqual(T.calibrated, 0);
            testCase.verifyEqual(T.ImpulseResponseType, 'Load file');
            testCase.verifyEqual(T.Pitch, 0.35884e-3, 'RelTol', 1e-4);
        end

        function invalidModelFilesError(testCase)
            folder = testCase.modelFolder();
            T.Type = 'X';
            Transducer = struct('IR_trans',1,'IR_rec',1,'Fs_trans',1,'Fs_rec',2);
            save(fullfile(folder,'X.mat'),'Transducer');
            testCase.verifyError(@() load_transducer_model(T, folder), ...
                'load_transducer_model:samplingRate');
            Transducer = struct('Pitch',1);
            save(fullfile(folder,'X.mat'),'Transducer');
            testCase.verifyError(@() load_transducer_model(T, folder), ...
                'load_transducer_model:noImpulseResponse');
            Other = 1; %#ok<NASGU>
            save(fullfile(folder,'X.mat'),'Other');
            testCase.verifyError(@() load_transducer_model(T, folder), ...
                'load_transducer_model:unknownFormat');
        end
    end

    methods (Access = private)
        function folder = modelFolder(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            tmp = testCase.applyFixture(TemporaryFolderFixture);
            folder = tmp.Folder;
        end
    end
end
